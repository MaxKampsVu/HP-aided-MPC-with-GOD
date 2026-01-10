#include "dmgod_online_evaluator.h"

namespace dmAsyncAsteriskGOD {

    OnlineEvaluator::OnlineEvaluator(int nP, int id, int security_param, std::shared_ptr<NetIOMP> network, PreprocCircuit<Field> preproc, 
        LevelOrderedCircuit circ, int threads, uint64_t seed) 
        : nP_(nP), id_(id), security_param_(security_param), rgen_(id, nP, seed), network_(std::move(network)), 
        preproc_(std::move(preproc)), circ_(std::move(circ)), wires_(circ.num_gates), q_val_(circ.num_gates), 
        q_sh_(circ.num_gates), start_recv_(false)
    {
        tpool_ = std::make_shared<ThreadPool>(threads);

        if (id_ == 0) {            
            static ZZ_pContext ZZ_p_ctx;
            ZZ_p_ctx.save();
            for(size_t pid = 1; pid <= nP_; pid++) {
                tpool_->enqueue([&,pid]() {                    
                    ZZ_p_ctx.restore();
                    {
                        std::unique_lock<std::mutex> lock(mtx_);
                        cv_start_recv_.wait(lock, [&]() { return start_recv_; });
                    }
                    for (size_t depth = 0; depth < circ_.gates_by_level.size(); depth++) {
                        size_t total_comm;
                        network_->recv(pid, &total_comm, sizeof(size_t));
                        network_->getRecvChannel(pid)->flush();
                        std::vector<Field> online_comm_to_HP(total_comm);
                        network_->recv(pid, online_comm_to_HP.data(), sizeof(Field) * total_comm);
                        network_->getRecvChannel(pid)->flush();
                        {
                            std::lock_guard<std::mutex> lock(mtx_);
                            message_buffer_[depth].push({pid, online_comm_to_HP});
                        }
                        cv_.notify_one();
                    }            
                });
            }
        }
    }

    OnlineEvaluator::~OnlineEvaluator() {
        tpool_.reset();
    }

    void OnlineEvaluator::setRandomInputs() {
        for (auto &g : circ_.gates_by_level[0]) {
            if (g->type == GateType::kInp) {
                //randomizeZZp(rgen_.all(), wires_[g->out], sizeof(Field));
                wires_[g->out] = 0;
            }
        }
    }

    void OnlineEvaluator::setInputs(const std::unordered_map<wire_t, Field> &inputs) {
        std::vector<Field> masked_values;
        for (auto &g : circ_.gates_by_level[0]) {
            if (g->type == GateType::kInp) {
                // Find the input provider pid for the current input gate
                auto *pre_input = static_cast<PreprocInput<Field> *>(preproc_.gates[g->out].get());
                auto pid = pre_input->pid;

                // HP waits for masked value from input provider 
                if (id_ == 0) {
                    network_->recv(pid, &q_val_[g->out], sizeof(Field));
                    network_->getRecvChannel(pid)->flush();
                    std::vector<std::future<void>> send_t;
                    for (size_t i = 1; i <= nP_; i++) {
                        if (i == pid)
                            continue;
                        // Forwards masked input to all parties except dealer 
                        send_t.push_back(tpool_->enqueue([&,i]() {
                            network_->send(i, &q_val_[g->out], sizeof(Field));
                            network_->getSendChannel(i)->flush();
                        }));
                    }
                    for(auto& t : send_t) {
                        if (t.valid()) {
                            t.wait();
                        }
                    }
                }
                // Receive mask to input gate from HP
                else {
                    if (pid == id_) {
                        q_val_[g->out] = pre_input->mask_value + inputs.at(g->out);
                        network_->send(0, &q_val_[g->out], sizeof(Field));
                        network_->getSendChannel(0)->flush();
                    }
                    else {
                        network_->recv(0, &q_val_[g->out], sizeof(Field));
                        network_->getRecvChannel(0)->flush();
                    }
                }
                
                wires_[g->out] = q_val_[g->out];
            }
        }
    }

    void OnlineEvaluator::evaluateGatesAtDepthPartySend(size_t depth, std::vector<Field> &mult_nonTP, std::vector<Field> &mac_components) {
        for (auto &gate : circ_.gates_by_level[depth]) {
            switch (gate->type) {
                case GateType::kMul: {
                    auto *g = static_cast<FIn2Gate*>(gate.get());
                    auto *pre_out = static_cast<PreprocMultGate<Field>*>(preproc_.gates[g->out].get());
                    q_val_[g->out] = 0;

                    auto &m_in1 = preproc_.gates[g->in1]->mask;
                    auto &m_in2 = preproc_.gates[g->in2]->mask;
                    auto &m_prod = pre_out->mask_prod;
                    auto &m_out = pre_out->mask;
                    
                    auto q_share = m_prod + m_out - m_in1 * wires_[g->in2] - m_in2 * wires_[g->in1];   
                    q_share.add(wires_[g->in1] * wires_[g->in2], id_);                    
                    
                    q_sh_[g->out] = q_share;
                    mult_nonTP.push_back(q_share.getValue());
                    mac_components.push_back(q_share.getMACComponent());
                    break;
                }

                case GateType::kDotprod: {
                    auto *g = static_cast<SIMDGate*>(gate.get());
                    auto *pre_out = static_cast<PreprocDotpGate<Field>*>(preproc_.gates[g->out].get());

                    auto &m_vec_prod = pre_out->mask_prod;
                    auto &m_out = pre_out->mask;
                    
                    auto q_share = m_vec_prod + m_out;
                    for(size_t i = 0; i < g->in1.size(); i++) {
                        const auto& m_ai = preproc_.gates[g->in1[i]]->mask;
                        const auto& m_bi = preproc_.gates[g->in2[i]]->mask;

                        q_share = q_share - m_ai * wires_[g->in2[i]] - m_bi * wires_[g->in1[i]];
                        q_share.add(wires_[g->in1[i]] * wires_[g->in2[i]], id_);                 
                    }                   
                    
                    q_sh_[g->out] = q_share;
                    mult_nonTP.push_back(q_share.getValue());
                    mac_components.push_back(q_share.getMACComponent());
                    break;
                }

                case GateType::kAdd:
                case GateType::kSub:
                case GateType::kConstAdd:
                case GateType::kConstMul:
                    break;

                default:
                    break;
            }
        }
    }


    void OnlineEvaluator::evaluateGatesAtDepthPartyRecv(size_t depth, std::vector<Field> mult_all) {
        size_t idx_mult = 0;

        for (auto &gate : circ_.gates_by_level[depth]) {
            switch (gate->type) {
                case GateType::kAdd: {
                    auto *g = static_cast<FIn2Gate *>(gate.get());
                    wires_[g->out] = wires_[g->in1] + wires_[g->in2];                        
                    q_val_[g->out] = 0;
                    break;
                }

                case GateType::kSub: {
                    auto *g = static_cast<FIn2Gate *>(gate.get());
                    wires_[g->out] = wires_[g->in1] - wires_[g->in2];
                    q_val_[g->out] = 0;
                    break;
                }

                case GateType::kConstAdd: {
                    auto *g = static_cast<ConstOpGate<Field> *>(gate.get());
                    wires_[g->out] = wires_[g->in] + g->cval;
                    break;
                }

                case GateType::kConstMul: {
                    auto *g = static_cast<ConstOpGate<Field> *>(gate.get());
                    wires_[g->out] = wires_[g->in] * g->cval;
                    break;
                }

                case GateType::kMul: {
                    auto *g = static_cast<FIn2Gate *>(gate.get());
                    q_val_[g->out] = mult_all[idx_mult];
                    wires_[g->out] = q_val_[g->out];
                    idx_mult++;
                    break;
                }

                case GateType::kDotprod: {
                    auto *g = static_cast<SIMDGate*>(gate.get());
                    q_val_[g->out] = mult_all[idx_mult];
                    wires_[g->out] = q_val_[g->out];
                    idx_mult++;
                    break;
                }
                
                default:
                    break;
            }
        }
    }

    void OnlineEvaluator::evaluateGatesAtDepth(size_t depth) {
        size_t mult_num = 0;
        std::vector<Field> mult_nonTP;
        std::vector<Field> mac_components;

        // For gates with constant do nothing 
        for (auto &gate : circ_.gates_by_level[depth]) {
            switch (gate->type) {
                case GateType::kInp: 
                case GateType::kAdd:
                case GateType::kSub:
                case GateType::kConstAdd:
                case GateType::kConstMul:
                    break;

                case GateType::kMul: {
                    mult_num++;
                    break;
                }

                case GateType::kDotprod: {
                    mult_num++;
                    break;
                }

                default:
                    break;
            }
        }

        // Number of multiplication gates at the layer 
        size_t total_comm_to_HP = 0;
        if(mult_num != 0) {
            total_comm_to_HP = mult_num + 4;
        }

        size_t total_comm = mult_num;
        std::vector<Field> mult_all(mult_num);        
        std::vector<Field> agg_values(total_comm);

        evaluateGatesAtDepthPartySend(depth, mult_nonTP, mac_components);
        if (id_ != 0) {            
            auto total_comm_plus_dig = total_comm;
            fieldDig layer_dig;

            if(mult_num != 0) {
                total_comm_plus_dig += 4;
                layer_dig = hashFields(mac_components);
            }

            std::vector<Field> online_comm_to_HP(total_comm_plus_dig);

            for (size_t i = 0; i < mult_num; i++) {
                online_comm_to_HP[i] = mult_nonTP[i];
            }


            for (size_t i = mult_num; i < total_comm_plus_dig; i++) {
                online_comm_to_HP[i] = layer_dig[i - mult_num];
            }

            network_->send(0, &total_comm_plus_dig, sizeof(size_t));
            network_->getSendChannel(0)->flush();
            network_->send(0, online_comm_to_HP.data(), sizeof(Field) * total_comm_plus_dig, true);
            network_->getSendChannel(0)->flush();
        }
        else {
            for (size_t i=0; i<mult_num; i++) {
                agg_values[i] = mult_nonTP[i];
            }

            std::vector<std::future<void>> send_t;    
            for(size_t pid = 1; pid <= nP_; pid++) {
                send_t.push_back(tpool_->enqueue([&,pid]() {
                    network_->send(pid, agg_values.data(), sizeof(Field) * total_comm);
                    network_->getSendChannel(pid)->flush();
                }));
            }
            for(auto& t : send_t) {
                if (t.valid()) {
                    t.wait();
                }
            }
        }
        if (id_ != 0) {
            network_->recv(0, agg_values.data(), sizeof(Field) * total_comm);
            network_->getRecvChannel(0)->flush();
            for(size_t i = 0; i < mult_num; i++) {
                mult_all[i] = agg_values[i] + mult_nonTP[i];
            }
        }
        else {
            {
                std::lock_guard<std::mutex> lock(mtx_);
                start_recv_ = true;
            }
            cv_start_recv_.notify_all();
            {
                std::unique_lock<std::mutex> lock(mtx_);
                cv_.wait(lock, [&]() { return message_buffer_[depth].size() >= 1; });
            }

            Message message = message_buffer_[depth].front();
            std::queue<Message> empty;
            std::swap(message_buffer_[depth], empty);

            // If message size is not 0 perform mac check 
            if(message.data.size() != 0) {
                // Compute macs from shares 
                auto tp_macs = std::vector<Field>(mac_components.size());
                for(size_t i = 0; i < mac_components.size(); i++) {
                    tp_macs[i] =  preproc_.tp_key * message.data[i] - mac_components[i];
                }
                // Hash the macs 
                fieldDig tp_layer_dig = hashFields(tp_macs);

                fieldDig p_layer_dig(message.data.end() - 4, 
                              message.data.end()
                );

                if(tp_layer_dig != p_layer_dig) {
                    
                    std::cout << ": Received inconsistent shares from Party " << message.receiver_id << std::endl;
                }
            }

            for(size_t i = 0; i < mult_num; i++) {
                mult_all[i] = agg_values[i] + message.data[i];
            }
        }

        evaluateGatesAtDepthPartyRecv(depth, mult_all);
        
    }

    std::vector<Field> OnlineEvaluator::getOutputs() {
        std::vector<Field> outvals(circ_.outputs.size());
        if (circ_.outputs.empty())
            return outvals;

        // HP sends there out mask to all other parties 
        if (id_ == 0) {
            std::vector<Field> output_masks(circ_.outputs.size());
            for (size_t i = 0; i < circ_.outputs.size(); ++i) {
                Field outmask = output_masks[i];
                auto wout = circ_.outputs[i];
                outmask = preproc_.gates[wout]->mask.getValue();
                output_masks[i] = outmask;                
            }
            std::vector<std::future<void>> send_t;
            for (size_t i = 1; i <= nP_; ++i) {
                send_t.push_back(tpool_->enqueue([&,i]() {
                    network_->send(i, output_masks.data(), output_masks.size() * sizeof(Field));
                    network_->getSendChannel(i)->flush();
                }));
            }
            for(auto& t : send_t) {
                if (t.valid()) {
                    t.wait();
                }
            }
            return outvals;
        }
        else {
            std::vector<Field> output_masks(circ_.outputs.size());
            network_->recv(0, output_masks.data(), output_masks.size() * sizeof(Field));
            network_->getRecvChannel(0)->flush();
            for (size_t i = 0; i < circ_.outputs.size(); ++i) {
                Field outmask = output_masks[i];
                auto wout = circ_.outputs[i];
                outmask += preproc_.gates[wout]->mask.getValue();
                outvals[i] = wires_[wout] - outmask;
            }
            return outvals;
        }
    }

    std::vector<Field> OnlineEvaluator::evaluateCircuit(const std::unordered_map<wire_t, Field> &inputs) {
        setInputs(inputs);
        for (size_t i = 0; i < circ_.gates_by_level.size(); ++i)
            evaluateGatesAtDepth(i);

        return getOutputs();
    }
}; // namespace dmAsyncAsteriskGOD