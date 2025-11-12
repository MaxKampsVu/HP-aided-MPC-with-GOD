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
            // Create threads for parties 
            for(size_t pid = 1; pid <= nP_; pid++) {
                tpool_->enqueue([&,pid]() {                    
                    ZZ_p_ctx.restore();
                    {
                        std::unique_lock<std::mutex> lock(mtx_);
                        cv_start_recv_.wait(lock, [&]() { return start_recv_; });
                    }
                    // TODO: remove the -1 to consider output gates!!!!!!!!!!
                    for (size_t depth = 0; depth <= circ_.gates_by_level.size() - 1; depth++) {
                        size_t total_comm;
                        network_->recv(pid, &total_comm, sizeof(size_t));
                        std::vector<Field> online_comm_to_HP(total_comm);
                        network_->recv(pid, online_comm_to_HP.data(), sizeof(Field) * total_comm);
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
                    // TODO: what is this? 
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
                    }
                }
                
                wires_[g->out] = q_val_[g->out];
            }
        }
    }

    void OnlineEvaluator::setRandomInputs() {
        for (auto &g : circ_.gates_by_level[0]) {
            if (g->type == GateType::kInp) {
                randomizeZZp(rgen_.all(), wires_[g->out], sizeof(Field));
            }
        }
    }

    void OnlineEvaluator::evaluateGatesAtDepthPartySend(size_t depth, std::vector<Field> &mult_nonTP) {
        for (auto &gate : circ_.gates_by_level[depth]) {
            switch (gate->type) {
                case GateType::kMul: {
                    auto *g = static_cast<FIn2Gate *>(gate.get());
                    auto *pre_out = static_cast<PreprocMultGate<Field> *>(preproc_.gates[g->out].get());
                    q_val_[g->out] = 0;

                    auto m_in1_val = preproc_.gates[g->in1]->mask.getValue();
                    auto m_in2_val = preproc_.gates[g->in2]->mask.getValue();
                    auto m_prod_val = pre_out->mask_prod.getValue();
                    auto m_out_val = pre_out->mask.getValue();


                    auto q_share = Field(0);

                    if (id_==0) {
                        auto q_share = wires_[g->in1] * wires_[g->in2] + m_prod_val + m_out_val - m_in1_val * wires_[g->in2] - m_in2_val * wires_[g->in1];    
                    } else {
                        auto q_share = m_prod_val + m_out_val - m_in1_val * wires_[g->in2] - m_in2_val * wires_[g->in1];    
                    }                        
                    
                    //q_sh_[g->out] = q_share;
                    mult_nonTP.push_back(q_share);
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

                // Set the reconstructed output of the mul gates to what HP has send 
                case GateType::kMul: {
                    auto *g = static_cast<FIn2Gate *>(gate.get());
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

                default:
                    break;
            }
        }

        // Number of multiplication gates at the layer 
        size_t total_comm = mult_num;
        std::vector<Field> mult_all(mult_num);        
        std::vector<Field> agg_values(total_comm);

        // Parties and TP prepare their output share for each multiplication gate at the current layer  
        evaluateGatesAtDepthPartySend(depth, mult_nonTP);
        if (id_ != 0) {            
            std::vector<Field> online_comm_to_HP(total_comm);

            // Write output share on multiplication gate for the layer to buffer 
            for (size_t i = 0; i < mult_num; i++) {
                online_comm_to_HP[i] = mult_nonTP[i];
            }
            // Send to HP size of buffer 
            network_->send(0, &total_comm, sizeof(size_t));
            network_->getSendChannel(0)->flush();
            // Send to HP the shares 
            network_->send(0, online_comm_to_HP.data(), sizeof(Field) * total_comm, true);
            network_->getSendChannel(0)->flush();
        }
        else {
            // HP prepares buffer for reconstructed output masks 
            for (size_t i=0; i<mult_num; i++) {
                agg_values[i] = mult_nonTP[i];
            }

            // Sends reconstructed output masks to parties the moment he has them 
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

        // Parties receive recombined shares from HP
        if (id_ != 0) {
            network_->recv(0, agg_values.data(), sizeof(Field) * total_comm);
            for(size_t i = 0; i < mult_num; i++) {
                mult_all[i] = agg_values[i] + mult_nonTP[i];
            }
        }
        // HP
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

            // HP reconstructs values 
            for(size_t i = 0; i < mult_num; i++) {
                mult_all[i] = agg_values[i] + message.data[i];
            }
        }

        evaluateGatesAtDepthPartyRecv(depth, mult_all);
    }

    std::vector<Field> OnlineEvaluator::getOutputs() {
        std::vector<Field> outvals(circ_.outputs.size());
        //if (circ_.outputs.empty())
            return outvals;

        // if (id_ == 0) {
        //     std::vector<Field> output_masks(circ_.outputs.size());
        //     for (size_t i = 0; i < circ_.outputs.size(); ++i) {
        //         auto wout = circ_.outputs[i];
        //         Field outmask = Field(0);
        //         std::vector<Field> valVec= preproc_.gates[wout]->mask.getValue();
        //         for (size_t j=0; j<nP_; j++) 
        //             outmask += valVec[j];
        //         output_masks[i] = outmask;
        //         outmask += preproc_.gates[wout]->mask_share_zero;                
        //         outvals[i] = wires_[wout] - outmask;
        //     }
        //     std::vector<std::future<void>> send_t;
        //     for (size_t i = 1; i <= nP_; ++i) {
        //         send_t.push_back(tpool_->enqueue([&,i]() {
        //             network_->send(i, output_masks.data(), output_masks.size() * sizeof(Field));
        //             network_->getSendChannel(i)->flush();
        //         }));
        //     }
        //     for(auto& t : send_t) {
        //         if (t.valid()) {
        //             t.wait();
        //         }
        //     }
        //     return outvals;
        // }
        // else {
        //     std::vector<Field> output_masks(circ_.outputs.size());
        //     network_->recv(0, output_masks.data(), output_masks.size() * sizeof(Field));
        //     for (size_t i = 0; i < circ_.outputs.size(); ++i) {
        //         Field outmask = output_masks[i];
        //         auto wout = circ_.outputs[i];
        //         outmask += preproc_.gates[wout]->mask.getValue()[1];
        //         outvals[i] = wires_[wout] - outmask;
        //     }
        //     return outvals;
        // }
    }

    std::vector<Field> OnlineEvaluator::evaluateCircuit(const std::unordered_map<wire_t, Field> &inputs) {
        setInputs(inputs);
        for (size_t i = 0; i < circ_.gates_by_level.size(); ++i)
            evaluateGatesAtDepth(i);

        return getOutputs();
    }
}; // namespace dmAsyncAsteriskGOD
