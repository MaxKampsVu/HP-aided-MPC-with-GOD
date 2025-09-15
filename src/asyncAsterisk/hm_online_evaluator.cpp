#include "hm_online_evaluator.h"

namespace hmAsyncAsterisk {
    OnlineEvaluator::OnlineEvaluator(int nP, int id, int security_param, std::shared_ptr<NetIOMP> network, PreprocCircuit<Field> preproc, 
        LevelOrderedCircuit circ, int threads, uint64_t seed) 
        : nP_(nP), th_((nP_-1)/2), id_(id), security_param_(security_param), rgen_(id, nP, seed), network_(std::move(network)), 
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
                    for (size_t depth = 0; depth <= circ_.gates_by_level.size(); depth++) {
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
                auto *pre_input = static_cast<PreprocInput<Field> *>(preproc_.gates[g->out].get());
                auto pid = pre_input->pid;

                if (id_ == 0) {
                    Field padded_q_val = Field(0);
                    network_->recv(pid, &padded_q_val, sizeof(Field));
                    std::vector<std::future<void>> send_t;
                    for (size_t i = 1; i <= nP_; i++) {
                        if (i == pid)
                            continue;
                        send_t.push_back(tpool_->enqueue([&,i]() {
                            network_->send(i, &padded_q_val, sizeof(Field));
                            network_->getSendChannel(i)->flush();
                        }));
                    }
                    for(auto& t : send_t) {
                        if (t.valid()) {
                            t.wait();
                        }
                    }
                }
                else {
                    Field r = Field(0);
                    randomizeZZp(rgen_.all_minus_0(), r, sizeof(Field));
                    if (pid == id_) {
                        q_val_[g->out] = pre_input->mask_value + inputs.at(g->out);
                        Field padded_q_val = q_val_[g->out] + r;
                        network_->send(0, &padded_q_val, sizeof(Field));
                        network_->getSendChannel(0)->flush();
                        wires_[g->out] = q_val_[g->out];
                    }
                    else {
                        Field padded_q_val = Field(0);
                        network_->recv(0, &padded_q_val, sizeof(Field));
                        q_val_[g->out] = padded_q_val - r;
                        wires_[g->out] = q_val_[g->out];
                    }
                }
            }
        }
    }

    void OnlineEvaluator::setRandomInputs() {
        for (auto &g : circ_.gates_by_level[0]) {
            if (g->type == GateType::kInp) {
                randomizeZZp(rgen_.all_minus_0(), wires_[g->out], sizeof(Field));
            }
        }
    }

    void OnlineEvaluator::evaluateGatesAtDepthPartySend(size_t depth, std::vector<Field> &mult_nonTP, std::vector<Field> &r_mult_pad) {
        for (auto &gate : circ_.gates_by_level[depth]) {
            switch (gate->type) {
                case GateType::kMul: {
                    auto *g = static_cast<FIn2Gate *>(gate.get());
                    q_val_[g->out] = 0;

                    if (id_ != 0) {
                        Field r_mul = Field(0);
                        if (std::find(circ_.outputs.begin(), circ_.outputs.end(),g->out)==circ_.outputs.end())
                            randomizeZZp(rgen_.all_minus_0(), r_mul, sizeof(Field));
                        r_mult_pad.push_back(r_mul);
                        auto &m_in1 = preproc_.gates[g->in1]->mask;
                        auto &m_in2 = preproc_.gates[g->in2]->mask;
                        auto *pre_out = static_cast<PreprocMultGate<Field> *>(preproc_.gates[g->out].get());
                        auto q_share = pre_out->mask + pre_out->mask_prod - m_in1 * wires_[g->in2] - m_in2 * wires_[g->in1];
                        q_share.add(wires_[g->in1] * wires_[g->in2] + r_mul);
                        mult_nonTP.push_back(q_share.valueAt());
                        q_sh_[g->out] = q_share;
                    }

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

    void OnlineEvaluator::evaluateGatesAtDepthPartyRecv(size_t depth, std::vector<Field> mult_all, std::vector<Field> r_mult_pad) {
        size_t idx_mult = 0;

        for (auto &gate : circ_.gates_by_level[depth]) {
            switch (gate->type) {
                case GateType::kAdd: {
                    auto *g = static_cast<FIn2Gate *>(gate.get());
                    if (id_ != 0)
                        wires_[g->out] = wires_[g->in1] + wires_[g->in2];
                    q_val_[g->out] = 0;
                    break;
                }

                case GateType::kSub: {
                    auto *g = static_cast<FIn2Gate *>(gate.get());
                    if (id_ != 0)
                        wires_[g->out] = wires_[g->in1] - wires_[g->in2];
                    q_val_[g->out] = 0;
                    break;
                }

                case GateType::kConstAdd: {
                    auto *g = static_cast<ConstOpGate<Field> *>(gate.get());
                    if (id_ != 0)
                        wires_[g->out] = wires_[g->in] + g->cval;
                    break;
                }

                case GateType::kConstMul: {
                    auto *g = static_cast<ConstOpGate<Field> *>(gate.get());
                    if (id_ != 0)
                        wires_[g->out] = wires_[g->in] * g->cval;
                    break;
                }

                case GateType::kMul: {
                    auto *g = static_cast<FIn2Gate *>(gate.get());
                    if (id_ != 0) {
                        q_val_[g->out] = mult_all[idx_mult];
                        wires_[g->out] = q_val_[g->out] - r_mult_pad[idx_mult];
                    }
                    else {
                        wires_[g->out] = mult_all[idx_mult];
                    }
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
        std::vector<Field> r_mult_pad;

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

        size_t total_comm = mult_num;

        if (id_ != 0) {
            evaluateGatesAtDepthPartySend(depth, mult_nonTP, r_mult_pad);
            std::vector<Field> online_comm_to_HP(total_comm);

            for (size_t i = 0; i < mult_num; i++) {
                online_comm_to_HP[i] = mult_nonTP[i];
            }
            network_->send(0, &total_comm, sizeof(size_t));
            network_->getSendChannel(0)->flush();
            network_->send(0, online_comm_to_HP.data(), sizeof(Field) * total_comm, true);
            network_->getSendChannel(0)->flush();
        }
        else {
            std::vector<Field> agg_values(total_comm, Field(0));
            {
                std::lock_guard<std::mutex> lock(mtx_);
                start_recv_ = true;
            }
            cv_start_recv_.notify_all();
            {
                std::unique_lock<std::mutex> lock(mtx_);
                cv_.wait(lock, [&]() { return message_buffer_[depth].size() >= th_ + 1; });
            }

            std::vector<Message> messages;
            for (size_t i = 0; i < th_ + 1; ++i) {
                messages.push_back(message_buffer_[depth].front());
                message_buffer_[depth].pop();
            }
            std::queue<Message> empty;
            std::swap(message_buffer_[depth], empty);

            std::vector<Field> evalPoints(th_ + 1);
            for (size_t i = 0; i < th_ + 1; ++i) {
                evalPoints[i] = Field(messages[i].receiver_id);
            }
            for (size_t i=0; i<total_comm; i++) {
                std::vector<Field> q_shares(th_ + 1);
                for (size_t j=0; j<=th_; j++) {
                    q_shares[j] = messages[j].data[i];
                }
                fieldPoly poly_masked_val = reconstructPolynomial(evalPoints, q_shares);
                agg_values[i] = coeff(poly_masked_val, 0);
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
            std::vector<Field> mult_all(mult_num);
            for(size_t i = 0; i < mult_num; i++) {
                mult_all[i] = agg_values[i];
            }
            evaluateGatesAtDepthPartyRecv(depth, mult_all, r_mult_pad);
        }

        if (id_ != 0) {
            std::vector<size_t> completed_pids(th_ + 1);
            std::vector<Field> agg_values(total_comm);
            network_->recv(0, agg_values.data(), sizeof(Field) * total_comm);
            std::vector<Field> mult_all(mult_num);
            for(size_t i = 0; i < mult_num; i++) {
                mult_all[i] = agg_values[i];
            }
            evaluateGatesAtDepthPartyRecv(depth, mult_all, r_mult_pad);
        }
    }

    bool OnlineEvaluator::MACVerification() {
        block cc_key[2];
        if (id_ == 0) {
            rgen_.self().random_block(cc_key, 2);
            std::vector<std::future<void>> send_t;
            for (size_t i = 1; i <= nP_; i++) {
                send_t.push_back(tpool_->enqueue([&,i]() {
                    network_->send(i, cc_key, 2 * sizeof(block));
                    network_->getSendChannel(i)->flush();
                }));
            }
            for(auto& t : send_t) {
                if (t.valid()) {
                    t.wait();
                }
            }
        }
        else {
            network_->recv(0, cc_key, 2 * sizeof(block));
        }
        PRG prg;
        prg.reseed(cc_key);
        Field res = Field(0);
        if (id_ != 0) {
            Field key = preproc_.gates[0]->mask.keySh();
            Field omega = Field(0);
            std::unordered_map<wire_t, Field> rho; 
            for (size_t i = 0; i < circ_.gates_by_level.size(); ++i) {
                for (auto &gate : circ_.gates_by_level[i]) {
                    switch (gate->type) {
                        case GateType::kMul: {
                            auto *g = static_cast<FIn2Gate *>(gate.get());
                            randomizeZZp(prg, rho[g->out], sizeof(Field));
                            omega += rho[g->out] * (q_val_[g->out] * key - q_sh_[g->out].tagAt());
                        }

                        case GateType::kConstAdd:
                        case GateType::kConstMul:
                        case GateType::kAdd:
                        case GateType::kSub:
                            break;
                        
                        default:
                            break;
                    }
                }
            }
            Field rho_zero;
            randomizeZZp(prg, rho_zero, sizeof(Field));
            fieldPoly random_poly = randomPolynomial(rgen_.all_minus_0(), th_, Field(0));
            omega += rho_zero*eval(random_poly, Field(id_));
            size_t total_comm = 1;
            network_->send(0, &total_comm, sizeof(size_t));
            network_->getSendChannel(0)->flush();
            network_->send(0, &omega, sizeof(Field), true);
            network_->getSendChannel(0)->flush();
        }
        else {
            size_t depth = circ_.gates_by_level.size();
            {
                std::unique_lock<std::mutex> lock(mtx_);
                cv_.wait(lock, [&]() { return message_buffer_[depth].size() >= th_ + 1; });
            }

            std::vector<Message> messages;
            for (size_t i = 0; i < th_ + 1; ++i) {
                messages.push_back(message_buffer_[depth].front());
                message_buffer_[depth].pop();
            }
            std::queue<Message> empty;
            std::swap(message_buffer_[depth], empty);

            std::vector<Field> evalPoints(th_ + 1);
            std::vector<Field> omega_shares(th_ + 1);
            for (size_t i = 0; i < th_ + 1; ++i) {
                evalPoints[i] = Field(messages[i].receiver_id);
                omega_shares[i] = messages[i].data[0];
            }
            fieldPoly poly_omega = reconstructPolynomial(evalPoints, omega_shares);
            res = coeff(poly_omega, 0);

            std::vector<std::future<void>> send_t;
            for (size_t i = 1; i <= nP_; i++) {
                send_t.push_back(tpool_->enqueue([&,i]() {
                    network_->send(i, &res, sizeof(Field));
                    network_->getSendChannel(i)->flush();
                }));
            }
            for(auto& t : send_t) {
                if (t.valid()) {
                    t.wait();
                }
            }
        }
        if (id_ != 0) {
            network_->recv(0, &res, sizeof(Field));
        }

        if (res == 0)
            return true;
        else
            return false;
    }

    std::vector<Field> OnlineEvaluator::getOutputs() {
        std::vector<Field> outvals(circ_.outputs.size());
        if (circ_.outputs.empty())
            return outvals;

        if (id_ == 0) {
            std::vector<Field> output_masks(circ_.outputs.size());
            for (size_t i = 0; i < circ_.outputs.size(); ++i) {
                auto wout = circ_.outputs[i];
                Field outmask = preproc_.gates[wout]->tpmask.secret();
                output_masks[i] = outmask;
                outvals[i] = wires_[wout] - outmask;
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
            for (size_t i = 0; i < circ_.outputs.size(); ++i) {
                Field outmask = output_masks[i];
                auto wout = circ_.outputs[i];
                outvals[i] = wires_[wout] - outmask;
            }
            return outvals;
        }
    }

    std::vector<Field> OnlineEvaluator::evaluateCircuit(const std::unordered_map<wire_t, Field> &inputs) {
        setInputs(inputs);
        for (size_t i = 0; i < circ_.gates_by_level.size(); ++i)
            evaluateGatesAtDepth(i);

        if (MACVerification())
            return getOutputs();
        else {
            std::cout << "Malicious Activity Detected!!! MAC verification failed!!!" << std::endl;
            std::vector<Field> abort(circ_.outputs.size(), Field(0));
            return abort;
        }
    }
}; // namespace hmAsyncAsterisk
