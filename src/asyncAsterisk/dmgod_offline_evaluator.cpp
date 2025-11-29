#include "dmgod_offline_evaluator.h"

namespace dmAsyncAsteriskGOD {
    template <typename T>
  void print_vector(const std::vector<T>& v) {
      std::cout << "[";

      for (size_t i = 0; i < v.size(); ++i) {
          std::cout << v[i];
          if (i + 1 < v.size())
              std::cout << ", ";
      }

      std::cout << "]\n";
  }

  OfflineEvaluator::OfflineEvaluator(int nP, int id, int security_param, std::shared_ptr<NetIOMP> network1, 
    std::shared_ptr<NetIOMP> network2, LevelOrderedCircuit circ, int threads, uint64_t seed, bool run_async) 
    : nP_(nP), id_(id), security_param_(security_param), rgen_(id, nP, seed), network_(std::move(network1)), 
    network_ot_(std::move(network2)), circ_(std::move(circ)), preproc_(circ.num_gates), start_ot_(2, false), 
    chunk_size_(50000), inputToOPE(2), run_async_(true)
  {

    // Threadpool setup 
    if (run_async_) {
        tpool_ = std::make_shared<ThreadPool>(threads);
        tpool_minus_one_ = nullptr;
    } else {
        // 1 thread for sender party, remaining threads for digest distribution
        tpool_ = std::make_shared<ThreadPool>(1);
        tpool_minus_one_ = std::make_shared<ThreadPool>(nP_ - 1);
    }

    // Setup OT provider instances 
    if(run_async_) {
      if (id_ == 0) {
        for (size_t i = 1; i <= nP_; i++) {
          ot_.emplace_back(std::make_unique<OTProviderHA>(id_, i, network_ot_->getRecvChannel(i)));
          network_ot_->getRecvChannel(i)->flush();
        }
      }
      else {
        ot_.emplace_back(std::make_unique<OTProviderHA>(id_, 0, network_ot_->getSendChannel(0)));
        network_ot_->getSendChannel(0)->flush();
      }
    } else {
      if (id_ == 0) {
        for (size_t pid = 1; pid <= nP_; pid++) {
          if (pid == SYNC_SENDER_PID_) {
            ot_.emplace_back(std::make_unique<OTProviderHA>(id_, SYNC_SENDER_PID_, network_ot_->getRecvChannel(SYNC_SENDER_PID_)));
            network_ot_->getRecvChannel(SYNC_SENDER_PID_)->flush();
          } else {
            std::make_unique<OTProviderHA>(id_, pid, network_ot_->getRecvChannel(pid));
            network_ot_->getSendChannel(pid)->flush();
          }
        }
      }
      else {
        // Create OT instance to HP 
        ot_.emplace_back(std::make_unique<OTProviderHA>(id_, 0, network_ot_->getSendChannel(0)));
        network_ot_->getSendChannel(0)->flush();
      }
    }

    // Worker thread setup 
    if (id_ != 0) 
        return;     

    static ZZ_pContext ZZ_p_ctx;
    ZZ_p_ctx.save();

    auto spawn_worker = [&](size_t pid, ThreadPool& pool) {
        pool.enqueue([&, pid]() {
            ZZ_p_ctx.restore();
            for (size_t count = 0; count < 2; count++) {
                std::vector<Field> sharesVec;
                std::vector<fieldDig> chunk_digs;

                // Wait for the signal to start this OT batch
                {
                    std::unique_lock<std::mutex> lock(mtx_);
                    cv_start_ot_[count].wait(lock, [&]() { return start_ot_[count]; });
                }

                sharesVec.resize(inputToOPE[count].size());
                fieldDig current_dig;

                // Chunked OPE
                for (size_t start = 0; start < inputToOPE[count].size(); start += chunk_size_) {
                    size_t end = std::min(start + chunk_size_, inputToOPE[count].size());
                    std::vector<Field> chunk(inputToOPE[count].begin() + start,
                                             inputToOPE[count].begin() + end);

                    auto out = ot_[pid - 1]->multiplyRecv(chunk, current_dig);
                    std::copy(out.begin(), out.end(), sharesVec.begin() + start);
                    chunk_digs.push_back(current_dig);
                }

                {
                    std::lock_guard<std::mutex> lock(mtx_);
                    offline_message_buffer_[count].push({pid, sharesVec, chunk_digs});
                }
                cv_.notify_one();
            }
        });
    };


    // Async mode: spawn worker for every party 
    if (run_async_) {
        for (size_t pid = 1; pid <= nP_; pid++)
            spawn_worker(pid, *tpool_);
        return;
    }

    // Sync mode: only spawn worker for SYNC_SENDER_PID_ 
    spawn_worker(SYNC_SENDER_PID_, *tpool_);
  }

  OfflineEvaluator::~OfflineEvaluator() {
    tpool_.reset();
    if(!run_async_) {
      tpool_minus_one_.reset();
    }
  }

  void OfflineEvaluator::keyGen()  {
    if(id_ == 0) {
      randomizeZZp(rgen_.p0(), key_, sizeof(Field));
      preproc_.setTPKey(key_);
    }
  }

  void OfflineEvaluator::randSS(int pid, RandGenPool& rgen, TwoShare<Field>& share, Field& mask_share_zero, bool isOutputWire) {
    // TP 
    if(pid == 0) {      
      Field valSh;
      randomizeZZp(rgen.self(), valSh, sizeof(Field));
      share.setValue(valSh);
    }
    // Parties 
    else {
      Field valSh;
      randomizeZZp(rgen.all_minus_0(), valSh, sizeof(Field));
      share.setValue(valSh);
    }
  }

  void OfflineEvaluator::randomShareSecret(int pid, RandGenPool& rgen, const TwoShare<Field>& share1, const TwoShare<Field>& share2, 
    TwoShare<Field>& prodShare, std::vector<Field>& inputToOPE) {
    auto share1_val = share1.getValue();
    auto share2_val = share2.getValue();
    // Swap shares to compute cross terms in OPE 
    if(pid == 0) {
      std::swap(share1_val, share2_val);
    }
    inputToOPE.push_back(share1_val);
    inputToOPE.push_back(share2_val);
    prodShare.setValue(Field(0));
  }  

  void OfflineEvaluator::randSSWithParty(int pid, int dealer, RandGenPool& rgen, TwoShare<Field>& share, Field& secret) {        
    secret = Field(0);
    Field valSh;
    // TP 
    if(pid == 0) {  
      // TP is not dealer 
      if (pid != dealer) {
          // Sample shareTP together with dealer
          randomizeZZp(rgen.pij(dealer, dealer), valSh, sizeof(Field));
          share.setValue(valSh);
      }
      // TP is dealer 
      else {
        // Sample shareTP alone  
        Field valSh;
        randomizeZZp(rgen.p0(), valSh, sizeof(Field));
        share.setValue(valSh);
        secret += valSh;
        // Sample shareP with all Parties 
        Field val;
        randomizeZZp(rgen.all(), val, sizeof(Field));
        secret += val;
      }
    }
    // Parties 
    else {
      // TP is dealer 
      if (dealer == 0) {
        // Sample shareP with TP
        randomizeZZp(rgen.all(), valSh, sizeof(Field));
        share.setValue(valSh);
      }
      // Party pid is not the dealer 
      else if (pid != dealer) {
        Field valSh;
        randomizeZZp(rgen.all_minus_0(), valSh, sizeof(Field));
        share.setValue(valSh);
      }
      // Party pid is the dealer 
      else {
        Field valSh;
        randomizeZZp(rgen.all_minus_0(), valSh, sizeof(Field));
        share.setValue(valSh);
        secret += valSh;
        // Sample shareTP together with TP
        Field val;
        randomizeZZp(rgen.pij(dealer, dealer), val, sizeof(Field));
        secret += val;
      }
    }
  }

  bool OfflineEvaluator::verifyOPEMsgs(std::vector<fieldDig> chunk_digs, Field sender_id) {
    constexpr size_t digest_size = 4;   
    const size_t num_chunks = chunk_digs.size();
    const size_t total_comm = 1 + digest_size * num_chunks;

    /**
     * recv_buf/send_buf layout:
     *  digest chunk 1 | digest chunk 2 | digest chunk 3 | ... | digest chunk n | sender_id 
     */

    if(id_ == 0) {
      std::vector<Field> send_buf;
      send_buf.reserve(total_comm);
      
      for (const auto& current_dig : chunk_digs) {
        send_buf.insert(send_buf.end(), current_dig.begin(), current_dig.end());
      }
      send_buf.push_back(sender_id);

      // Send send_buf to each party 
      std::vector<std::future<void>> send_t;    
      auto& pool = run_async_ ? *tpool_ : *tpool_minus_one_;

      for (size_t pid = 1; pid <= nP_; ++pid) {

          // In sync mode, skip sending to SYNC_SENDER_PID_
          if (!run_async_ && pid == SYNC_SENDER_PID_) {
              continue;
          }

          send_t.push_back(pool.enqueue([&, pid]() {
              network_->send(pid, send_buf.data(), sizeof(Field) * total_comm);
              network_->getSendChannel(pid)->flush();
          }));
      }
      
      for(auto& t : send_t) {
          if (t.valid()) {
              t.wait();
          }
      }
    }
    else {
      if(!run_async_ && Field(id_) == SYNC_SENDER_PID_) {
        return true;
      }
      // Receive buffer
      std::vector<Field> recv_buf(total_comm);
      network_->recv(0, recv_buf.data(), recv_buf.size() * sizeof(Field));
      network_->getRecvChannel(0)->flush();

      // Extract sender from buffer 
      sender_id = recv_buf.at(recv_buf.size() - 1);

      for (size_t start = 0; start < total_comm - 1; start += digest_size) {
        // Extract received digest from buffer
        const fieldDig recv_digest(
            recv_buf.begin() + start,
            recv_buf.begin() + start + digest_size
        );

        // My digest for chunk i
        const fieldDig my_digest = chunk_digs[start / digest_size];

        // Check mismatch
        if (Field(id_) != sender_id && my_digest != recv_digest) {
            std::cout << "Party " << id_
                      << " has identified Party " << sender_id
                      << " as a cheater during OPE!" << std::endl;
            return false;
        }
      }
    }
    return true;
  }

  void OfflineEvaluator::runOPE(std::vector<Field>& inputToOPE, std::vector<Field>& outputOfOPE, size_t count) {
    std::vector<fieldDig> chunk_digs;
    Field sender_id = Field(-1);

    if(id_ != 0) {
      for (size_t start = 0; start < inputToOPE.size(); start += chunk_size_) {
          size_t end = std::min(start + chunk_size_, inputToOPE.size());
          std::vector<Field> chunk(inputToOPE.begin() + start, inputToOPE.begin() + end);
          fieldDig current_dig;
          auto chunk_output = ot_[0]->multiplySend(chunk, rgen_.all_minus_0(), current_dig, run_async_, id_);
          outputOfOPE.insert(outputOfOPE.end(), chunk_output.begin(), chunk_output.end());
          chunk_digs.push_back(current_dig);
      }
    } else {
      {
        std::lock_guard<std::mutex> lock(mtx_);
        start_ot_[count] = true;
      }
      cv_start_ot_[count].notify_all();
      {
        std::unique_lock<std::mutex> lock(mtx_);
        cv_.wait(lock, [&]() { return offline_message_buffer_[count].size() >= 1; });
      }
      Offline_Message OPE_res = offline_message_buffer_[count].front();
      auto receiver_pid = OPE_res.receiver_id;
      std::queue<Offline_Message> empty;
      std::swap(offline_message_buffer_[count], empty);
      outputOfOPE = OPE_res.data;
      chunk_digs = OPE_res.chunk_digs;
      sender_id = Field(receiver_pid);
    }
    verifyOPEMsgs(chunk_digs, sender_id);
  }

  void OfflineEvaluator::multSS(const Field& share1_val, const Field& share2_val, Field& output_val, 
    const std::vector<Field>& outputOfOPE, size_t& idx_outputOfOPE) {
      output_val = share1_val * share2_val + outputOfOPE[idx_outputOfOPE] + outputOfOPE[idx_outputOfOPE+1];
      idx_outputOfOPE += 2;
  }

  void OfflineEvaluator::prepareMaskValues(const std::unordered_map<wire_t,int>& input_pid_map) {
    std::vector<Field> buffer; 
    size_t idx_buffer=0; 
    size_t buffer_num=0;

    for (const auto& level : circ_.gates_by_level) {
      for (const auto& gate : level) {
        switch (gate->type) {
          case GateType::kInp: {
            auto pregate = std::make_unique<PreprocInput<Field>>();      
            auto pid = input_pid_map.at(gate->out);
            pregate->pid = pid;
            // Give mask of input gate in the clear to dealer 
            randSSWithParty(id_, pid, rgen_, pregate->mask, pregate->mask_value);      
            preproc_.gates[gate->out] = std::move(pregate);                
            break;
          }

          case GateType::kAdd: {
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            const auto& mask_in1 = preproc_.gates[g->in1]->mask;
            const auto& mask_in2 = preproc_.gates[g->in2]->mask;
            preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask_in1 + mask_in2));    
            break;
          }

          case GateType::kConstAdd: {
            const auto* g = static_cast<ConstOpGate<Field>*>(gate.get());
            const auto& mask = preproc_.gates[g->in]->mask;
            preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask));
            break;
          }

          case GateType::kConstMul: {
            const auto* g = static_cast<ConstOpGate<Field>*>(gate.get());
            const auto& mask = preproc_.gates[g->in]->mask * g->cval;
            preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask));
            break;
          }

          case GateType::kSub: {
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            const auto& mask_in1 = preproc_.gates[g->in1]->mask;
            const auto& mask_in2 = preproc_.gates[g->in2]->mask;
            preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask_in1 - mask_in2));
            break;
          }

          case GateType::kMul: {
            // Create the output wire mask share and initialize it later 
            preproc_.gates[gate->out] = std::make_unique<PreprocMultGate<Field>>();
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            const auto& mask_in1 = preproc_.gates[g->in1]->mask;
            const auto& mask_in2 = preproc_.gates[g->in2]->mask;
            TwoShare<Field> mask_out; 
            Field mask_share_zero = Field(0);
            bool isOutputWire = false;
            if (std::find(circ_.outputs.begin(), circ_.outputs.end(),g->out)!=circ_.outputs.end())
              isOutputWire = true;
            // Generate a random mask for the output wire 
            randSS(id_, rgen_, mask_out, mask_share_zero, isOutputWire);
            TwoShare<Field> mask_product;
            // Compute cross terms with HP for mask_product = mask_in1 * mask_in2  
            randomShareSecret(id_, rgen_, mask_in1, mask_in2, mask_product, inputToOPE[0]);
            preproc_.gates[gate->out] = std::move(std::make_unique<PreprocMultGate<Field>> (mask_out, mask_product));
            break;
          }

          default: {
            break;
          }
        }
      }
    }
  
    // Run MSSR OLE for multiplication triples 
    std::vector<Field> outputOfOPE;
    size_t idx_outputOfOPE = 0;
    
    runOPE(inputToOPE[0], outputOfOPE, 0);

    // After OLEs compute output mask on multiplication gates 
    for (const auto& level : circ_.gates_by_level) {
      for (const auto& gate : level) {
        switch (gate->type) {
          case GateType::kMul: {
            auto* g = static_cast<FIn2Gate*>(gate.get()); 
            auto *pre_mul = static_cast<PreprocMultGate<Field> *>(preproc_.gates[gate->out].get());
            // Compute the share of the masked output 
            auto mask_in1_in2_product_val = pre_mul->mask_prod.getValue();
            auto mask_in1_val = preproc_.gates[g->in1]->mask.getValue();
            auto mask_in2_val = preproc_.gates[g->in2]->mask.getValue();
            // Compute product share 
            multSS(mask_in1_val, mask_in2_val, mask_in1_in2_product_val, outputOfOPE, idx_outputOfOPE);
            pre_mul->mask_prod.setValue(mask_in1_in2_product_val);
            break;
          }
  
          default: {
            break;
          }
        }
      }
    }
    outputOfOPE.clear();
    outputOfOPE.shrink_to_fit();
  }

  void OfflineEvaluator::prepareMaskMACs() {
    for (const auto& level : circ_.gates_by_level) {
      for (const auto& gate : level) {
        switch (gate->type) {
          case GateType::kInp: {
            if (id_!=0) {      
              auto *pre_input = static_cast<PreprocInput<Field> *>(preproc_.gates[gate->out].get());
              auto mask = pre_input->mask.getValue();
              inputToOPE[1].push_back(mask);
            }
            else {
              inputToOPE[1].push_back(key_);
            }
            break;
          }

          case GateType::kMul: {
            if (id_!=0) {     
              auto *pre_gate = preproc_.gates[gate->out].get();
              auto *pre_mul = static_cast<PreprocMultGate<Field> *>(preproc_.gates[gate->out].get());
              auto mask = pre_gate->mask.getValue();
              auto mask_prod = pre_mul->mask_prod.getValue(); 
              inputToOPE[1].push_back(mask);
              inputToOPE[1].push_back(mask_prod);
            }
            else {
              inputToOPE[1].push_back(key_);
              inputToOPE[1].push_back(key_);
            }
            break;
          }

          default: {
            break;
          }
        }
      }
    }

    std::vector<Field> outputOfOPE;
    size_t idx_outputOfOPE = 0;

    runOPE(inputToOPE[1], outputOfOPE, 1);

    for (const auto& level : circ_.gates_by_level) {
      for (const auto& gate : level) {
        switch (gate->type) {
          case GateType::kInp: {
            auto *pre_input = static_cast<PreprocInput<Field> *>(preproc_.gates[gate->out].get());
            pre_input->mask.setMACComponent(outputOfOPE[idx_outputOfOPE++]);                
            break;
          }

          case GateType::kAdd: {
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            const auto& mask_in1 = preproc_.gates[g->in1]->mask;
            const auto& mask_in2 = preproc_.gates[g->in2]->mask;
            preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask_in1 + mask_in2));    
            break;
          }

          case GateType::kConstAdd: {
            const auto* g = static_cast<ConstOpGate<Field>*>(gate.get());
            const auto& mask = preproc_.gates[g->in]->mask + g->cval;
            preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask));
            break;
          }

          case GateType::kConstMul: {
            const auto* g = static_cast<ConstOpGate<Field>*>(gate.get());
            const auto& mask = preproc_.gates[g->in]->mask * g->cval;
            preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask));
            break;
          }

          case GateType::kSub: {
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            const auto& mask_in1 = preproc_.gates[g->in1]->mask;
            const auto& mask_in2 = preproc_.gates[g->in2]->mask;
            preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask_in1 - mask_in2));
            break;
          }

          case GateType::kMul: {
            auto *pre_gate = preproc_.gates[gate->out].get();
            auto *pre_mul = static_cast<PreprocMultGate<Field> *>(preproc_.gates[gate->out].get());
            pre_gate->mask.setMACComponent(outputOfOPE[idx_outputOfOPE++]);
            pre_mul->mask_prod.setMACComponent(outputOfOPE[idx_outputOfOPE++]);
            break;
          }
  
          default: {
            break;
          }
        }
      }
    }
    outputOfOPE.clear();
    outputOfOPE.shrink_to_fit();      
  }
    
  void OfflineEvaluator::setWireMasks(const std::unordered_map<wire_t,int>& input_pid_map) {      
    keyGen();
    prepareMaskValues(input_pid_map);
    prepareMaskMACs(); 
  }

  PreprocCircuit<Field> OfflineEvaluator::run(const std::unordered_map<wire_t, int>& input_pid_map) {
      setWireMasks(input_pid_map);
      return std::move(preproc_);    
  }

}; // namespace dmAsyncAsteriskGOD