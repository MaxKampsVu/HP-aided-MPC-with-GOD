#include "dmgod_offline_evaluator.h"

namespace dmAsyncAsteriskGOD {
  OfflineEvaluator::OfflineEvaluator(int nP, int id, int security_param, std::shared_ptr<NetIOMP> network1, 
    std::shared_ptr<NetIOMP> network2, LevelOrderedCircuit circ, int threads, uint64_t seed) 
    : nP_(nP), id_(id), security_param_(security_param), rgen_(id, nP, seed), network_(std::move(network1)), 
    network_ot_(std::move(network2)), circ_(std::move(circ)), preproc_(circ.num_gates), start_ot_(2, false), 
    chunk_size_(50000), inputToOPE(2), run_async_(true)
  {
    run_async_ ? setupASync(threads) : setupSync();
  }

  void OfflineEvaluator::setupASync(int threads) {
    tpool_ = std::make_shared<ThreadPool>(threads);
    if (id_ == 0) {
      // Create OT instance for each Party
      for (size_t i = 1; i <= nP_; i++) {
        ot_.emplace_back(std::make_unique<OTProviderHA>(id_, i, network_ot_->getRecvChannel(i)));
        network_ot_->getRecvChannel(i)->flush();
      }
    }
    else {
      // Create OT instance to HP 
      ot_.emplace_back(std::make_unique<OTProviderHA>(id_, 0, network_ot_->getSendChannel(0)));
      network_ot_->getSendChannel(0)->flush();
    }

    if (id_ == 0) {
      static ZZ_pContext ZZ_p_ctx;
      ZZ_p_ctx.save();
      // Start a thread for every party 
      for (size_t pid = 1; pid <= nP_; pid++) {
        tpool_->enqueue([&, pid]() {
          ZZ_p_ctx.restore();
          // Count = 0 for multiplication triples // Count = 1 for MAC generation
          for (size_t count=0; count < 2; count++) {
            std::vector<Field> sharesVec;
            {
              // Create a lock for each OT 
              std::unique_lock<std::mutex> lock(mtx_);
              cv_start_ot_[count].wait(lock, [&]() { return start_ot_[count]; }); 
            }            
            sharesVec.resize(inputToOPE[count].size());
            chunk_ot_dig_pid_vec.resize(std::ceil(inputToOPE[count].size() / chunk_size_)); 
            // Do the OPEs in chunks of chunk_size 
            for (size_t start = 0; start < inputToOPE[count].size(); start += chunk_size_) {
              fieldDig chunk_ot_dig;
              size_t end = std::min(start + chunk_size_, inputToOPE[count].size());
              std::vector<Field> chunk(inputToOPE[count].begin() + start, inputToOPE[count].begin() + end);
              auto chunk_output = ot_[pid - 1]->multiplyRecv(chunk, chunk_ot_dig);
              std::copy(chunk_output.begin(), chunk_output.end(), sharesVec.begin() + start);
              chunk_ot_dig_pid_vec.push_back(std::make_pair(chunk_ot_dig, pid));
            }
            {
              std::lock_guard<std::mutex> lock(mtx_); 
              offline_message_buffer_[count].push({pid, sharesVec});
            }
            cv_.notify_one();
          } 
        });
      }
    } 
  } 

  void OfflineEvaluator::setupSync() {
    tpool_ = std::make_shared<ThreadPool>(1);
    if (id_ == 0) {
      // Create OT instance for each Part
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

    if (id_ == 0) {
      static ZZ_pContext ZZ_p_ctx;
      ZZ_p_ctx.save();

      // Start a thread only for sender party 
      auto pid = SYNC_SENDER_PID_;

      tpool_->enqueue([&, pid]() {
        ZZ_p_ctx.restore();
        // Count = 0 for multiplication triples // Count = 1 for MAC generation
        for (size_t count=0; count < 2; count++) {
          std::vector<Field> sharesVec;
          {
            // Create a lock for each OT 
            std::unique_lock<std::mutex> lock(mtx_);
            cv_start_ot_[count].wait(lock, [&]() { return start_ot_[count]; }); 
          }            
          sharesVec.resize(inputToOPE[count].size());
          chunk_ot_dig_pid_vec.resize(std::ceil(inputToOPE[count].size() / chunk_size_));  // TODO: make sure number of chunks is allocated correctly 
          // Do the OPEs in chunks of chunk_size 
          for (size_t start = 0; start < inputToOPE[count].size(); start += chunk_size_) {
            size_t end = std::min(start + chunk_size_, inputToOPE[count].size());
            std::vector<Field> chunk(inputToOPE[count].begin() + start, inputToOPE[count].begin() + end);
            fieldDig chunk_ot_dig;
            auto chunk_output = ot_[0]->multiplyRecv(chunk, chunk_ot_dig);
            chunk_ot_dig_pid_vec.push_back(std::make_pair(chunk_ot_dig, pid));
            std::copy(chunk_output.begin(), chunk_output.end(), sharesVec.begin() + start);
          }
          {
            std::lock_guard<std::mutex> lock(mtx_); 
            offline_message_buffer_[count].push({pid, sharesVec});
          }

          cv_.notify_one();
        } 
      });
    } 
  }

  OfflineEvaluator::~OfflineEvaluator() {
    tpool_.reset();
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

  // bool OfflineEvaluator::digestCheck(fieldDig ot_dig) {
  //   constexpr size_t honest_abort_message_length = 5;
  //     auto buf = std::vector<Field>(honest_abort_message_length); 
  //     network_->recv(0, buf.data(), buf.size() * sizeof(Field));
  //     auto hp_receiver_id = buf[honest_abort_message_length - 1];
        
  //     // If someone elses message was used for the OPE and their digest does not match mine, identify them as a cheater
  //     ot_dig = std::vector<Field>{Field(0), Field(0), Field(0), Field(0)}; // TODO: remove dummy input when ot_dig generation works properly 
  //     if(id_ != hp_receiver_id && !std::equal(buf.begin(), buf.end() - 1, ot_dig.begin())) {
  //       std::cout << "Party " << id_ << " has identified Party" << hp_receiver_id << " as a cheater!" << std::endl;
  //     }
  // }

  // bool OfflineEvaluator::sendDigest(fieldDig ot_dig) {
  //   constexpr size_t honest_abort_message_length = 5;
  //     auto buf = std::vector<Field>(honest_abort_message_length); 
  //     network_->recv(0, buf.data(), buf.size() * sizeof(Field));
  //     auto hp_receiver_id = buf[honest_abort_message_length - 1];
        
  //     // If someone elses message was used for the OPE and their digest does not match mine, identify them as a cheater
  //     ot_dig = std::vector<Field>{Field(0), Field(0), Field(0), Field(0)}; // TODO: remove dummy input when ot_dig generation works properly 
  //     if(id_ != hp_receiver_id && !std::equal(buf.begin(), buf.end() - 1, ot_dig.begin())) {
  //       std::cout << "Party " << id_ << " has identified Party" << hp_receiver_id << " as a cheater!" << std::endl;
  //     }
  // }

  void OfflineEvaluator::runOPE(std::vector<Field>& inputToOPE, std::vector <Field>& outputOfOPE, size_t count) {
    run_async_ ? runOPEASync(inputToOPE, outputOfOPE, count) : runOPESync(inputToOPE, outputOfOPE, count);
  }

  void OfflineEvaluator::runOPEASync(std::vector<Field>& inputToOPE, std::vector <Field>& outputOfOPE, size_t count) {
    constexpr size_t honest_abort_message_length = 5;
    if(id_ != 0) {
      for (size_t start = 0; start < inputToOPE.size(); start += chunk_size_) {
          // Complete OPE and compute digest of my message to HP 
          size_t end = std::min(start + chunk_size_, inputToOPE.size());
          std::vector<Field> chunk(inputToOPE.begin() + start, inputToOPE.begin() + end);
          fieldDig ot_dig;
          auto chunk_output = ot_[0]->multiplySend(chunk, rgen_.all_minus_0(), ot_dig);
          outputOfOPE.insert(outputOfOPE.end(), chunk_output.begin(), chunk_output.end());
      }
    } else {
      // HP completes OPE with parties by proceding with first message it receives  
      {
        std::lock_guard<std::mutex> lock(mtx_);
        start_ot_[count] = true;
      }
      cv_start_ot_[count].notify_one();
      {
        std::unique_lock<std::mutex> lock(mtx_);
        cv_.wait(lock, [&]() { return offline_message_buffer_[count].size() >= 1; });
      }
      Offline_Message OPE_res = offline_message_buffer_[count].front();
      auto receiver_pid = OPE_res.receiver_id;
      std::queue<Offline_Message> empty;
      std::swap(offline_message_buffer_[count], empty);
      outputOfOPE = OPE_res.data;
    }
  }


  void OfflineEvaluator::runOPESync(std::vector<Field>& inputToOPE, std::vector <Field>& outputOfOPE, size_t count) {
    constexpr size_t honest_abort_message_length = 4;
    if(id_ != 0) {
      for (size_t start = 0; start < inputToOPE.size(); start += chunk_size_) {
          size_t end = std::min(start + chunk_size_, inputToOPE.size());
          std::vector<Field> chunk(inputToOPE.begin() + start, inputToOPE.begin() + end);
          fieldDig ot_dig;
          // Complete OPE and compute digest of my message to HP 
          if(id_ == SYNC_SENDER_PID_) {
            auto chunk_output = ot_[0]->multiplySend(chunk, rgen_.all_minus_0(), ot_dig);
            outputOfOPE.insert(outputOfOPE.end(), chunk_output.begin(), chunk_output.end());
          } else {
            auto chunk_output = ot_[0]->multiplySendOffline(chunk, rgen_.all_minus_0(), ot_dig);
            outputOfOPE.insert(outputOfOPE.end(), chunk_output.begin(), chunk_output.end());
          }
      }
    } 
    else {
      // HP completes OPE with parties by proceding with first message it receives  
      {
        std::lock_guard<std::mutex> lock(mtx_);
        start_ot_[count] = true;
      }
      cv_start_ot_[count].notify_one();
      {
        std::unique_lock<std::mutex> lock(mtx_);
        cv_.wait(lock, [&]() { return offline_message_buffer_[count].size() >= 1; });
      }
      Offline_Message OPE_res = offline_message_buffer_[count].front();
      auto receiver_pid = OPE_res.receiver_id;
      std::queue<Offline_Message> empty;
      std::swap(offline_message_buffer_[count], empty);
      outputOfOPE = OPE_res.data;
    }
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
            if (id_ == nP_) {
              buffer_num += 2;
            }
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

    if (id_ != nP_) {
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
