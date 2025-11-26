#include "dmgod_offline_evaluator.h"

namespace dmAsyncAsteriskGOD {
  OfflineEvaluator::OfflineEvaluator(int nP, int id, int security_param, std::shared_ptr<NetIOMP> network1, 
    std::shared_ptr<NetIOMP> network2, LevelOrderedCircuit circ, int threads, uint64_t seed) 
    : nP_(nP), id_(id), security_param_(security_param), rgen_(id, nP, seed), network_(std::move(network1)), 
    network_ot_(std::move(network2)), circ_(std::move(circ)), preproc_(circ.num_gates), start_ot_(2, false), 
    chunk_size_(50000), inputToOPE(2), run_async_(false)
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
            chunk_dig_pid_.resize(std::ceil(inputToOPE[count].size() / chunk_size_)); 
            // Do the OPEs in chunks of chunk_size 
            for (size_t start = 0; start < inputToOPE[count].size(); start += chunk_size_) {
              fieldDig chunk_ot_dig;
              size_t end = std::min(start + chunk_size_, inputToOPE[count].size());
              std::vector<Field> chunk(inputToOPE[count].begin() + start, inputToOPE[count].begin() + end);
              auto chunk_output = ot_[pid - 1]->multiplyRecv(chunk, chunk_ot_dig);
              std::copy(chunk_output.begin(), chunk_output.end(), sharesVec.begin() + start);
              chunk_dig_pid_.push_back(std::make_pair(chunk_ot_dig, pid));
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
    tpool_ = std::make_shared<ThreadPool>(1); // Threadpool for OPE with one party 
    tpool_ope_sync_ = std::make_shared<ThreadPool>(nP_ - 1); // Threadpool to send OPE digest back to other parties for honest abort 
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
          chunk_dig_pid_.resize(std::ceil(inputToOPE[count].size() / chunk_size_));  // TODO: make sure number of chunks is allocated correctly 
          // Do the OPEs in chunks of chunk_size 
          for (size_t start = 0; start < inputToOPE[count].size(); start += chunk_size_) {
            size_t end = std::min(start + chunk_size_, inputToOPE[count].size());
            std::vector<Field> chunk(inputToOPE[count].begin() + start, inputToOPE[count].begin() + end);
            fieldDig chunk_dig;
            auto chunk_output = ot_[0]->multiplyRecv(chunk, chunk_dig);
            chunk_dig_pid_.push_back(std::make_pair(chunk_dig, pid));
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

  bool OfflineEvaluator::verifyOPEMsgsASync() {
    constexpr size_t kChunkMsgLength = 5;   // 1 sender_id + 4 digest fields
    constexpr size_t kDigestLength   = 4;
    const size_t num_chunks = chunk_dig_pid_.size();
    const size_t total_comm = kChunkMsgLength * num_chunks;

    /**
     * recv_buf/send_buf layout:
     *   [ sender_id | digest(4) ] ... (for each chunk)
     */

    if(id_ == 0) {
      std::vector<Field> send_buf(total_comm);

      // Prepare send_buf 
      size_t idx = 0;
      for (auto& [chunk_dig, sender_pid] : chunk_dig_pid_) {
          send_buf[idx++] = Field(sender_pid);
          for (auto& dig_elem : chunk_dig) {
            send_buf[idx++] = dig_elem;
          }
      }

      // Send send_buf to each party 
      std::vector<std::future<void>> send_t;    
      for(size_t pid = 1; pid <= nP_; pid++) {
          send_t.push_back(tpool_->enqueue([&,pid]() {
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
      // Receive buffer
      std::vector<Field> recv_buf(total_comm);
      network_->recv(0, recv_buf.data(), recv_buf.size() * sizeof(Field));

      for (size_t i = 0; i < num_chunks; i++) {
        const size_t offset = i * kChunkMsgLength;

        // Extract sender from buffer 
        const Field sender_id = recv_buf[offset];

        // Extract received digest from buffer
        fieldDig recv_digest(
            recv_buf.begin() + offset + 1,
            recv_buf.begin() + offset + 1 + kDigestLength
        );

        // My digest for chunk i
        const fieldDig& my_digest = chunk_dig_pid_[i].first;

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

  bool OfflineEvaluator::verifyOPEMsgsSync() {
    constexpr size_t kChunkMsgLength = 4;   // 1 sender_id + 4 digest fields
    constexpr size_t kDigestLength   = 4;
    const size_t num_chunks = chunk_dig_pid_.size();
    const size_t total_comm = kChunkMsgLength * num_chunks;
    /**
     * recv_buf/send_buf layout:
     *   [ sender_id | digest(4) ] ... (for each chunk)
     */

    if(id_ == 0) {
      std::vector<Field> send_buf(total_comm);

      // Prepare send_buf 
      size_t idx = 0;
      for (auto& [chunk_dig, sender_pid] : chunk_dig_pid_) {
          for (auto& dig_elem : chunk_dig)
              send_buf[idx++] = dig_elem;
      }

      // Send send_buf to each party 
      std::vector<std::future<void>> send_t;    
      for(size_t pid = 1; pid <= nP_; pid++) {
          if(pid != SYNC_SENDER_PID_) { // Don't send to designated sender in OPE 
            send_t.push_back(tpool_ope_sync_->enqueue([&,pid]() {
                network_->send(pid, send_buf.data(), sizeof(Field) * total_comm);
                network_->getSendChannel(pid)->flush();
            }));
          }
      }
      for(auto& t : send_t) {
          if (t.valid()) {
              t.wait();
          }
      }
    }
    else if(id_ != SYNC_SENDER_PID_) {
      // Receive buffer
      std::vector<Field> recv_buf(total_comm);
      network_->recv(0, recv_buf.data(), recv_buf.size() * sizeof(Field));

      for (size_t i = 0; i < num_chunks; i++) {
        const size_t offset = i * kChunkMsgLength;

        // Extract received digest from buffer
        fieldDig recv_digest(
            recv_buf.begin() + offset,
            recv_buf.begin() + offset + kDigestLength
        );

        // My digest for chunk i
        const fieldDig& my_digest = chunk_dig_pid_[i].first;

        // Check mismatch
        if (my_digest != recv_digest) {
            std::cout << "Party " << id_
                      << " has identified Party " << SYNC_SENDER_PID_
                      << " as a cheater!" << std::endl;
            return false;
        }
      }
    }
    return true;
  }

  void OfflineEvaluator::runOPE(std::vector<Field>& inputToOPE, std::vector<Field>& outputOfOPE, size_t count) {
    run_async_ ? runOPEASync(inputToOPE, outputOfOPE, count) : runOPESync(inputToOPE, outputOfOPE, count);
  }

  void OfflineEvaluator::runOPEASync(std::vector<Field>& inputToOPE, std::vector <Field>& outputOfOPE, size_t count) {
    if(id_ != 0) {
      chunk_dig_pid_.resize(std::ceil(inputToOPE.size() / chunk_size_));
      for (size_t start = 0; start < inputToOPE.size(); start += chunk_size_) {
          // Complete OPE and compute digest of my message to HP 
          size_t end = std::min(start + chunk_size_, inputToOPE.size());
          std::vector<Field> chunk(inputToOPE.begin() + start, inputToOPE.begin() + end);
          fieldDig chunk_dig;
          auto chunk_output = ot_[0]->multiplySend(chunk, rgen_.all_minus_0(), chunk_dig);
          outputOfOPE.insert(outputOfOPE.end(), chunk_output.begin(), chunk_output.end());
          chunk_dig_pid_.push_back(std::make_pair(chunk_dig, id_));
      }
    } else {
      // HP completes OPE with parties by proceding with first message it receives  
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
    }

    verifyOPEMsgsASync();
  }


  void OfflineEvaluator::runOPESync(std::vector<Field>& inputToOPE, std::vector <Field>& outputOfOPE, size_t count) {
    constexpr size_t honest_abort_message_length = 4;

    if(id_ != 0) {
      chunk_dig_pid_.resize(std::ceil(inputToOPE.size() / chunk_size_));
      for (size_t start = 0; start < inputToOPE.size(); start += chunk_size_) {
          size_t end = std::min(start + chunk_size_, inputToOPE.size());
          std::vector<Field> chunk(inputToOPE.begin() + start, inputToOPE.begin() + end);
          fieldDig chunk_dig;
          // Complete OPE and compute digest of my message to HP 
          if(id_ == SYNC_SENDER_PID_) {
            auto chunk_output = ot_[0]->multiplySend(chunk, rgen_.all_minus_0(), chunk_dig);
            outputOfOPE.insert(outputOfOPE.end(), chunk_output.begin(), chunk_output.end());
          } else {
            auto chunk_output = ot_[0]->multiplySendOffline(chunk, rgen_.all_minus_0(), chunk_dig);
            outputOfOPE.insert(outputOfOPE.end(), chunk_output.begin(), chunk_output.end());
          }
          chunk_dig_pid_.push_back(std::make_pair(chunk_dig, id_));
      }
    } 
    else {
      // HP completes OPE with parties by proceding with first message it receives  
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
      std::queue<Offline_Message> empty;
      std::swap(offline_message_buffer_[count], empty);
      outputOfOPE = OPE_res.data;
    }

    verifyOPEMsgsSync();
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