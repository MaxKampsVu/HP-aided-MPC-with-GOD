#include "dm_offline_evaluator.h"

namespace dmAsyncAsterisk {
  OfflineEvaluator::OfflineEvaluator(int nP, int id, int security_param, std::shared_ptr<NetIOMP> network1, 
    std::shared_ptr<NetIOMP> network2, LevelOrderedCircuit circ, int threads, uint64_t seed) 
    : nP_(nP), id_(id), security_param_(security_param), rgen_(id, nP, seed), network_(std::move(network1)), 
    network_ot_(std::move(network2)), circ_(std::move(circ)), preproc_(circ.num_gates), start_ot_(2, false), 
    chunk_size_(50000), inputToOPE(2)  
  {
    tpool_ = std::make_shared<ThreadPool>(threads);
    // HP does the following 
    if (id_ == 0) {
      // Creare recv channel for every Party
      for (size_t i = 1; i <= nP_; i++) {
        ot_.emplace_back(std::make_unique<OTProvider>(id_, i, network_ot_->getRecvChannel(i)));
        network_ot_->getRecvChannel(i)->flush();
      }
    }
    // Parties do the following 
    else {
      // Create send send channel to HP 
      ot_.emplace_back(std::make_unique<OTProvider>(id_, 0, network_ot_->getSendChannel(0)));
      network_ot_->getSendChannel(0)->flush();
    }
    
    // HP does the following: 
    if (id_ == 0) {
      static ZZ_pContext ZZ_p_ctx;
      ZZ_p_ctx.save();
      // Start a thread for every party 
      for (size_t pid = 1; pid <= nP_; pid++) {
        tpool_->enqueue([&, pid]() {
          ZZ_p_ctx.restore();
          // Start an OLE for all pairs of shares (x0, y0) for multiplication gates 
          for (size_t count=0; count < 1; count++) {
            std::vector<Field> sharesVec;
            {
              // Create a lock for each OT 
              std::unique_lock<std::mutex> lock(mtx_);
              // Start the ot for x0 (and in the next iteration for y0)
              cv_start_ot_[count].wait(lock, [&]() { return start_ot_[count]; });
            }            
            sharesVec.resize(inputToOPE[count].size()); 
            // Do the OLEs in chunks of chunk_size 
            for (size_t start = 0; start < inputToOPE[count].size(); start += chunk_size_) {
              size_t end = std::min(start + chunk_size_, inputToOPE[count].size());
              std::vector<Field> chunk(inputToOPE[count].begin() + start, inputToOPE[count].begin() + end);
              // initiate ot, with ot provider pid and multiply my chunk with whatever ot provider sends 
              std::vector<Field> chunk_output = ot_[pid - 1]->multiplyRecv(chunk);
              // Copy chunk_output to sharesVec 
              std::copy(chunk_output.begin(), chunk_output.end(), sharesVec.begin() + start);
            }
            {
              std::lock_guard<std::mutex> lock(mtx_);
              // Record with whomst sharesVec was computed 
              offline_message_buffer_[count].push({pid, sharesVec});
            }
            // TODO: Look up what this does 
            cv_.notify_one();
          }
          //TODO: What does this do? 
          for(size_t count=0; count < 3; count++) {
            size_t total_comm;
            network_->recv(pid, &total_comm, sizeof(size_t));
            std::vector<Field> offline_comm_to_HP(total_comm);
            network_->recv(pid, offline_comm_to_HP.data(), sizeof(Field) * total_comm);
            {
                std::lock_guard<std::mutex> lock(mtx_);
                offline_message_buffer_[count+2].push({pid, offline_comm_to_HP});
            }
            cv_.notify_one();
          }
        });
      }
    }
  }

  OfflineEvaluator::~OfflineEvaluator() {
    tpool_.reset();
  }

  PreprocCircuit<Field> OfflineEvaluator::dummy(int nP, int id, const LevelOrderedCircuit& circ, 
    const std::unordered_map<wire_t, int>& input_pid_map, PRG& prg) {
    PreprocCircuit<Field> preproc(circ.num_gates);
    DummyShare<Field> mask, mask_prod;

    std::vector<Field> key_sh(nP+1), values(nP+1), tags(nP+1);
    Field key = Field(0), value = Field(0), tag = Field(0);
    for (size_t i=0; i<nP+1; i++) {
      randomizeZZp(prg, key_sh[i], sizeof(Field));
      key += key_sh[i];
    }
    mask.setKeySh(key_sh);
    mask.setKey(key);

    for (size_t i=0; i<nP+1; i++) {
      randomizeZZp(prg, values[i], sizeof(Field));
      value += values[i];
    }
    mask.setValues(values);
    mask.setValue(value);

    tag = value * key;
    mask.setTag(tag);
    for (size_t i=0; i<nP; i++) {
      randomizeZZp(prg, tags[i], sizeof(Field));
      tag -= tags[i];
    }
    tags[nP] = tag;
    mask.setTags(tags);

    mask_prod.setKeySh(key_sh);
    mask_prod.setKey(key);

    Field secret = value * value;
    tag = secret * key;    
    mask_prod.setValue(secret);
    mask_prod.setTag(tag);
    for (size_t i=0; i<nP; i++) {
      randomizeZZp(prg, values[i], sizeof(Field));
      randomizeZZp(prg, tags[i], sizeof(Field));
      secret -= values[i];
      tag -= tags[i];
    }
    values[nP] = secret;
    tags[nP] = tag;
    mask_prod.setValues(values);
    mask_prod.setTags(tags);

    for (const auto& level : circ.gates_by_level) {
      for (const auto& gate : level) {
        switch (gate->type) {
          case GateType::kInp: {
            auto pregate = std::make_unique<PreprocInput<Field>>();      
            auto pid = input_pid_map.at(gate->out);
            pregate->pid = pid;
            if (id == pid) {
              pregate->mask_value = mask.secretVal();
            }
            pregate->mask = mask.getRSS(id);
            preproc.gates[gate->out] = std::move(pregate);                
            break;
          }

          case GateType::kAdd: {
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            DummyShare<Field> temp = mask + mask;
            preproc.gates[gate->out] = std::make_unique<PreprocGate<Field>>(temp.getRSS(id));    
            break;
          }
  
          case GateType::kConstAdd: {
            const auto* g = static_cast<ConstOpGate<Field>*>(gate.get());
            preproc.gates[gate->out] = std::make_unique<PreprocGate<Field>>(mask.getRSS(id));
            break;
          }
  
          case GateType::kConstMul: {
            const auto* g = static_cast<ConstOpGate<Field>*>(gate.get());
            DummyShare<Field> temp = mask * g->cval;
            preproc.gates[gate->out] = std::make_unique<PreprocGate<Field>>(temp.getRSS(id));
            break;
          }
  
          case GateType::kSub: {
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            DummyShare<Field> temp = mask - mask;
            preproc.gates[gate->out] = std::make_unique<PreprocGate<Field>>(temp.getRSS(id));
            break;
          }

          case GateType::kMul: {
            preproc.gates[gate->out] = std::make_unique<PreprocMultGate<Field>>();
            preproc.gates[gate->out] = std::move(std::make_unique<PreprocMultGate<Field>> (mask.getRSS(id), mask_prod.getRSS(id)));
            break;
          }

          default: {
            break;
          }
        }
      }
    }
    return preproc;
  }

  void OfflineEvaluator::keyGen()  {
    if(id_ == 0) {
      key_sh_.resize(nP_);
      for(size_t i = 1; i <= nP_; i++) {
        randomizeZZp(rgen_.pi(i), key_sh_[i-1], sizeof(Field));
      }
    }
    else {
      key_sh_.resize(2);
      randomizeZZp(rgen_.p0(), key_sh_[0], sizeof(Field));
      randomizeZZp(rgen_.all_minus_0(), key_sh_[1], sizeof(Field));
    }
  }

  void OfflineEvaluator::randomShare(int nP, int pid, RandGenPool& rgen, RepShare<Field>& share, const std::vector<Field>& keySh, 
    Field& mask_share_zero, bool isOutputWire) {
    if(pid == 0) {      
      for(size_t i = 1; i <= nP; i++) {
        Field valSh;
        randomizeZZp(rgen.pi(i), valSh, sizeof(Field));
        share.pushValues(valSh);
        share.setKeySh(keySh[i-1]);
      }
      if (isOutputWire)
        randomizeZZp(rgen.all(), mask_share_zero, sizeof(Field));
    }
    else {
      Field valSh;
      randomizeZZp(rgen.p0(), valSh, sizeof(Field));
      share.pushValues(valSh);
      share.setKeySh(keySh[0]);
      if (isOutputWire)
        randomizeZZp(rgen.all(), valSh, sizeof(Field));
      else
        randomizeZZp(rgen.all_minus_0(), valSh, sizeof(Field));
      share.pushValues(valSh);
      share.setKeySh(keySh[1]);
    }
  }

  // Start OLE to compute share1.p0 * sum(share1.p1, ..., share1.pn) + share1.p0 * sum(share1.p1, ..., share1.pn)
  void OfflineEvaluator::randomShareSecret(int nP, int pid, RandGenPool& rgen, const RepShare<Field>& share1, const RepShare<Field>& share2, 
    RepShare<Field>& prodShare, const std::vector<Field>& keySh, std::vector<Field>& inputToOPE) {
    auto share1_vals = share1.getValues();
    auto share2_vals = share2.getValues();
    if (pid!=0) {   
      // Parties participate with their common share in OPE    
      inputToOPE.push_back(share1_vals[1]);
      inputToOPE.push_back(share2_vals[1]);
    }
    // HP 
    else {
      // Sum values for OPE  
      Field a = Field(0), b = Field(0);
      for (size_t i=0; i < nP; i++) {
        a += share1_vals[i];
        b += share2_vals[i];
      }
      // Input for OPE is the sum of shares 
      inputToOPE.push_back(b);
      inputToOPE.push_back(a);
    }
    // prodShare will be initialized after OLE 
    std::vector<Field> vals(share1_vals.size());
    prodShare.setValues(vals);
    prodShare.setKeySh(keySh);
  }  

  void OfflineEvaluator::randomShareWithParty(int nP, int pid, int dealer, RandGenPool& rgen, RepShare<Field>& share, Field& secret, 
    const std::vector<Field>& keySh) {        
      secret = Field(0);
      Field valSh;
      if(pid == 0) {  
      if (pid != dealer){
        for(size_t i = 1; i <= nP; i++) {
          randomizeZZp(rgen.pij(i, dealer), valSh, sizeof(Field));
          share.pushValues(valSh);
          share.setKeySh(keySh[i-1]);
        }
      }
      else {
        for(size_t i = 1; i <= nP; i++) {
          randomizeZZp(rgen.pi(i), valSh, sizeof(Field));
          share.pushValues(valSh);
          secret += valSh;
          share.setKeySh(keySh[i-1]);
        }
        Field val;
        randomizeZZp(rgen.all(), val, sizeof(Field));
        secret += val;
      }
    }
    else {
      if (dealer == 0) {
        randomizeZZp(rgen.p0(), valSh, sizeof(Field));
        share.pushValues(valSh);
        share.setKeySh(keySh[0]);
        randomizeZZp(rgen.all(), valSh, sizeof(Field));
        share.pushValues(valSh);
        share.setKeySh(keySh[1]);
      }
      else if (pid != dealer) {
        randomizeZZp(rgen.pj(dealer), valSh, sizeof(Field));
        share.pushValues(valSh);
        share.setKeySh(keySh[0]);
        randomizeZZp(rgen.all_minus_0(), valSh, sizeof(Field));
        share.pushValues(valSh);
        share.setKeySh(keySh[1]);
      }
      else {
        for(size_t j = 1; j <= nP; j++) {
          randomizeZZp(rgen.pj(j), valSh, sizeof(Field));
          secret += valSh;
          if (pid == j)
            share.pushValues(valSh);
        }
        share.setKeySh(keySh[0]);
        randomizeZZp(rgen.all_minus_0(), valSh, sizeof(Field));
        secret += valSh;
        share.pushValues(valSh);
        share.setKeySh(keySh[1]);
      }
    }
  }

  void OfflineEvaluator::runOPE(std::vector<Field>& inputToOPE, std::vector <Field>& outputOfOPE, size_t count) {
    if (id_ != 0) {
        for (size_t start = 0; start < inputToOPE.size(); start += chunk_size_) {
            size_t end = std::min(start + chunk_size_, inputToOPE.size());
            std::vector<Field> chunk(inputToOPE.begin() + start, inputToOPE.begin() + end);
            std::vector<Field> chunk_output = ot_[0]->multiplySend(chunk, rgen_.all_minus_0());
            outputOfOPE.insert(outputOfOPE.end(), chunk_output.begin(), chunk_output.end());
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
      std::queue<Offline_Message> empty;
      std::swap(offline_message_buffer_[count], empty);
      outputOfOPE = OPE_res.data;
    }
  }

  void OfflineEvaluator::multiply(const std::vector<Field>& vec_a, const std::vector<Field>& vec_b, std::vector<Field>& vec_c, 
    const std::vector<Field>& outputOfOPE, std::vector<Field>& buffer, size_t& idx_outputOfOPE, size_t& idx_buffer) {
    Field buffer_elem = Field(0);
    // Multiply mt shares one by one a0 * b0 + a1 * b1 + ...
    for (size_t i=0; i < vec_a.size(); i++) {
      vec_c[i] = vec_a[i] * vec_b[i];
    }
    // HP does 
    if (id_==0) {
      Field t = Field(0);
      // Compute the cross terms I can compute e.g. a1*b2
      for (size_t i=0; i < nP_; i++) {
        for (size_t j=0; j < nP_; j++) {
          if (i!=j)
            t+= vec_a[i] * vec_b[j];
        }
      }
      for (size_t i=1; i < nP_; i++) {
        Field c;
        randomizeZZp(rgen_.pi(i), c, sizeof(Field));
        vec_c[i-1] += c;
        t -= c;
      }
      vec_c[nP_-1] += t;
      buffer_elem += t;
    }
    // Parties other then last party do 
    else if (id_<nP_) {
      Field c;
      randomizeZZp(rgen_.p0(), c, sizeof(Field));
      vec_c[0] += c;
    }

    // Hp does 
    if (id_==0) {
      Field t = outputOfOPE[idx_outputOfOPE] + outputOfOPE[idx_outputOfOPE+1];
      idx_outputOfOPE += 2;
      for (size_t i=1; i < nP_; i++) {
        Field c;
        randomizeZZp(rgen_.pi(i), c, sizeof(Field));
        vec_c[i-1] += c;
        t -= c;
      }
      vec_c[nP_-1] += t;
      buffer_elem += t;
    }
    // Parties other then HP do 
    else if (id_<nP_) {
      vec_c[1] += outputOfOPE[idx_outputOfOPE] + outputOfOPE[idx_outputOfOPE+1];
      idx_outputOfOPE += 2;
      Field c;
      randomizeZZp(rgen_.p0(), c, sizeof(Field));
      vec_c[0] += c;
    }
    // Last party does 
    else {
      vec_c[1] += outputOfOPE[idx_outputOfOPE] + outputOfOPE[idx_outputOfOPE+1];
      idx_outputOfOPE += 2;
    }
    //HP does 
    if (id_==0)
      buffer.push_back(buffer_elem);
    //other parties do 
    else if (id_==nP_)
      vec_c[0] += buffer[idx_buffer++];
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
            randomShareWithParty(nP_, id_, pid, rgen_, pregate->mask, pregate->mask_value, key_sh_);      
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
            preproc_.gates[gate->out] = std::make_unique<PreprocMultGate<Field>>();
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            const auto& mask_in1 = preproc_.gates[g->in1]->mask;
            const auto& mask_in2 = preproc_.gates[g->in2]->mask;
            RepShare<Field> rand_mask, rand_ver;
            Field mask_share_zero = Field(0);
            bool isOutputWire = false;
            if (std::find(circ_.outputs.begin(), circ_.outputs.end(),g->out)!=circ_.outputs.end())
              isOutputWire = true;
            randomShare(nP_, id_, rgen_, rand_mask, key_sh_, mask_share_zero, isOutputWire);
            randomShare(nP_, id_, rgen_, rand_ver, key_sh_, mask_share_zero, false);
            RepShare<Field> mask_product, ver_product;
            // Create mask
            randomShareSecret(nP_, id_, rgen_, mask_in1, mask_in2, mask_product, key_sh_, inputToOPE[0]);
            // TODO: what is rand_ver 
            randomShareSecret(nP_, id_, rgen_, rand_ver, mask_in2, ver_product, key_sh_, inputToOPE[0]);
            preproc_.gates[gate->out] = std::move(std::make_unique<PreprocMultGate<Field>> (rand_mask, mask_product, mask_share_zero, rand_ver, ver_product));
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
    
    // Run the OPEs to obtain the product shares for the multiplication gates 
    runOPE(inputToOPE[0], outputOfOPE, 0);

    if (id_ != nP_) {
      for (const auto& level : circ_.gates_by_level) {
        for (const auto& gate : level) {
          switch (gate->type) {
            case GateType::kMul: {
              auto* g = static_cast<FIn2Gate*>(gate.get());
              auto *pre_mul = static_cast<PreprocMultGate<Field> *>(preproc_.gates[gate->out].get());
              auto valVec = pre_mul->mask_prod.getValues();
              auto mask_in1_val = preproc_.gates[g->in1]->mask.getValues();
              auto mask_in2_val = preproc_.gates[g->in2]->mask.getValues();
              // Compute the mask on the ouput wire of the multiplication gate 
              multiply(mask_in1_val, mask_in2_val, valVec, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
              pre_mul->mask_prod.setValues(valVec);
              valVec = pre_mul->ver_prod.getValues();
              auto share1 = pre_mul->ver.getValues();
              auto share2 = preproc_.gates[g->in2]->mask.getValues();
              multiply(share1, share2, valVec, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
              pre_mul->ver_prod.setValues(valVec);
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

      if (id_ == 0) {
        buffer_num = buffer.size();
        network_->send(nP_, buffer.data(), sizeof(Field) * buffer_num);
        network_->getSendChannel(nP_)->flush();        
        buffer.clear();
        buffer.shrink_to_fit();
      }
    }
    else {
      buffer.resize(buffer_num);
      network_->recv(0, buffer.data(), sizeof(Field) * buffer_num);

      for (const auto& level : circ_.gates_by_level) {
        for (const auto& gate : level) {
          switch (gate->type) {
            case GateType::kMul: {
              auto* g = static_cast<FIn2Gate*>(gate.get());
              auto *pre_mul = static_cast<PreprocMultGate<Field> *>(preproc_.gates[gate->out].get());
              auto valVec = pre_mul->mask_prod.getValues();
              auto mask_in1_val = preproc_.gates[g->in1]->mask.getValues();
              auto mask_in2_val = preproc_.gates[g->in2]->mask.getValues();
              multiply(mask_in1_val, mask_in2_val, valVec, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
              pre_mul->mask_prod.setValues(valVec);
              valVec = pre_mul->ver_prod.getValues();
              auto share1 = pre_mul->ver.getValues();
              auto share2 = preproc_.gates[g->in2]->mask.getValues();
              multiply(share1, share2, valVec, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
              pre_mul->ver_prod.setValues(valVec);
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
      buffer.clear();
      buffer.shrink_to_fit();
    }
  }

  void OfflineEvaluator::prepareMaskTags() {
    std::vector<Field> buffer; 
    size_t idx_buffer=0; 
    size_t buffer_num=0;

    for (const auto& level : circ_.gates_by_level) {
      for (const auto& gate : level) {
        switch (gate->type) {
          case GateType::kInp: {
            if (id_ == nP_) {
              buffer_num += 1;
            }
            auto *pre_input = static_cast<PreprocInput<Field> *>(preproc_.gates[gate->out].get());
            auto share1 = pre_input->mask.getValues();
            auto share2 = pre_input->mask.getKeySh();
            if (id_!=0) {      
              inputToOPE[1].push_back(share1[1]);
              inputToOPE[1].push_back(share2[1]);
            }
            else {
              Field a = Field(0), b = Field(0);
              for (size_t i=0; i < nP_; i++) {
                a += share1[i];
                b += share2[i];
              }
              inputToOPE[1].push_back(b);
              inputToOPE[1].push_back(a);
            }
            break;
          }

          case GateType::kMul: {
            if (id_ == nP_) {
              buffer_num += 4;
            }
            auto *pre_mul = static_cast<PreprocMultGate<Field> *>(preproc_.gates[gate->out].get());
            std::vector<Field> share1, share2;
            share1 = pre_mul->mask.getValues();
            share2 = pre_mul->mask.getKeySh();
            if (id_!=0) {      
              inputToOPE[1].push_back(share1[1]);
              inputToOPE[1].push_back(share2[1]);
            }
            else {
              Field a = Field(0), b = Field(0);
              for (size_t i=0; i < nP_; i++) {
                a += share1[i];
                b += share2[i];
              }
              inputToOPE[1].push_back(b);
              inputToOPE[1].push_back(a);
            }
            share1 = pre_mul->mask_prod.getValues();
            share2 = pre_mul->mask_prod.getKeySh();
            if (id_!=0) {      
              inputToOPE[1].push_back(share1[1]);
              inputToOPE[1].push_back(share2[1]);
            }
            else {
              Field a = Field(0), b = Field(0);
              for (size_t i=0; i < nP_; i++) {
                a += share1[i];
                b += share2[i];
              }
              inputToOPE[1].push_back(b);
              inputToOPE[1].push_back(a);
            }
            share1 = pre_mul->ver.getValues();
            share2 = pre_mul->ver.getKeySh();
            if (id_!=0) {      
              inputToOPE[1].push_back(share1[1]);
              inputToOPE[1].push_back(share2[1]);
            }
            else {
              Field a = Field(0), b = Field(0);
              for (size_t i=0; i < nP_; i++) {
                a += share1[i];
                b += share2[i];
              }
              inputToOPE[1].push_back(b);
              inputToOPE[1].push_back(a);
            }
            share1 = pre_mul->ver_prod.getValues();
            share2 = pre_mul->ver_prod.getKeySh();
            if (id_!=0) {      
              inputToOPE[1].push_back(share1[1]);
              inputToOPE[1].push_back(share2[1]);
            }
            else {
              Field a = Field(0), b = Field(0);
              for (size_t i=0; i < nP_; i++) {
                a += share1[i];
                b += share2[i];
              }
              inputToOPE[1].push_back(b);
              inputToOPE[1].push_back(a);
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
              std::vector<Field> tags;
              auto share1 = pre_input->mask.getValues();;
              auto share2 = pre_input->mask.getKeySh();
              tags.resize(share1.size());
              multiply(share1, share2, tags, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
              pre_input->mask.setTags(tags);                
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
              auto *pre_mul = static_cast<PreprocMultGate<Field> *>(preproc_.gates[gate->out].get());
              std::vector<Field> share1, share2, tags1, tags2, tags3, tags4;
              share1 = pre_mul->mask.getValues();
              share2 = pre_mul->mask.getKeySh();
              tags1.resize(share1.size());
              multiply(share1, share2, tags1, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
              pre_mul->mask.setTags(tags1);
              share1 = pre_mul->mask_prod.getValues();
              share2 = pre_mul->mask_prod.getKeySh();
              tags2.resize(share1.size());
              multiply(share1, share2, tags2, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
              pre_mul->mask_prod.setTags(tags2);
              share1 = pre_mul->ver.getValues();
              share2 = pre_mul->ver.getKeySh();
              tags3.resize(share1.size());
              multiply(share1, share2, tags3, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
              pre_mul->ver.setTags(tags3);
              share1 = pre_mul->ver_prod.getValues();
              share2 = pre_mul->ver_prod.getKeySh();
              tags4.resize(share1.size());
              multiply(share1, share2, tags4, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
              pre_mul->ver_prod.setTags(tags4);            
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

      if (id_ == 0) {
        buffer_num = buffer.size();
        network_->send(nP_, buffer.data(), sizeof(Field) * buffer_num);
        network_->getSendChannel(nP_)->flush();
        buffer.clear();
        buffer.shrink_to_fit();
      }      
    }
    else {
      buffer.resize(buffer_num);
      network_->recv(0, buffer.data(), sizeof(Field) * buffer_num);

      for (const auto& level : circ_.gates_by_level) {
        for (const auto& gate : level) {
          switch (gate->type) {
            case GateType::kInp: {
              auto *pre_input = static_cast<PreprocInput<Field> *>(preproc_.gates[gate->out].get());
              std::vector<Field> tags;
              auto share1 = pre_input->mask.getValues();
              auto share2 = pre_input->mask.getKeySh();
              tags.resize(share1.size());
              multiply(share1, share2, tags, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
              pre_input->mask.setTags(tags);                
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
              auto *pre_mul = static_cast<PreprocMultGate<Field> *>(preproc_.gates[gate->out].get());
              std::vector<Field> share1, share2, tags1, tags2, tags3, tags4;
              share1 = pre_mul->mask.getValues();
              share2 = pre_mul->mask.getKeySh();
              tags1.resize(share1.size());
              multiply(share1, share2, tags1, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
              pre_mul->mask.setTags(tags1);
              share1 = pre_mul->mask_prod.getValues();
              share2 = pre_mul->mask_prod.getKeySh();
              tags2.resize(share1.size());
              multiply(share1, share2, tags2, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
              pre_mul->mask_prod.setTags(tags2);
              share1 = pre_mul->ver.getValues();
              share2 = pre_mul->ver.getKeySh();
              tags3.resize(share1.size());
              multiply(share1, share2, tags3, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
              pre_mul->ver.setTags(tags3);
              share1 = pre_mul->ver_prod.getValues();
              share2 = pre_mul->ver_prod.getKeySh();
              tags4.resize(share1.size());
              multiply(share1, share2, tags4, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
              pre_mul->ver_prod.setTags(tags4);
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
      buffer.clear();
      buffer.shrink_to_fit();
    }
  }
    
  void OfflineEvaluator::setWireMasks(const std::unordered_map<wire_t,int>& input_pid_map) {      
    keyGen();
    prepareMaskValues(input_pid_map);
    //prepareMaskTags(); 
  }

  bool OfflineEvaluator::TripleSacrifice() {
    bool res = true;
    std::vector<RepShare<Field>> v_rep_vec, w_rep_vec;
    std::vector<Field> r_vec, v_vec, w_vec;
    std::vector<Field> buffer;
    size_t mult_num = 0, idx_buffer = 0, idx = 0;

    for (const auto& level : circ_.gates_by_level) {
      for (const auto& gate : level) {
        switch (gate->type) {
          case GateType::kMul: {
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            auto *pre_mul = static_cast<PreprocMultGate<Field> *>(preproc_.gates[gate->out].get());
            auto mask_in1 = preproc_.gates[g->in1]->mask;
            Field r;
            randomizeZZp(rgen_.all(), r, sizeof(Field));
            r_vec.push_back(r);
            RepShare<Field> v = mask_in1 * r - pre_mul->ver;
            v_rep_vec.push_back(v);
            Field sum = Field(0);
            if (id_ == 0) {
              for (size_t i = 0; i < nP_; i++) {
                sum += v.getValues()[i];
              }
              buffer.push_back(sum);
            }
            else {
              buffer.push_back(v.getValues()[1]);
            }
            mult_num++;              
            break;
          }

          default: {
            break;
          }
        }
      }
    }

    if (id_ == 0){
      std::vector<std::future<void>> send_t;
      for (size_t i = 1; i <= nP_; i++) {
        send_t.push_back(tpool_->enqueue([&,i]() {
          network_->send(i, buffer.data(), sizeof(Field) * mult_num);
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
      network_->send(0, &mult_num, sizeof(size_t));
      network_->getSendChannel(0)->flush();
      network_->send(0, buffer.data(), sizeof(Field) * mult_num, true);
      network_->getSendChannel(0)->flush();
    }

    if (id_ == 0) {
      size_t count = 2;
      {
        std::unique_lock<std::mutex> lock(mtx_);
        cv_.wait(lock, [&]() { return offline_message_buffer_[count].size() >= 1; });
      }
      Offline_Message trip_messages = offline_message_buffer_[count].front();
      std::queue<Offline_Message> empty;
      std::swap(offline_message_buffer_[count], empty);

      std::vector<Field> v_val = trip_messages.data;
      for (size_t i = 0; i < mult_num; i++) {
        v_vec.push_back(v_val[i] + buffer[i]);
      }        
    }
    else {
      std::vector<Field> v_val(mult_num);
      network_->recv(0, v_val.data(), sizeof(Field) * mult_num);
      for (size_t i = 0; i < mult_num; i++) {
        v_vec.push_back(v_val[i] + buffer[i]);
      }
    }

    buffer.clear();
    buffer.shrink_to_fit();
    idx_buffer = 0;
    idx = 0;

    for (const auto& level : circ_.gates_by_level) {
      for (const auto& gate : level) {
        switch (gate->type) {
          case GateType::kMul: {
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            auto *pre_mul = static_cast<PreprocMultGate<Field> *>(preproc_.gates[gate->out].get());
            auto mask_in2 = preproc_.gates[g->in2]->mask;
            RepShare<Field> w = mask_in2 * v_vec[idx] - pre_mul->mask_prod * r_vec[idx] + pre_mul->ver_prod;
            idx++;
            w_rep_vec.push_back(w);
            Field sum = Field(0);
            if (id_ == 0) {
              for (size_t i = 0; i < nP_; i++) {
                sum += w.getValues()[i];
              }
              buffer.push_back(sum);
            }
            else {
              buffer.push_back(w.getValues()[1]);
            }              
            break;
          }

          default: {
            break;
          }
        }
      }
    }

    if (id_ == 0){
      std::vector<std::future<void>> send_t;
      for (size_t i = 1; i <= nP_; i++) {
        send_t.push_back(tpool_->enqueue([&,i]() {
          network_->send(i, buffer.data(), sizeof(Field) * mult_num);
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
      network_->send(0, &mult_num, sizeof(size_t));
      network_->getSendChannel(0)->flush();
      network_->send(0, buffer.data(), sizeof(Field) * mult_num, true);
      network_->getSendChannel(0)->flush();
    }

    if (id_ == 0) {
      size_t count = 3;
      {
        std::unique_lock<std::mutex> lock(mtx_);
        cv_.wait(lock, [&]() { return offline_message_buffer_[count].size() >= 1; });
      }
      Offline_Message trip_messages = offline_message_buffer_[count].front();
      std::queue<Offline_Message> empty;
      std::swap(offline_message_buffer_[count], empty);

      std::vector<Field> w_val = trip_messages.data;
      for (size_t i = 0; i < mult_num; i++) {
        w_vec.push_back(w_val[i] + buffer[i]);
      }        
    }
    else {
      std::vector<Field> w_val(mult_num);
      network_->recv(0, w_val.data(), sizeof(Field) * mult_num);
      for (size_t i = 0; i < mult_num; i++) {
        w_vec.push_back(w_val[i] + buffer[i]);
      }
    }

    buffer.clear();
    buffer.shrink_to_fit();

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

    Field omega = Field(0);

    if (id_ != 0) {
      Field rho;
      for (size_t i = 0; i < mult_num; i++) {
        randomizeZZp(rgen_.all(), rho, sizeof(Field));
        omega += rho * (v_vec[i] * v_rep_vec[i].getKeySh()[1] - v_rep_vec[i].getTags()[1]);
        randomizeZZp(rgen_.all(), rho, sizeof(Field));
        omega += rho * (w_vec[i] * w_rep_vec[i].getKeySh()[1] - w_rep_vec[i].getTags()[1]);
      }
      size_t count = 2;
      size_t total_comm = 1;
      network_->send(0, &total_comm, sizeof(size_t));
      network_->getSendChannel(0)->flush();
      network_->send(0, &omega, sizeof(Field), true);
      network_->getSendChannel(0)->flush();
    }
    else {
      size_t count = 4;
      {
        std::unique_lock<std::mutex> lock(mtx_);
        cv_.wait(lock, [&]() { return offline_message_buffer_[count].size() >= 1; });
      }
      Offline_Message trip_messages = offline_message_buffer_[count].front();
      std::queue<Offline_Message> empty;
      std::swap(offline_message_buffer_[count], empty);

      omega = trip_messages.data[0];
      
      Field rho_v, rho_w;
      for (size_t i = 0; i < mult_num; i++) {
        randomizeZZp(rgen_.all(), rho_v, sizeof(Field));
        randomizeZZp(rgen_.all(), rho_w, sizeof(Field));
        for (size_t j = 0; j < nP_; j++) {
          omega += rho_v * (v_vec[i] * v_rep_vec[i].getKeySh()[j] - v_rep_vec[i].getTags()[j]);
          omega += rho_w * (w_vec[i] * w_rep_vec[i].getKeySh()[j] - w_rep_vec[i].getTags()[j]);
        }
      }
    }   

    if (id_ == 0) {
      if (omega != 0) {
        res = false;
      }
      else {
        for (size_t i=0; i<mult_num; i++) {
          if (w_vec[i] != 0) {
            res = false;
            break;
          }
        }
      }

      std::vector<std::future<void>> send_t;
      for (size_t i = 1; i <= nP_; i++) {
        send_t.push_back(tpool_->enqueue([&,i]() {
          network_->send(i, &res, sizeof(bool));
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
      network_->recv(0, &res, sizeof(bool));
    }    

    return res;
  }

  PreprocCircuit<Field> OfflineEvaluator::run(const std::unordered_map<wire_t, int>& input_pid_map) {
    setWireMasks(input_pid_map);
    if (!TripleSacrifice()) {
      std::cout << "Malicious Activity Detected!!! Triple verification failed!!!" << std::endl;
      exit(0);
    }
    return std::move(preproc_);    
  }
}; // namespace dmAsyncAsterisk
