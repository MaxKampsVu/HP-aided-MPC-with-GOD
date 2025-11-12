#include "dmgod_offline_evaluator.h"

namespace dmAsyncAsteriskGOD {
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
          // TODO: do OLE on tags later  
          // for (size_t count=0; count < 2; count++) {
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
            // I guess notify different waiting thread to aquire lock  
            cv_.notify_one();
          }
        });
      }
    }
  }

  OfflineEvaluator::~OfflineEvaluator() {
    tpool_.reset();
  }

  void OfflineEvaluator::keyGen()  {
    std::cout << id_ << ": I did the key gen" << std::endl;
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

  void OfflineEvaluator::randomShare(int pid, RandGenPool& rgen, TwoShare<Field>& share, Field& mask_share_zero, bool isOutputWire) {
    // TP 
    if(pid == 0) {      
      Field valSh;
      randomizeZZp(rgen.p0(), valSh, sizeof(Field));
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
    inputToOPE.push_back(share1_val);
    inputToOPE.push_back(share2_val);
    prodShare.setValue(Field(0));
  }  

  void OfflineEvaluator::randomShareWithParty(int pid, int dealer, RandGenPool& rgen, TwoShare<Field>& share, Field& secret) {        
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

  void OfflineEvaluator::multSS(const Field& share1, const Field& share2, Field& output, 
    const std::vector<Field>& outputOfOPE, size_t& idx_outputOfOPE) {
      output = share1 * share2 + outputOfOPE[idx_outputOfOPE] + outputOfOPE[idx_outputOfOPE+1];
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
            randomShareWithParty(id_, pid, rgen_, pregate->mask, pregate->mask_value);      
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
            TwoShare<Field> rand_mask; 
            Field mask_share_zero = Field(0);
            bool isOutputWire = false;
            if (std::find(circ_.outputs.begin(), circ_.outputs.end(),g->out)!=circ_.outputs.end())
              isOutputWire = true;
            randomShare(id_, rgen_, rand_mask, mask_share_zero, isOutputWire);
            TwoShare<Field> mask_product;
            std::cout << id_ << ": starts randomShareSecret" << std::endl;
            randomShareSecret(id_, rgen_, mask_in1, mask_in2, mask_product, inputToOPE[0]);
            auto dummy = TwoShare<Field>();
            preproc_.gates[gate->out] = std::move(std::make_unique<PreprocMultGate<Field>> (rand_mask, mask_product, mask_share_zero, dummy, dummy));
            break;
          }

          default: {
            break;
          }
        }
      }
    }
  
    // Output after MSSR OLE is one Field element 
    std::vector<Field> outputOfOPE;
    size_t idx_outputOfOPE = 0;
    
    std::cout << id_ << ": I started runOPE" << std::endl;
    runOPE(inputToOPE[0], outputOfOPE, 0);
    std::cout << id_ << ": I finished runOPE" << std::endl;

    // This is only for Castor where the last Party receives a different share 
    // // if (id_ != nP_) {
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

//   void OfflineEvaluator::prepareMaskTags() {
//     std::vector<Field> buffer; 
//     size_t idx_buffer=0; 
//     size_t buffer_num=0;

//     for (const auto& level : circ_.gates_by_level) {
//       for (const auto& gate : level) {
//         switch (gate->type) {
//           case GateType::kInp: {
//             if (id_ == nP_) {
//               buffer_num += 1;
//             }
//             auto *pre_input = static_cast<PreprocInput<Field> *>(preproc_.gates[gate->out].get());
//             auto share1 = pre_input->mask.getValues();
//             auto share2 = pre_input->mask.getKeySh();
//             if (id_!=0) {      
//               inputToOPE[1].push_back(share1[1]);
//               inputToOPE[1].push_back(share2[1]);
//             }
//             else {
//               Field a = Field(0), b = Field(0);
//               for (size_t i=0; i < nP_; i++) {
//                 a += share1[i];
//                 b += share2[i];
//               }
//               inputToOPE[1].push_back(b);
//               inputToOPE[1].push_back(a);
//             }
//             break;
//           }

//           case GateType::kMul: {
//             if (id_ == nP_) {
//               buffer_num += 4;
//             }
//             auto *pre_mul = static_cast<PreprocMultGate<Field> *>(preproc_.gates[gate->out].get());
//             std::vector<Field> share1, share2;
//             share1 = pre_mul->mask.getValues();
//             share2 = pre_mul->mask.getKeySh();
//             if (id_!=0) {      
//               inputToOPE[1].push_back(share1[1]);
//               inputToOPE[1].push_back(share2[1]);
//             }
//             else {
//               Field a = Field(0), b = Field(0);
//               for (size_t i=0; i < nP_; i++) {
//                 a += share1[i];
//                 b += share2[i];
//               }
//               inputToOPE[1].push_back(b);
//               inputToOPE[1].push_back(a);
//             }
//             share1 = pre_mul->mask_prod.getValues();
//             share2 = pre_mul->mask_prod.getKeySh();
//             if (id_!=0) {      
//               inputToOPE[1].push_back(share1[1]);
//               inputToOPE[1].push_back(share2[1]);
//             }
//             else {
//               Field a = Field(0), b = Field(0);
//               for (size_t i=0; i < nP_; i++) {
//                 a += share1[i];
//                 b += share2[i];
//               }
//               inputToOPE[1].push_back(b);
//               inputToOPE[1].push_back(a);
//             }
//             share1 = pre_mul->ver.getValues();
//             share2 = pre_mul->ver.getKeySh();
//             if (id_!=0) {      
//               inputToOPE[1].push_back(share1[1]);
//               inputToOPE[1].push_back(share2[1]);
//             }
//             else {
//               Field a = Field(0), b = Field(0);
//               for (size_t i=0; i < nP_; i++) {
//                 a += share1[i];
//                 b += share2[i];
//               }
//               inputToOPE[1].push_back(b);
//               inputToOPE[1].push_back(a);
//             }
//             share1 = pre_mul->ver_prod.getValues();
//             share2 = pre_mul->ver_prod.getKeySh();
//             if (id_!=0) {      
//               inputToOPE[1].push_back(share1[1]);
//               inputToOPE[1].push_back(share2[1]);
//             }
//             else {
//               Field a = Field(0), b = Field(0);
//               for (size_t i=0; i < nP_; i++) {
//                 a += share1[i];
//                 b += share2[i];
//               }
//               inputToOPE[1].push_back(b);
//               inputToOPE[1].push_back(a);
//             }
//             break;
//           }

//           default: {
//             break;
//           }
//         }
//       }
//     }

//     std::vector<Field> outputOfOPE;
//     size_t idx_outputOfOPE = 0;

//     runOPE(inputToOPE[1], outputOfOPE, 1);

//     if (id_ != nP_) {
//       for (const auto& level : circ_.gates_by_level) {
//         for (const auto& gate : level) {
//           switch (gate->type) {
//             case GateType::kInp: {
//               auto *pre_input = static_cast<PreprocInput<Field> *>(preproc_.gates[gate->out].get());
//               std::vector<Field> tags;
//               auto share1 = pre_input->mask.getValues();;
//               auto share2 = pre_input->mask.getKeySh();
//               tags.resize(share1.size());
//               multiply(share1, share2, tags, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
//               pre_input->mask.setTags(tags);                
//               break;
//             }

//             case GateType::kAdd: {
//               const auto* g = static_cast<FIn2Gate*>(gate.get());
//               const auto& mask_in1 = preproc_.gates[g->in1]->mask;
//               const auto& mask_in2 = preproc_.gates[g->in2]->mask;
//               preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask_in1 + mask_in2));    
//               break;
//             }
  
//             case GateType::kConstAdd: {
//               const auto* g = static_cast<ConstOpGate<Field>*>(gate.get());
//               const auto& mask = preproc_.gates[g->in]->mask;
//               preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask));
//               break;
//             }
  
//             case GateType::kConstMul: {
//               const auto* g = static_cast<ConstOpGate<Field>*>(gate.get());
//               const auto& mask = preproc_.gates[g->in]->mask * g->cval;
//               preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask));
//               break;
//             }
  
//             case GateType::kSub: {
//               const auto* g = static_cast<FIn2Gate*>(gate.get());
//               const auto& mask_in1 = preproc_.gates[g->in1]->mask;
//               const auto& mask_in2 = preproc_.gates[g->in2]->mask;
//               preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask_in1 - mask_in2));
//               break;
//             }

//             case GateType::kMul: {
//               auto *pre_mul = static_cast<PreprocMultGate<Field> *>(preproc_.gates[gate->out].get());
//               std::vector<Field> share1, share2, tags1, tags2, tags3, tags4;
//               share1 = pre_mul->mask.getValues();
//               share2 = pre_mul->mask.getKeySh();
//               tags1.resize(share1.size());
//               multiply(share1, share2, tags1, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
//               pre_mul->mask.setTags(tags1);
//               share1 = pre_mul->mask_prod.getValues();
//               share2 = pre_mul->mask_prod.getKeySh();
//               tags2.resize(share1.size());
//               multiply(share1, share2, tags2, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
//               pre_mul->mask_prod.setTags(tags2);
//               share1 = pre_mul->ver.getValues();
//               share2 = pre_mul->ver.getKeySh();
//               tags3.resize(share1.size());
//               multiply(share1, share2, tags3, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
//               pre_mul->ver.setTags(tags3);
//               share1 = pre_mul->ver_prod.getValues();
//               share2 = pre_mul->ver_prod.getKeySh();
//               tags4.resize(share1.size());
//               multiply(share1, share2, tags4, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
//               pre_mul->ver_prod.setTags(tags4);            
//               break;
//             }
    
//             default: {
//               break;
//             }
//           }
//         }
//       }
//       outputOfOPE.clear();
//       outputOfOPE.shrink_to_fit();

//       if (id_ == 0) {
//         buffer_num = buffer.size();
//         network_->send(nP_, buffer.data(), sizeof(Field) * buffer_num);
//         network_->getSendChannel(nP_)->flush();
//         buffer.clear();
//         buffer.shrink_to_fit();
//       }      
//     }
//     else {
//       buffer.resize(buffer_num);
//       network_->recv(0, buffer.data(), sizeof(Field) * buffer_num);

//       for (const auto& level : circ_.gates_by_level) {
//         for (const auto& gate : level) {
//           switch (gate->type) {
//             case GateType::kInp: {
//               auto *pre_input = static_cast<PreprocInput<Field> *>(preproc_.gates[gate->out].get());
//               std::vector<Field> tags;
//               auto share1 = pre_input->mask.getValues();
//               auto share2 = pre_input->mask.getKeySh();
//               tags.resize(share1.size());
//               multiply(share1, share2, tags, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
//               pre_input->mask.setTags(tags);                
//               break;
//             }

//             case GateType::kAdd: {
//               const auto* g = static_cast<FIn2Gate*>(gate.get());
//               const auto& mask_in1 = preproc_.gates[g->in1]->mask;
//               const auto& mask_in2 = preproc_.gates[g->in2]->mask;
//               preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask_in1 + mask_in2));    
//               break;
//             }
  
//             case GateType::kConstAdd: {
//               const auto* g = static_cast<ConstOpGate<Field>*>(gate.get());
//               const auto& mask = preproc_.gates[g->in]->mask;
//               preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask));
//               break;
//             }
  
//             case GateType::kConstMul: {
//               const auto* g = static_cast<ConstOpGate<Field>*>(gate.get());
//               const auto& mask = preproc_.gates[g->in]->mask * g->cval;
//               preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask));
//               break;
//             }
  
//             case GateType::kSub: {
//               const auto* g = static_cast<FIn2Gate*>(gate.get());
//               const auto& mask_in1 = preproc_.gates[g->in1]->mask;
//               const auto& mask_in2 = preproc_.gates[g->in2]->mask;
//               preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask_in1 - mask_in2));
//               break;
//             }

//             case GateType::kMul: {
//               auto *pre_mul = static_cast<PreprocMultGate<Field> *>(preproc_.gates[gate->out].get());
//               std::vector<Field> share1, share2, tags1, tags2, tags3, tags4;
//               share1 = pre_mul->mask.getValues();
//               share2 = pre_mul->mask.getKeySh();
//               tags1.resize(share1.size());
//               multiply(share1, share2, tags1, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
//               pre_mul->mask.setTags(tags1);
//               share1 = pre_mul->mask_prod.getValues();
//               share2 = pre_mul->mask_prod.getKeySh();
//               tags2.resize(share1.size())buffer;
//               multiply(share1, share2, tags2, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
//               pre_mul->mask_prod.setTags(tags2);
//               share1 = pre_mul->ver.getValues();
//               share2 = pre_mul->ver.getKeySh();
//               tags3.resize(share1.size());
//               multiply(share1, share2, tags3, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
//               pre_mul->ver.setTags(tagsbuffer3);
//               share1 = pre_mul->ver_prod.getValues();
//               share2 = pre_mul->ver_prod.getKeySh();
//               tags4.resize(share1.size());
//               multiply(share1, share2, tags4, outputOfOPE, buffer, idx_outputOfOPE, idx_buffer);
//               pre_mul->ver_prod.setTags(tags4);
//               break;
//             }

//             default: {
//               break;
//             }
//           }
//         }
//       }
//       outputOfOPE.clear();
//       outputOfOPE.shrink_to_fit();
//       buffer.clear();
//       buffer.shrink_to_fit();
//     }
//   }
    
  void OfflineEvaluator::setWireMasks(const std::unordered_map<wire_t,int>& input_pid_map) {      
    keyGen();
    prepareMaskValues(input_pid_map);
    //prepareMaskTags(); 

  }

  PreprocCircuit<Field> OfflineEvaluator::run(const std::unordered_map<wire_t, int>& input_pid_map) {
      setWireMasks(input_pid_map);
      return std::move(preproc_);    
  }

}; // namespace dmAsyncAsteriskGOD
