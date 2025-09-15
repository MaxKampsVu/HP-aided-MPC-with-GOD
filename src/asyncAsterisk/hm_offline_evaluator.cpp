#include "hm_offline_evaluator.h"

namespace hmAsyncAsterisk {
  OfflineEvaluator::OfflineEvaluator(int nP, int id, int security_param, std::shared_ptr<NetIOMP> network, LevelOrderedCircuit circ, 
    int threads, uint64_t seed) 
    : nP_(nP), th_((nP_-1)/2), id_(id), security_param_(security_param), rgen_(id, nP, seed), network_(std::move(network)), 
    circ_(std::move(circ)), preproc_(circ.num_gates) 
  {
    tpool_ = std::make_shared<ThreadPool>(threads);
  }

  PreprocCircuit<Field> OfflineEvaluator::dummy(int nP, int id, const LevelOrderedCircuit& circ, 
    const std::unordered_map<wire_t, int>& input_pid_map, PRG& prg) {
    int th = (nP-1)/2;
    PreprocCircuit<Field> preproc(circ.num_gates);
    DummyShare<Field> mask, mask_prod;

    std::vector<Field> key_sh, values, tags, tagPoints, valPoints;
    std::vector<Field> evalPoints1(th + 1), evalPoints2(th + 1);
    std::iota(evalPoints1.begin(), evalPoints1.end(), 1);
    std::iota(evalPoints2.begin(), evalPoints2.end(), 0);

    for (size_t i=0; i<=th; i++) {
      Field temp;
      randomizeZZp(prg, temp, sizeof(Field));
      key_sh.push_back(temp);
    }    
    fieldPoly poly_key = reconstructPolynomial(evalPoints1, key_sh);
    for (size_t i=th+1; i<nP; i++) {
      key_sh.push_back(eval(poly_key, Field(i+1)));
    }
    mask.setKeySh(key_sh);
    mask.setKey(coeff(poly_key, 0));

    for (size_t i=0; i<=th; i++) {
      Field temp;
      randomizeZZp(prg, temp, sizeof(Field));
      values.push_back(temp);
    }
    fieldPoly poly_val = reconstructPolynomial(evalPoints1, values);
    for (size_t i=th+1; i<nP; i++) {
      values.push_back(eval(poly_val, Field(i+1)));
    }
    mask.setValues(values);
    mask.setValue(coeff(poly_val, 0));

    tagPoints.push_back((poly_key, 0) * coeff(poly_val, 0));
    for (size_t i=0; i<th; i++) {
      Field temp;
      randomizeZZp(prg, temp, sizeof(Field));
      tagPoints.push_back(temp);
      tags.push_back(temp);
    }
    fieldPoly poly_tag = reconstructPolynomial(evalPoints2, tagPoints);
    for (size_t i=th; i<nP; i++) {
      tags.push_back(eval(poly_tag, Field(i+1)));
    }
    mask.setTags(tags);
    mask.setTag(tagPoints[0]);

    mask_prod.setKeySh(key_sh);
    mask_prod.setKey(mask.secretKey());

    Field secret = mask.secretVal() * mask.secretVal();
    valPoints.push_back(secret);
    for (size_t i=0; i<th; i++) {
      Field temp;
      randomizeZZp(prg, temp, sizeof(Field));
      valPoints.push_back(temp);
      values.push_back(temp);
    }
    fieldPoly poly_val_prod = reconstructPolynomial(evalPoints2, valPoints);
    for (size_t i=th; i<nP; i++) {
      values.push_back(eval(poly_val_prod, Field(i+1)));
    }
    mask_prod.setValues(values);
    mask_prod.setValue(secret);

    Field tag = secret * mask_prod.secretKey();
    tagPoints.clear();
    tagPoints.shrink_to_fit();
    tags.clear();
    tags.shrink_to_fit();
    tagPoints.push_back(tag);
    for (size_t i=0; i<th; i++) {
      Field temp;
      randomizeZZp(prg, temp, sizeof(Field));
      tagPoints.push_back(temp);
      tags.push_back(temp);
    }
    fieldPoly poly_tag_prod = reconstructPolynomial(evalPoints2, tagPoints);
    for (size_t i=th; i<nP; i++) {
      tags.push_back(eval(poly_tag_prod, Field(i+1)));
    }
    mask_prod.setTags(tags);
    mask_prod.setTag(tag);

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
            pregate->mask = mask.getASS(id);
            pregate->tpmask = mask.getTSS(id);
            preproc.gates[gate->out] = std::move(pregate);                
            break;
          }

          case GateType::kAdd: {
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            DummyShare<Field> temp = mask + mask;
            preproc.gates[gate->out] = std::make_unique<PreprocGate<Field>>(temp.getASS(id), temp.getTSS(id));    
            break;
          }
  
          case GateType::kConstAdd: {
            const auto* g = static_cast<ConstOpGate<Field>*>(gate.get());
            preproc.gates[gate->out] = std::make_unique<PreprocGate<Field>>(mask.getASS(id), mask.getTSS(id));
            break;
          }
  
          case GateType::kConstMul: {
            const auto* g = static_cast<ConstOpGate<Field>*>(gate.get());
            DummyShare<Field> temp = mask * g->cval;
            preproc.gates[gate->out] = std::make_unique<PreprocGate<Field>>(temp.getASS(id), temp.getTSS(id));
            break;
          }
  
          case GateType::kSub: {
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            DummyShare<Field> temp = mask - mask;
            preproc.gates[gate->out] = std::make_unique<PreprocGate<Field>>(temp.getASS(id), temp.getTSS(id));
            break;
          }

          case GateType::kMul: {
            preproc.gates[gate->out] = std::make_unique<PreprocMultGate<Field>>();
            preproc.gates[gate->out] = std::move(std::make_unique<PreprocMultGate<Field>> (mask.getASS(id), mask.getTSS(id), mask_prod.getASS(id), mask_prod.getTSS(id)));
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

  void OfflineEvaluator::keyGen(std::vector<Field>& keySh, Field& key) {
    if(id_ == 0)  {
      key = 0;
      keySh[0] = 0;
      std::vector<Field> keyShs;
      randomizeZZp(rgen_.self(), key, sizeof(Field));
      key_sh_ = key;
      keyShs.push_back(key);
      for(int i = 1; i <= th_; i++) {
        randomizeZZp(rgen_.pi(i), keySh[i], sizeof(Field));
        keyShs.push_back(keySh[i]);
      }            

      std::vector<Field> evalPoints(th_ + 1);
      std::iota(evalPoints.begin(), evalPoints.end(), 0);
      fieldPoly poly_key = reconstructPolynomial(evalPoints, keyShs);
      std::vector<std::future<void>> send_keysh_t;
      for(size_t i = th_+1; i <= nP_; i++) {
        keySh[i] = eval(poly_key, Field(i));
        send_keysh_t.push_back(tpool_->enqueue([&,i](){
          network_->send(i, &keySh[i], sizeof(Field));
          network_->getSendChannel(i)->flush();
        }));
      }
      for(auto& t : send_keysh_t) {
        if (t.valid()) {
          t.wait();
        }
      }
    }
    else if (id_ <= th_) {
      randomizeZZp(rgen_.p0(), key, sizeof(Field));
      key_sh_ = key;
    }
    else {
      network_->recv(0, &key, sizeof(Field));
      key_sh_ = key;
    }
  }

  void OfflineEvaluator::randomShare(int nP, int pid, RandGenPool& rgen, NetIOMP& network, AuthShamirShare<Field>& share, 
    TPShamirShare<Field>& tpShare,  Field key, std::vector<Field> keySh, std::vector<std::vector<Field>>& rand_sh, 
    std::vector<size_t>& idx_rand_sh) {

    Field secret = Field(0);
    if(pid == 0) {
        randomizeZZp(rgen.self(), secret, sizeof(Field));            
    }
    randomShareSecret(nP, pid, rgen, network, share, tpShare, secret, key, keySh, rand_sh, idx_rand_sh);
  }

  void OfflineEvaluator::randomShareSecret(int nP, int pid, RandGenPool& rgen, NetIOMP& network, AuthShamirShare<Field>& share, 
    TPShamirShare<Field>& tpShare, Field secret, Field key, std::vector<Field> keySh, std::vector<std::vector<Field>>& rand_sh_sec, 
    std::vector<size_t>& idx_rand_sh_sec) {

    int th = (nP-1)/2;
    
    Field val = Field(0);
    Field tag = Field(0);

    if(pid == 0) {
      share.pushValue(Field(0));
      share.pushTag(Field(0));
      share.setKey(keySh[0]);
      tpShare.pushValues(Field(0));
      tpShare.pushTags(Field(0));
      tpShare.setKeySh(keySh[0]);
      tpShare.setKey(key);

      for(int i = 1; i <= th; i++) {
          randomizeZZp(rgen.pi(i), val, sizeof(Field));
          tpShare.pushValues(val);
          tpShare.setKeySh(keySh[i]);
          randomizeZZp(rgen.pi(i), tag, sizeof(Field));
          tpShare.pushTags(tag);
      }
      
      tag = key * secret;
      std::vector<Field> evalPoints(th + 1);
      std::iota(evalPoints.begin(), evalPoints.end(), 0);
      std::vector<Field> vals, tags;
      vals.push_back(secret);
      tags.push_back(tag);
      for(size_t i = 1; i <= th; i++) {
          vals.push_back(tpShare.commonValueWithParty(i));
          tags.push_back(tpShare.commonTagWithParty(i));
      }
      fieldPoly poly_val = reconstructPolynomial(evalPoints, vals);
      fieldPoly poly_tag = reconstructPolynomial(evalPoints, tags);

      for(size_t i = th+1; i <= nP; i++) {
          Field val = eval(poly_val, Field(i));
          tpShare.pushValues(val);
          rand_sh_sec[i-th-1].push_back(val);
          tpShare.setKeySh(keySh[i]);
          Field tag = eval(poly_tag, Field(i));
          tpShare.pushTags(tag);
          rand_sh_sec[i-th-1].push_back(tag);
      }
    }
    else if(pid <= th) {
      share.setKey(key);
      randomizeZZp(rgen.p0(), val, sizeof(Field));
      share.pushValue(val);
      randomizeZZp(rgen.p0(), tag, sizeof(Field));
      share.pushTag(tag);
    }
    else {
      share.setKey(key);
      int index = pid-th-1;
      Field val = rand_sh_sec[index][idx_rand_sh_sec[index]];
      idx_rand_sh_sec[index]++;
      share.pushValue(val);
      Field tag = rand_sh_sec[index][idx_rand_sh_sec[index]];
      idx_rand_sh_sec[index]++;
      share.pushTag(tag);
    }
  }  

  void OfflineEvaluator::randomShareWithParty(int nP, int pid, int dealer, RandGenPool& rgen, NetIOMP& network, AuthShamirShare<Field>& share, 
    TPShamirShare<Field>& tpShare, Field& secret, Field key, std::vector<Field> keySh, std::vector<std::vector<Field>>& rand_sh_party, 
    std::vector<size_t>& idx_rand_sh_party) {
      
    if(pid == 0) {
      if(dealer != 0) {
        randomizeZZp(rgen.pi(dealer), secret, sizeof(Field));
      }
      else {
        randomizeZZp(rgen.self(), secret, sizeof(Field));
      }            
    }
    else if(pid == dealer) {
      randomizeZZp(rgen.p0(), secret, sizeof(Field));
    }
    randomShareSecret(nP, pid, rgen, network, share, tpShare, secret, key, keySh, rand_sh_party, idx_rand_sh_party);
  }

  void OfflineEvaluator::setWireMasksParty(const std::unordered_map<wire_t,int>& input_pid_map, Field key, std::vector<Field> keySh, 
    std::vector<std::vector<Field>>& rand_sh, std::vector<std::vector<Field>>& rand_sh_sec, std::vector<std::vector<Field>>& rand_sh_party) {    
        
    std::vector<size_t> idx_rand_sh(nP_-th_, 0), idx_rand_sh_sec(nP_-th_, 0), idx_rand_sh_party(nP_-th_, 0);

    for (const auto& level : circ_.gates_by_level) {
      for (const auto& gate : level) {
        switch (gate->type) {
          case GateType::kInp: {
            auto pregate = std::make_unique<PreprocInput<Field>>();      
            auto pid = input_pid_map.at(gate->out);
            pregate->pid = pid;
            randomShareWithParty(nP_, id_, pid, rgen_, *network_, pregate->mask, pregate->tpmask, pregate->mask_value, key, keySh, rand_sh_party, idx_rand_sh_party);      
            preproc_.gates[gate->out] = std::move(pregate);                
            break;
          }
  
          case GateType::kAdd: {
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            const auto& mask_in1 = preproc_.gates[g->in1]->mask;
            const auto& tpmask_in1 = preproc_.gates[g->in1]->tpmask;
            const auto& mask_in2 = preproc_.gates[g->in2]->mask;
            const auto& tpmask_in2 = preproc_.gates[g->in2]->tpmask;
            preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask_in1 + mask_in2), (tpmask_in1 + tpmask_in2));
            break;
          }
  
          case GateType::kConstAdd: {
            const auto* g = static_cast<ConstOpGate<Field>*>(gate.get());
            const auto& mask = preproc_.gates[g->in]->mask;
            const auto& tpmask = preproc_.gates[g->in]->tpmask;
            preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask), (tpmask));
            break;
          }
  
          case GateType::kConstMul: {
            const auto* g = static_cast<ConstOpGate<Field>*>(gate.get());
            const auto& mask = preproc_.gates[g->in]->mask * g->cval;
            const auto& tpmask = preproc_.gates[g->in]->tpmask * g->cval;
            preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask), (tpmask));
            break;
          }
  
          case GateType::kSub: {
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            const auto& mask_in1 = preproc_.gates[g->in1]->mask;
            const auto& tpmask_in1 = preproc_.gates[g->in1]->tpmask;
            const auto& mask_in2 = preproc_.gates[g->in2]->mask;
            const auto& tpmask_in2 = preproc_.gates[g->in2]->tpmask;
            preproc_.gates[gate->out] = std::make_unique<PreprocGate<Field>>((mask_in1 - mask_in2),(tpmask_in1 - tpmask_in2));
            break;
          }
  
          case GateType::kMul: {
            preproc_.gates[gate->out] = std::make_unique<PreprocMultGate<Field>>();
            const auto* g = static_cast<FIn2Gate*>(gate.get());
            const auto& mask_in1 = preproc_.gates[g->in1]->mask;
            const auto& tpmask_in1 = preproc_.gates[g->in1]->tpmask;
            const auto& mask_in2 = preproc_.gates[g->in2]->mask;
            const auto& tpmask_in2 = preproc_.gates[g->in2]->tpmask;
            Field tp_prod;
            if(id_ == 0) {tp_prod = tpmask_in1.secret() * tpmask_in2.secret();}
            TPShamirShare<Field> tprand_mask;
            AuthShamirShare<Field> rand_mask;
            randomShare(nP_, id_, rgen_, *network_, rand_mask, tprand_mask, key, keySh, rand_sh, idx_rand_sh);
            TPShamirShare<Field> tpmask_product;
            AuthShamirShare<Field> mask_product; 
            randomShareSecret(nP_, id_, rgen_, *network_, mask_product, tpmask_product, tp_prod, key, keySh, rand_sh_sec, idx_rand_sh_sec);
            preproc_.gates[gate->out] = std::move(std::make_unique<PreprocMultGate<Field>> (rand_mask, tprand_mask, mask_product, tpmask_product));
            break;
          }
  
          default: {
            break;
          }
        }
      }
    }
  }

  void OfflineEvaluator::setWireMasks(const std::unordered_map<wire_t, int>& input_pid_map) {

    // key setup
    std::vector<Field> keySh(nP_ + 1);
    Field key = Field(0);
    keyGen(keySh, key);   
      
    std::vector<std::vector<Field>> rand_sh(nP_-th_), rand_sh_sec(nP_-th_), rand_sh_party(nP_-th_);
      
    if(id_ <= th_) {
      setWireMasksParty(input_pid_map, key, keySh, rand_sh, rand_sh_sec, rand_sh_party);
    
      if(id_ == 0) {
        std::vector<std::future<void>> send_t;
        for (size_t i = 0; i < nP_-th_; i++) {
          size_t rand_sh_num, rand_sh_sec_num, rand_sh_party_num, arith_comm;
          send_t.push_back(tpool_->enqueue([&,i]() {
            rand_sh_num = rand_sh[i].size();
            rand_sh_sec_num = rand_sh_sec[i].size();
            rand_sh_party_num = rand_sh_party[i].size();
            arith_comm = rand_sh_num + rand_sh_sec_num + rand_sh_party_num;

            std::vector<size_t> lengths(4);
            lengths[0] = rand_sh_num;
            lengths[1] = rand_sh_sec_num;
            lengths[2] = rand_sh_party_num;
            network_->send(i+th_+1, lengths.data(), sizeof(size_t) * 3);
            network_->getSendChannel(i+th_+1)->flush();
           
            std::vector<Field> offline_arith_comm(arith_comm);

            for(size_t j = 0; j < rand_sh_num; j++) {
              offline_arith_comm[j] = rand_sh[i][j];
            }
            for(size_t j = 0; j < rand_sh_sec_num; j++) {
              offline_arith_comm[rand_sh_num + j] = rand_sh_sec[i][j];
            }
            for(size_t j = 0; j < rand_sh_party_num; j++) {
              offline_arith_comm[rand_sh_sec_num + rand_sh_num + j] = rand_sh_party[i][j];
            }

            network_->send(i+th_+1, offline_arith_comm.data(), sizeof(Field) * arith_comm);
            network_->getSendChannel(i+th_+1)->flush();
          }));
        }

        for(auto& t : send_t) {
          if (t.valid()) {
            t.wait();
          }
        }
      }
    }
    else {      
      std::vector<size_t> lengths(4);      
      network_->recv(0, lengths.data(), sizeof(size_t) * 3);      
      size_t rand_sh_num = lengths[0];
      size_t rand_sh_sec_num = lengths[1];
      size_t rand_sh_party_num = lengths[2];
      size_t arith_comm = rand_sh_num + rand_sh_sec_num + rand_sh_party_num;

      std::vector<Field> offline_arith_comm(arith_comm);
      network_->recv(0, offline_arith_comm.data(), sizeof(Field) * arith_comm);  

      rand_sh[id_-th_-1].resize(rand_sh_num);
      for(size_t i = 0; i < rand_sh_num; i++) {
        rand_sh[id_-th_-1][i] = offline_arith_comm[i];
      }
      
      rand_sh_sec[id_-th_-1].resize(rand_sh_sec_num);
      for(size_t i = 0; i < rand_sh_sec_num; i++) {
        rand_sh_sec[id_-th_-1][i] = offline_arith_comm[rand_sh_num + i];
      }
      
      rand_sh_party[id_-th_-1].resize(rand_sh_party_num);
      for(size_t i = 0; i < rand_sh_party_num; i++) {
        rand_sh_party[id_-th_-1][i] = offline_arith_comm[rand_sh_num + rand_sh_sec_num + i];
      }      
      
      setWireMasksParty(input_pid_map, key, keySh, rand_sh, rand_sh_sec, rand_sh_party);
    }    
  }

  PreprocCircuit<Field> OfflineEvaluator::run(
    const std::unordered_map<wire_t, int>& input_pid_map) {
    setWireMasks(input_pid_map);
    return std::move(preproc_);    
  }
}; // namespace hmAsyncAsterisk
