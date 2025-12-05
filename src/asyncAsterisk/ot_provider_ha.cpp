#include "ot_provider_ha.h"

namespace asyncAsterisk {
  constexpr size_t ot_bsize = emp::ot_bsize;

  OTProviderHA::OTProviderHA(int my_id, int other_id, NetIO* io) : ios_{io} {
    if (my_id == 0) {
      std::string filename = "./ot_data/s" + std::to_string(other_id) + "_r" + std::to_string(my_id) + "_receiver.bin";
      ot_ = std::make_unique<FerretCOT<NetIO>>(BOB, 1, ios_.data(), true, true, ferret_b13, filename);
    }
    else {
      std::string filename = "./ot_data/s" + std::to_string(my_id) + "_r" + std::to_string(other_id) + "_sender.bin";
      ot_ = std::make_unique<FerretCOT<NetIO>>(ALICE, 1, ios_.data(), true, true, ferret_b13, filename);
    }
  }

  void OTProviderHA::sendFieldElements(const Field* data, size_t length) {
    std::vector<uint8_t> serialized(length);
    size_t num = (length + FIELDSIZE - 1) / FIELDSIZE;
    for (size_t i = 0; i < num; ++i) {
      BytesFromZZ(serialized.data() + i * FIELDSIZE, conv<ZZ>(data[i]), FIELDSIZE);
    }
    ios_[0]->send_data(serialized.data(), serialized.size());
    ios_[0]->flush();
  }
        
  void OTProviderHA::recvFieldElements(Field* data, size_t length) {
    std::vector<uint8_t> serialized(length);
    ios_[0]->recv_data(serialized.data(), serialized.size());
    size_t num = (length + FIELDSIZE - 1) / FIELDSIZE;
    for (size_t i = 0; i < num; ++i) {
      data[i] = conv<Field>(ZZFromBytes(serialized.data() + i * FIELDSIZE, FIELDSIZE));
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

  std::vector<Field> blockToFields(emp::block b) {
    alignas(16) uint64_t tmp[2];
    _mm_storeu_si128(reinterpret_cast<__m128i*>(tmp), b);
    return {Field(tmp[0]), Field(tmp[1])};
  }


  std::vector<Field> blockArrayToFields(emp::block* arr, size_t size) {
    std::vector<Field> flat_result;

    for (size_t i = 0; i < size; ++i) {
      std::vector<Field> tmp = blockToFields(arr[i]);
      flat_result.insert(flat_result.end(), tmp.begin(), tmp.end());
    }

    return flat_result;
  }
  
  void OTProviderHA::send(const Field* data0, const Field* data1, size_t length, PRG& prg, fieldDig& ot_dig) {
    auto* data = new emp::block[length];    
    ot_->send_cot(data, length);
    emp::block s = {0, 0};
    prg.random_block(&s, 1); 
    ios_[0]->send_block(&s, 1);
    ot_->mitccrh.setS(s);
    ios_[0]->flush();

    block pad[2 * ot_bsize];
    std::vector<Field> upad(2 * length);
    auto* tpad = reinterpret_cast<uint64_t*>(pad);
    for (size_t i = 0; i < length; i += ot_bsize) {
      for (size_t j = i; j < std::min(i + ot_bsize, length); j++) {
        pad[2 * (j - i)] = data[j];
        pad[2 * (j - i) + 1] = data[j] ^ ot_->Delta; // TODO: make sure delta is the same 
      }
      hashFields(blockArrayToFields(pad, 2 * ot_bsize)); // TODO: append once deterministic 
      ot_->mitccrh.hash<ot_bsize, 2>(pad);
      for (size_t j = i; j < std::min(i + ot_bsize, length); j++) {
        upad[2 * j] = Field(tpad[4 * (j - i)]) + data0[j];
        upad[2 * j + 1] = Field(tpad[4 * (j - i) + 2]) + data1[j];
        
      }
    }
    hashFields(upad); // TODO: append once deterministic 
    ot_dig = hashFields(blockToFields(s)); // TODO: append upad 
    sendFieldElements(upad.data(), 2 * sizeof(Field) * length);
    delete[] data;
  }

  void OTProviderHA::recv(Field* rdata, const bool* r, size_t length, fieldDig& ot_dig) {
    auto* data = new emp::block[length];
    ot_->recv_cot(data, r, length);    
    emp::block s;
    ios_[0]->recv_block(&s, 1);
    ot_->mitccrh.setS(s);

    block pad[ot_bsize];
    std::vector<Field> res(2 * length);
    recvFieldElements(res.data(), 2 * sizeof(Field) * length);
    auto* tpad = reinterpret_cast<uint64_t*>(pad);
    for (size_t i = 0; i < length; i += ot_bsize) {
      memcpy(pad, data + i, std::min(ot_bsize, length - i) * sizeof(block));
      ot_->mitccrh.hash<ot_bsize, 1>(pad);
      for (size_t j = 0; j < ot_bsize and j < length - i; j++) {
        rdata[i + j] = res[2 * (i + j) + r[i + j]] - Field(tpad[2 * j]);
      }
    }
    ot_dig = hashFields(blockToFields(s)); // TODO: append upad 
    delete[] data;
  }

  std::vector<Field> OTProviderHA::multiplySendOnline(const std::vector<Field>& inputs, PRG& prg, fieldDig& ot_dig) {
    size_t num_bits = sizeof(Field) * 8;
    size_t num_blocks = num_bits * inputs.size();
  
    std::vector<Field> vrand(num_blocks);
    for (size_t i = 0; i < num_blocks; i++) {
      randomizeZZp(prg, vrand[i], sizeof(Field));
    }
  
    std::vector<Field> inp_0(num_blocks);
    std::vector<Field> inp_1(num_blocks);
    std::vector<Field> shares(inputs.size(), Field(0));
    size_t idx = 0;
    for (size_t i = 0; i < inputs.size(); i++) {
      const auto& input = inputs[i];
      auto& share = shares[i];
      for (size_t j = 0; j < num_bits; j++) {
        auto val = vrand[idx];
        share -= val;
        inp_0[idx] = val;
        inp_1[idx] = power(Field(2), conv<ZZ>(j))*input + val;
        idx++;
      }
    }

    send(inp_0.data(), inp_1.data(), num_blocks, prg, ot_dig);
    return shares;
  }

  std::vector<Field> OTProviderHA::multiplySendOffline(const std::vector<Field> inputs, PRG& prg, fieldDig& ot_dig) {
    // Compute my share without sending anything to HP 
    size_t num_bits = sizeof(Field) * 8;
    size_t num_blocks = num_bits * inputs.size();
  
    std::vector<Field> vrand(num_blocks);
    for (size_t i = 0; i < num_blocks; i++) {
      randomizeZZp(prg, vrand[i], sizeof(Field));
    }
  
    std::vector<Field> inp_0(num_blocks);
    std::vector<Field> inp_1(num_blocks);
    std::vector<Field> shares(inputs.size(), Field(0));
    size_t idx = 0;
    for (size_t i = 0; i < inputs.size(); i++) {
      const auto& input = inputs[i];
      auto& share = shares[i];
      for (size_t j = 0; j < num_bits; j++) {
        auto val = vrand[idx];
        share -= val;
        inp_0[idx] = val;
        inp_1[idx] = power(Field(2), conv<ZZ>(j))*input + val;
        idx++;
      }
    }

    // Compute the digest on expected message to HP 
    auto length = num_blocks;
    auto data0 = inp_0.data();
    auto data1 = inp_1.data();
    auto* data = new emp::block[length];    
    emp::block s;
    prg.random_block(&s, 1); 

    block pad[2 * ot_bsize];
    std::vector<Field> upad(2 * length);
    auto* tpad = reinterpret_cast<uint64_t*>(pad);

    for (size_t i = 0; i < length; i += ot_bsize) {
      for (size_t j = i; j < std::min(i + ot_bsize, length); j++) {
        pad[2 * (j - i)] = data[j];
        pad[2 * (j - i) + 1] = data[j] ^ ot_->Delta; // TODO: make sure Delta is the same 
      }
      hashFields(blockArrayToFields(pad, 2 * ot_bsize)); // TODO: append once deterministic 
      for (size_t j = i; j < std::min(i + ot_bsize, length); j++) {
        upad[2 * j] = Field(tpad[4 * (j - i)]) + data0[j];
        upad[2 * j + 1] = Field(tpad[4 * (j - i) + 2]) + data1[j];
      }
    }

    hashFields(upad); // TODO: append once deterministic 
    ot_dig = hashFields(blockToFields(s)); // TODO: append upad 
    delete[] data;
    return shares;
  }


  std::vector<Field> OTProviderHA::multiplySend(const std::vector<Field>& inputs, PRG& prg, fieldDig& ot_dig, bool run_async, size_t pid) {
    if(run_async) {
      return multiplySendOnline(inputs, prg, ot_dig);
    } 
    else {
      return pid == SYNC_SENDER_PID ? multiplySendOnline(inputs, prg, ot_dig) : multiplySendOffline(inputs, prg, ot_dig);
    }
  }
      
  std::vector<Field> OTProviderHA::multiplyRecv(const std::vector<Field>& inputs, fieldDig& ot_dig) {
  size_t num_bits = sizeof(Field) * 8;
    size_t num_blocks = num_bits * inputs.size();

    std::unique_ptr<bool[]> choice_bits(new bool[num_blocks]);
    size_t idx = 0;
    for (auto input : inputs) {
      uint64_t input_64_t = conv<uint64_t>(input);
      for (size_t j = 0; j < num_bits; ++j) {
        bool lsb = ((input_64_t >> j) & 1U) == 1;
        choice_bits[idx++] = lsb;
      }
    }

    std::vector<Field> recv_blocks(num_blocks);
    recv(recv_blocks.data(), choice_bits.get(), num_blocks, ot_dig);

    auto dig = std::vector<Field>{Field(1), Field(2), Field(3), Field(4)};
    auto empty = std::vector<Field>(0); 

    std::vector<Field> shares(inputs.size(), Field(0));
    idx = 0;
    for (size_t i = 0; i < inputs.size(); ++i) {
      for (size_t j = 0; j < num_bits; ++j) {
        shares[i] += recv_blocks[idx];
        idx++;
      }
    }

    return shares;
  }

}; // namespace asyncAsterisk