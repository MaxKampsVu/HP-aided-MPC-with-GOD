#include "ot_provider.h"

namespace asyncAsterisk {
  constexpr size_t ot_bsize = emp::ot_bsize;

  OTProvider::OTProvider(int my_id, int other_id, NetIO* io) : ios_{io} {

    if (my_id == 0) {
      std::string filename = "./ot_data/s" + std::to_string(other_id) + "_r" + std::to_string(my_id) + "_receiver.bin";
      ot_ = std::make_unique<FerretCOT<NetIO>>(BOB, 1, ios_.data(), true, true, ferret_b13, filename);
    }
    else {
      std::string filename = "./ot_data/s" + std::to_string(my_id) + "_r" + std::to_string(other_id) + "_sender.bin";
      ot_ = std::make_unique<FerretCOT<NetIO>>(ALICE, 1, ios_.data(), true, true, ferret_b13, filename);
    }
  }

  void OTProvider::sendFieldElements(const Field* data, size_t length) {
    std::vector<uint8_t> serialized(length);
    size_t num = (length + FIELDSIZE - 1) / FIELDSIZE;
    for (size_t i = 0; i < num; ++i) {
      BytesFromZZ(serialized.data() + i * FIELDSIZE, conv<ZZ>(data[i]), FIELDSIZE);
    }
    ios_[0]->send_data(serialized.data(), serialized.size());
    ios_[0]->flush();
  }
        
  void OTProvider::recvFieldElements(Field* data, size_t length) {
    std::vector<uint8_t> serialized(length);
    ios_[0]->recv_data(serialized.data(), serialized.size());
    size_t num = (length + FIELDSIZE - 1) / FIELDSIZE;
    for (size_t i = 0; i < num; ++i) {
      data[i] = conv<Field>(ZZFromBytes(serialized.data() + i * FIELDSIZE, FIELDSIZE));
    }
  }

  void OTProvider::send(const Field* data0, const Field* data1, size_t length) {
    auto* data = new emp::block[length];    
    ot_->send_cot(data, length);
    emp::block s;
    ot_->prg.random_block(&s, 1);
    ios_[0]->send_block(&s, 1);
    ot_->mitccrh.setS(s);
    ios_[0]->flush();

    block pad[2 * ot_bsize];
    std::vector<Field> upad(2 * length);
    auto* tpad = reinterpret_cast<uint64_t*>(pad);
    for (size_t i = 0; i < length; i += ot_bsize) {
      for (size_t j = i; j < std::min(i + ot_bsize, length); j++) {
        pad[2 * (j - i)] = data[j];
        pad[2 * (j - i) + 1] = data[j] ^ ot_->Delta;
      }
      ot_->mitccrh.hash<ot_bsize, 2>(pad);
      for (size_t j = i; j < std::min(i + ot_bsize, length); j++) {
        upad[2 * j] = Field(tpad[4 * (j - i)]) + data0[j];
        upad[2 * j + 1] = Field(tpad[4 * (j - i) + 2]) + data1[j];
        
      }
    }
    sendFieldElements(upad.data(), 2 * sizeof(Field) * length);
    delete[] data;
  }

  void OTProvider::recv(Field* rdata, const bool* r, size_t length) {
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
    delete[] data;
  }

  std::vector<Field> OTProvider::multiplySend(const std::vector<Field>& inputs, PRG& prg) {
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

    send(inp_0.data(), inp_1.data(), num_blocks);

    return shares;
  }
      
  std::vector<Field> OTProvider::multiplyRecv(const std::vector<Field>& inputs) {
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
    recv(recv_blocks.data(), choice_bits.get(), num_blocks);

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
