#pragma once

#include <emp-ot/emp-ot.h>

#include "helpers.h"

using namespace emp;
using namespace utils;

namespace asyncAsterisk {
    class OTProvider {
        std::array<NetIO*, 1> ios_;
        std::unique_ptr<FerretCOT<NetIO>> ot_;

        public:
        OTProvider(int my_id, int other_id, NetIO* io);

        void sendFieldElements(const Field* data, size_t length);
        void recvFieldElements(Field* data, size_t length);
        void send(const Field* data0, const Field* data1, size_t length);
        void recv(Field* data, const bool* r, size_t length);
        std::vector<Field> multiplySend(const std::vector<Field>& inputs, PRG& prg);
        std::vector<Field> multiplyRecv(const std::vector<Field>& inputs);
    };    
}; // namespace asyncAsterisk