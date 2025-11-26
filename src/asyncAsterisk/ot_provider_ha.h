#pragma once

#include <emp-ot/emp-ot.h>

#include "helpers.h"

using namespace emp;
using namespace utils;

namespace asyncAsterisk {
    // Forward declaration
    class OTProviderHAWrapper;

    class OTProviderHA {
        friend class OTProviderHAWrapper;
        std::array<NetIO*, 1> ios_;
        std::unique_ptr<FerretCOT<NetIO>> ot_;

        public:
        OTProviderHA(int my_id, int other_id, NetIO* io);

        void sendFieldElements(const Field* data, size_t length);
        void recvFieldElements(Field* data, size_t length);
        void send(const Field* data0, const Field* data1, size_t length, PRG& prg, fieldDig& ot_dig);
        void recv(Field* data, const bool* r, size_t length, fieldDig& ot_dig);
        std::vector<Field> multiplySend(const std::vector<Field>& inputs, PRG& prg, fieldDig& ot_dig);
        std::vector<Field> multiplySendOffline(const std::vector<Field> inputs, PRG& prg, fieldDig& ot_dig);
        std::vector<Field> multiplyRecv(const std::vector<Field>& inputs, fieldDig& ot_dig);
    }; 
}; // namespace asyncAsterisk