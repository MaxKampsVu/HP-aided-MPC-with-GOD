#pragma once

#include <emp-tool/emp-tool.h>
#include "hm_preproc.h"
#include "rand_gen_pool.h"
#include "netmp.h"
#include "types.h"
#include "circuit.h"

using namespace io;
using namespace utils;
using namespace asyncAsterisk;

namespace hmAsyncAsterisk {
    class OfflineEvaluator {
        int nP_;  
        int th_;
        int id_;
        int security_param_;
        Field key_sh_;
        RandGenPool rgen_;
        std::shared_ptr<NetIOMP> network_;
        LevelOrderedCircuit circ_;
        std::shared_ptr<ThreadPool> tpool_;
        PreprocCircuit<Field> preproc_;

        public:
        OfflineEvaluator(int nP, int id, int security_param, std::shared_ptr<NetIOMP> network, LevelOrderedCircuit circ, 
            int threads, uint64_t seed = 200);

        static PreprocCircuit<Field> dummy(int nP, int id, const LevelOrderedCircuit& circ, 
            const std::unordered_map<wire_t, int>& input_pid_map, PRG& prg);
        void keyGen(std::vector<Field>& keySh, Field& key);
        static void randomShare(int nP, int pid, RandGenPool& rgen, NetIOMP& network, AuthShamirShare<Field>& share, TPShamirShare<Field>& tpShare, Field key,  
            std::vector<Field> keySh, std::vector<std::vector<Field>>& rand_sh, std::vector<size_t>& idx_rand_sh);
        static void randomShareSecret(int nP, int pid, RandGenPool& rgen, NetIOMP& network, AuthShamirShare<Field>& share, TPShamirShare<Field>& tpShare, 
            Field secret, Field key, std::vector<Field> keySh, std::vector<std::vector<Field>>& rand_sh_sec, std::vector<size_t>& idx_rand_sh_sec);
        static void randomShareWithParty(int nP, int pid, int dealer, RandGenPool& rgen, NetIOMP& network, AuthShamirShare<Field>& share, TPShamirShare<Field>& tpShare, 
            Field& secret, Field key, std::vector<Field> keySh, std::vector<std::vector<Field>>& rand_sh_party, std::vector<size_t>& idx_rand_sh_party);
        void setWireMasksParty(const std::unordered_map<wire_t, int>& input_pid_map, Field key, std::vector<Field> keySh, std::vector<std::vector<Field>>& rand_sh, 
            std::vector<std::vector<Field>>& rand_sh_sec, std::vector<std::vector<Field>>& rand_sh_party);
        void setWireMasks(const std::unordered_map<wire_t, int>& input_pid_map);
        PreprocCircuit<Field> run(const std::unordered_map<wire_t, int>& input_pid_map);
    };
}; // namespace hmAsyncAsterisk