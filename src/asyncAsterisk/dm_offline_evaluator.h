#pragma once

#include <emp-tool/emp-tool.h>
#include "dm_preproc.h"
#include "rand_gen_pool.h"
#include "ot_provider.h"
#include "netmp.h"
#include "types.h"
#include "circuit.h"

using namespace io;
using namespace utils;
using namespace asyncAsterisk;

namespace dmAsyncAsterisk {
    typedef struct Offline_Message {
        size_t receiver_id;
        std::vector<Field> data;
    } Offline_Message;

    class OfflineEvaluator {
        int nP_;  
        int id_;
        int security_param_;
        std::vector<Field> key_sh_;
        RandGenPool rgen_;
        std::shared_ptr<NetIOMP> network_;
        std::shared_ptr<NetIOMP> network_ot_;
        LevelOrderedCircuit circ_;
        std::shared_ptr<ThreadPool> tpool_;
        PreprocCircuit<Field> preproc_;
        std::vector<std::unique_ptr<OTProvider>> ot_;
        std::unordered_map<size_t, std::queue<Offline_Message>> offline_message_buffer_;
        std::size_t chunk_size_;
        std::mutex mtx_;
        std::condition_variable cv_;
        std::array<std::condition_variable, 2> cv_start_ot_;
        std::vector<bool> start_ot_;
        std::vector<std::vector<Field>> inputToOPE;

        public:
        OfflineEvaluator(int nP, int id, int security_param, std::shared_ptr<NetIOMP> network1, std::shared_ptr<NetIOMP> network2, 
            LevelOrderedCircuit circ, int threads, uint64_t seed = 200);
        ~OfflineEvaluator();

        static PreprocCircuit<Field> dummy(int nP, int id, const LevelOrderedCircuit& circ, 
            const std::unordered_map<wire_t, int>& input_pid_map, PRG& prg);
        void keyGen();
        static void randomShare(int nP, int pid, RandGenPool& rgen, RepShare<Field>& share, const std::vector<Field>& keySh, 
            Field& mask_share_zero, bool isOutputWire);
        static void randomShareSecret(int nP, int pid, RandGenPool& rgen, const RepShare<Field>& share1, const RepShare<Field>& share2, 
            RepShare<Field>& prodShare, const std::vector<Field>& keySh, std::vector<Field>& inputToOPE);
        static void randomShareWithParty(int nP, int pid, int dealer, RandGenPool& rgen, RepShare<Field>& share, Field& secret, 
            const std::vector<Field>& keySh);
        void runOPE(std::vector<Field>& inputToOPE, std::vector <Field>& outputOfOPE, size_t count);
        void multiply(const std::vector<Field>& vec_a, const std::vector<Field>& vec_b, std::vector<Field>& vec_c, 
            const std::vector<Field>& outputOfOPE, std::vector<Field>& buffer, size_t& idx_outputOfOPE, size_t& idx_buffer);
        void prepareMaskValues(const std::unordered_map<wire_t,int>& input_pid_map);
        void prepareMaskTags();
        void setWireMasks(const std::unordered_map<wire_t, int>& input_pid_map);
        bool TripleSacrifice();
        PreprocCircuit<Field> run(const std::unordered_map<wire_t, int>& input_pid_map);
    };
}; // namespace dmAsyncAsterisk
