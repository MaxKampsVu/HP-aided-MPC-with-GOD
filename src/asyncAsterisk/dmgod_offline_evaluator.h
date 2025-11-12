#pragma once

#include <emp-tool/emp-tool.h>
#include "dmgod_preproc.h"
#include "rand_gen_pool.h"
#include "ot_provider.h"
#include "netmp.h"
#include "types.h"
#include "circuit.h"

using namespace io;
using namespace utils;
using namespace asyncAsterisk;

namespace dmAsyncAsteriskGOD {
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

        // static PreprocCircuit<Field> dummy(int nP, int id, const LevelOrderedCircuit& circ, 
        //     const std::unordered_map<wire_t, int>& input_pid_map, PRG& prg);
        void keyGen();
        static void randomShare(int pid, RandGenPool& rgen, TwoShare<Field>& share, Field& mask_share_zero, bool isOutputWire);
        static void randomShareSecret(int pid, RandGenPool& rgen, const TwoShare<Field>& share1, const TwoShare<Field>& share2, 
            TwoShare<Field>& prodShare, std::vector<Field>& inputToOPE);
        static void randomShareWithParty(int pid, int dealer, RandGenPool& rgen, TwoShare<Field>& share, Field& secret);
        void runOPE(std::vector<Field>& inputToOPE, std::vector <Field>& outputOfOPE, size_t count);
        void multSS(const Field& share1, const Field& share2, Field& output, const std::vector<Field>& outputOfOPE, size_t& idx_outputOfOPE);
        void prepareMaskValues(const std::unordered_map<wire_t,int>& input_pid_map);
        // void prepareMaskTags();
        void setWireMasks(const std::unordered_map<wire_t, int>& input_pid_map);
        PreprocCircuit<Field> run(const std::unordered_map<wire_t, int>& input_pid_map);
    };
}; // namespace dmAsyncAsteriskGOD
