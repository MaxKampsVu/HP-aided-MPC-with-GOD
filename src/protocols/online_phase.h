#pragma once

#include "preproc.h"
#include "rand_gen_pool.h"
#include "netmp.h"
#include "circuit.h"
#include "sharing.h"
#include <set>
#include <queue>
#include <mutex>
#include <condition_variable>

using namespace io;
using namespace utils;
using namespace dmGOD;

namespace dmGOD {
    typedef struct Message {
        int receiver_id; // Changed to int to match party IDs
        std::vector<Field> data;
    } Message;

    class OnlineEvaluator {
        int nP_;
        int id_;
        int security_param_;
        int num_cheaters_;
        std::set<int> excluded_parties_;
        RandGenPool rgen_;
        std::shared_ptr<NetIOMP> network_;
        PreprocCircuit<Field> preproc_;
        LevelOrderedCircuit circ_;
        std::vector<Field> wires_;
        std::vector<Field> q_val_;
        std::vector<TwoShare<Field>> q_sh_;
        std::shared_ptr<ThreadPool> tpool_;

        // Storage for messages received by Party 0
        std::unordered_map<size_t, std::queue<Message>> message_buffer_;

        // Set to track parties that failed MAC verification

        std::mutex mtx_;
        std::condition_variable cv_;
        std::condition_variable cv_start_recv_;
        static int counter_tp;
        static int counter_p;
        bool start_recv_;

    public:
        // Updated constructor signature to include num_cheaters
        OnlineEvaluator(int nP, int id, int security_param, std::shared_ptr<NetIOMP> network,
                        PreprocCircuit<Field> preproc, LevelOrderedCircuit circ,
                        int threads, uint64_t seed = 200, int num_cheaters = 0);
        ~OnlineEvaluator();

        void setInputs(const std::unordered_map<wire_t, Field> &inputs);
        void setRandomInputs();

        void evaluateGatesAtDepthPartySend(size_t depth, std::vector<Field> &mult_nonTP, std::vector<Field> &mac_components);
        void evaluateGatesAtDepthPartyRecv(size_t depth, std::vector<Field> mult_all);
        void evaluateGatesAtDepth(size_t depth);

        // MAC verification logic is now integrated into evaluateGatesAtDepth
        // per the logic in the .cpp file, but kept if you need external calls.
        bool MACVerification();

        std::vector<Field> getOutputs();
        std::vector<Field> evaluateCircuit(const std::unordered_map<wire_t, Field> &inputs);
    };
}; // namespace dmGOD