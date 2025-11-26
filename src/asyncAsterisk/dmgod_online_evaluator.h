#pragma once

#include "dmgod_preproc.h"
#include "rand_gen_pool.h"
#include "netmp.h"
#include "circuit.h"
#include "dmgod_sharing.h"

using namespace io;
using namespace utils;
using namespace asyncAsterisk;

namespace dmAsyncAsteriskGOD {
    typedef struct Message {
        size_t receiver_id;
        std::vector<Field> data;
    } Message;

    class OnlineEvaluator {
        int nP_;
        int id_;
        int security_param_;
        RandGenPool rgen_;
        std::shared_ptr<NetIOMP> network_;
        PreprocCircuit<Field> preproc_;
        LevelOrderedCircuit circ_;
        std::vector<Field> wires_;
        std::vector<Field> q_val_;
        std::vector<TwoShare<Field>> q_sh_;
        std::shared_ptr<ThreadPool> tpool_;
        std::unordered_map<size_t, std::queue<Message>> message_buffer_;
        std::mutex mtx_;
        std::condition_variable cv_;
        std::condition_variable cv_start_recv_;
        static int counter_tp;
        static int counter_p;
        bool start_recv_;

        public:
        OnlineEvaluator(int nP, int id, int security_param, std::shared_ptr<NetIOMP> network, PreprocCircuit<Field> preproc, 
            LevelOrderedCircuit circ, int threads, uint64_t seed = 200);
        ~OnlineEvaluator();

        void setInputs(const std::unordered_map<wire_t, Field> &inputs);
        void setRandomInputs();

        void evaluateGatesAtDepthPartySend(size_t depth, std::vector<Field> &mult_nonTP, std::vector<Field> &mac_components);
        void evaluateGatesAtDepthPartyRecv(size_t depth, std::vector<Field> mult_all);
        void evaluateGatesAtDepth(size_t depth);
        bool MACVerification();
        std::vector<Field> getOutputs();
        std::vector<Field> evaluateCircuit(const std::unordered_map<wire_t, Field> &inputs);
    };
}; // namespace dmAsyncAsterisk