#pragma once

#include <emp-tool/emp-tool.h>

using namespace emp;

namespace asyncAsterisk {
  // Collection of PRGs.
  class RandGenPool {
    int id_;

    PRG k_p0;
    PRG k_self;
    PRG k_all_minus_0;
    PRG k_all;
    std::vector<PRG> k_pi;
    std::vector<std::vector<PRG>> k_pij;
    std::vector<PRG> k_pj;
    

  public:
    explicit RandGenPool(int my_id, int num_parties, uint64_t seed = 200);
    
    PRG& self();
    PRG& all_minus_0();
    PRG& all();
    PRG& p0();
    PRG& pi(int i);
    PRG& pij(int i, int j);
    PRG& pj(int j);
  };
};  // namespace asyncAsterisk
