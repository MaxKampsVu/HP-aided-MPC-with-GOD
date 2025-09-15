#include "rand_gen_pool.h"

namespace asyncAsterisk {
  RandGenPool::RandGenPool(int my_id, int num_parties,  uint64_t seed) : id_{my_id}, k_pi(num_parties), k_pij(num_parties), k_pj(num_parties) { 
    auto seed_block = makeBlock(seed, 0); 
    k_self.reseed(&seed_block, 0);
    k_all.reseed(&seed_block, 0);
    k_all_minus_0.reseed(&seed_block, 0);
    k_p0.reseed(&seed_block, 0);
    for(int i = 1; i <= num_parties; i++) {k_pi[i-1].reseed(&seed_block, 0);}
    for(int i = 1; i <= num_parties; i++) {
      k_pij[i-1].resize(num_parties);
      for(int j = 1; j <= num_parties; j++) {
        if (i==j)
          k_pij[i-1][j-1] = k_pi[i];
        else
          k_pij[i-1][j-1].reseed(&seed_block, 0);
      }
    }
    for(int j = 1; j <= num_parties; j++) {
      if (id_ == j)
        k_pj[j-1] = k_p0;
      else
        k_pj[j-1].reseed(&seed_block, 0);
    }
  }
  //all keys will be the same.  for different keys look at emp toolkit

  PRG& RandGenPool::self() { return k_self; }
  PRG& RandGenPool::all() { return k_all; }
  PRG& RandGenPool::all_minus_0() { return k_all_minus_0; }
  PRG& RandGenPool::p0() {  return k_p0; }
  PRG& RandGenPool::pi(int i) { return k_pi[i-1]; }
  PRG& RandGenPool::pij(int i, int j) { return k_pij[i-1][j-1]; }
  PRG& RandGenPool::pj(int j) { return k_pj[j-1]; }
};  // namespace asyncAsterisk
