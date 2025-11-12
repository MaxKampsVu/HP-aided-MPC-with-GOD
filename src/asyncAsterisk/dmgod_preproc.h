#pragma once

#include "dmgod_sharing.h"

namespace dmAsyncAsteriskGOD {
  template <class R>
  struct PreprocGate {
    TwoShare<R> mask; // Secret shared mask for the output wire of the gate.

    PreprocGate() = default;
    explicit PreprocGate(const TwoShare<R>& mask) : mask(mask) {}
    virtual ~PreprocGate() = default;
  };

  template <class R>
  using preprocg_ptr_t = std::unique_ptr<PreprocGate<R>>;

  template <class R>
  struct PreprocInput : public PreprocGate<R> {
    int pid{}; // id of input provider 
    R mask_value{}; // mask known only to input provider 

    PreprocInput() = default;
    PreprocInput(const TwoShare<R>& mask, int pid, R mask_value = 0) 
        : PreprocGate<R>(mask), pid(pid), mask_value(mask_value) {}
    PreprocInput(const PreprocInput<R>& pregate) 
        : PreprocGate<R>(pregate.mask), pid(pregate.pid), mask_value(pregate.mask_value) {}
  };

  template <class R>
  struct PreprocMultGate : public PreprocGate<R> {
    // Secret shared product of inputs masks.
    TwoShare<R> mask_prod{};

    PreprocMultGate() = default;
    PreprocMultGate(const TwoShare<R>& mask, const TwoShare<R>& mask_prod)
        : PreprocGate<R>(mask), mask_prod(mask_prod) {}
  };

  // Preprocessed data for the circuit.
  template <class R>
  struct PreprocCircuit {
    std::vector<preprocg_ptr_t<R>> gates;

    PreprocCircuit() = default;
    PreprocCircuit(size_t num_gates) : gates(num_gates) {}
  };
  
};  // namespace dmAsyncAsterisksGOD
