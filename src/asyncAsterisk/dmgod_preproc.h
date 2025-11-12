#pragma once

#include "dmgod_sharing.h"

namespace dmAsyncAsteriskGOD {
  template <class R>
  struct PreprocGate {
    // Secret shared mask for the output wire of the gate.
    TwoShare<R> mask{};
    R mask_share_zero{};

    PreprocGate() = default;
    explicit PreprocGate(const TwoShare<R>& mask) : mask(mask) {}
    explicit PreprocGate(const TwoShare<R>& mask, const R& mask_share_zero) : mask(mask), mask_share_zero(mask_share_zero) {}
    virtual ~PreprocGate() = default;
  };

  template <class R>
  using preprocg_ptr_t = std::unique_ptr<PreprocGate<R>>;

  template <class R>
  struct PreprocInput : public PreprocGate<R> {
    int pid{};
    R mask_value{};

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
    TwoShare<R> ver{};
    TwoShare<R> ver_prod{};

    PreprocMultGate() = default;
    PreprocMultGate(const TwoShare<R>& mask, const TwoShare<R>& mask_prod)
        : PreprocGate<R>(mask), mask_prod(mask_prod) {}
    PreprocMultGate(const TwoShare<R>& mask, const TwoShare<R>& mask_prod, const R& mask_share_zero, const TwoShare<R>& ver, const TwoShare<R>& ver_prod)
        : PreprocGate<R>(mask, mask_share_zero), mask_prod(mask_prod), ver(ver), ver_prod(ver_prod) {}
  };

  // Preprocessed data for the circuit.
  template <class R>
  struct PreprocCircuit {
    std::vector<preprocg_ptr_t<R>> gates;

    PreprocCircuit() = default;
    PreprocCircuit(size_t num_gates) : gates(num_gates) {}
  };
  
};  // namespace dmAsyncAsterisksGOD
