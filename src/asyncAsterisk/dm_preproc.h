#pragma once

#include "dm_sharing.h"

namespace dmAsyncAsterisk {
  template <class R>
  struct PreprocGate {
    // Secret shared mask for the output wire of the gate.
    RepShare<R> mask{};
    R mask_share_zero{};

    PreprocGate() = default;
    explicit PreprocGate(const RepShare<R>& mask) : mask(mask) {}
    explicit PreprocGate(const RepShare<R>& mask, const R& mask_share_zero) : mask(mask), mask_share_zero(mask_share_zero) {}
    virtual ~PreprocGate() = default;
  };

  template <class R>
  using preprocg_ptr_t = std::unique_ptr<PreprocGate<R>>;

  template <class R>
  struct PreprocInput : public PreprocGate<R> {
    int pid{};
    R mask_value{};

    PreprocInput() = default;
    PreprocInput(const RepShare<R>& mask, int pid, R mask_value = 0) 
        : PreprocGate<R>(mask), pid(pid), mask_value(mask_value) {}
    PreprocInput(const PreprocInput<R>& pregate) 
        : PreprocGate<R>(pregate.mask), pid(pregate.pid), mask_value(pregate.mask_value) {}
  };

  template <class R>
  struct PreprocMultGate : public PreprocGate<R> {
    // Secret shared product of inputs masks.
    RepShare<R> mask_prod{};
    RepShare<R> ver{};
    RepShare<R> ver_prod{};

    PreprocMultGate() = default;
    PreprocMultGate(const RepShare<R>& mask, const RepShare<R>& mask_prod)
        : PreprocGate<R>(mask), mask_prod(mask_prod) {}
    PreprocMultGate(const RepShare<R>& mask, const RepShare<R>& mask_prod, const R& mask_share_zero, const RepShare<R>& ver, const RepShare<R>& ver_prod)
        : PreprocGate<R>(mask, mask_share_zero), mask_prod(mask_prod), ver(ver), ver_prod(ver_prod) {}
  };

  // Preprocessed data for the circuit.
  template <class R>
  struct PreprocCircuit {
    std::vector<preprocg_ptr_t<R>> gates;

    PreprocCircuit() = default;
    PreprocCircuit(size_t num_gates) : gates(num_gates) {}
  };
};  // namespace dmAsyncAsterisk
