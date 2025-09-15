#pragma once

#include "hm_sharing.h"

namespace hmAsyncAsterisk {
  template <class R>
  struct PreprocGate {
    // Secret shared mask for the output wire of the gate.
    AuthShamirShare<R> mask{};
    TPShamirShare<R> tpmask{};

    PreprocGate() = default;
    explicit PreprocGate(const AuthShamirShare<R>& mask, const TPShamirShare<R>& tpmask) : mask(mask), tpmask(tpmask) {}
    virtual ~PreprocGate() = default;
  };

  template <class R>
  using preprocg_ptr_t = std::unique_ptr<PreprocGate<R>>;

  template <class R>
  struct PreprocInput : public PreprocGate<R> {
    int pid{};
    R mask_value{};

    PreprocInput() = default;
    PreprocInput(const AuthShamirShare<R>& mask, const TPShamirShare<R>& tpmask, int pid, R mask_value = 0) 
        : PreprocGate<R>(mask, tpmask), pid(pid), mask_value(mask_value) {}
    PreprocInput(const PreprocInput<R>& pregate) 
        : PreprocGate<R>(pregate.mask, pregate.tpmask), pid(pregate.pid), mask_value(pregate.mask_value) {}
  };

  template <class R>
  struct PreprocMultGate : public PreprocGate<R> {
    // Secret shared product of inputs masks.
    AuthShamirShare<R> mask_prod{};
    TPShamirShare<R> tpmask_prod{};

    PreprocMultGate() = default;
    PreprocMultGate(const AuthShamirShare<R>& mask, const TPShamirShare<R>& tpmask,
                    const AuthShamirShare<R>& mask_prod, const TPShamirShare<R>& tpmask_prod)
        : PreprocGate<R>(mask, tpmask), mask_prod(mask_prod), tpmask_prod(tpmask_prod) {}
  };

  // Preprocessed data for the circuit.
  template <class R>
  struct PreprocCircuit {
    std::vector<preprocg_ptr_t<R>> gates;

    PreprocCircuit() = default;
    PreprocCircuit(size_t num_gates) : gates(num_gates) {}
  };
};  // namespace hmAsyncAsterisk
