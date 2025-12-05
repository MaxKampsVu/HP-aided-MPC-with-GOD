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
    R tp_key;

    PreprocCircuit() = default;
    PreprocCircuit(size_t num_gates) : gates(num_gates){}

    void setTPKey(Field key) {
      tp_key = key;
    }
  };

  template <class R>
  struct PreprocDotpGate : public PreprocGate<R> {
    TwoShare<R> mask_prod{};

    PreprocDotpGate() = default;
    PreprocDotpGate(const TwoShare<R>& mask, const TwoShare<R>& mask_prod)
        : PreprocGate<R>(mask), mask_prod(mask_prod) {}
  };

  template <class R>
struct PreprocMult3Gate : public PreprocGate<R> {
  // Secret shared product of inputs masks.
  TwoShare<R> mask_ab{};
  TwoShare<R> mask_ac{};
  TwoShare<R> mask_bc{};
  TwoShare<R> mask_abc{};

  PreprocMult3Gate() = default;

  PreprocMult3Gate(const TwoShare<R>& mask,
                  const TwoShare<R>& mask_ab,
                  const TwoShare<R>& mask_ac,
                  const TwoShare<R>& mask_bc)
      : PreprocGate<R>(mask), mask_ab(mask_ab),
                              mask_ac(mask_ac),
                              mask_bc(mask_bc) {}

  void setLength2Terms(TwoShare<R> m_ab, TwoShare<R> m_ac, TwoShare<R> m_bc) {
    mask_ab = m_ab;
    mask_ac = m_ac;
    mask_bc = m_bc;
  }
  
  void setLength3Terms(TwoShare<R> m_abc) {
    mask_abc = m_abc;
  }

  PreprocMult3Gate(const TwoShare<R>& mask,
                  const TwoShare<R>& mask_ab,
                  const TwoShare<R>& mask_ac,
                  const TwoShare<R>& mask_bc,
                  const TwoShare<R>& mask_abc)
      : PreprocGate<R>(mask), mask_ab(mask_ab),
                              mask_ac(mask_ac),
                              mask_bc(mask_bc),
                              mask_abc(mask_abc){}
};


template <class R>
struct PreprocMult4Gate : public PreprocGate<R> {
  // Secret shared product of inputs masks.
  TwoShare<R> mask_abcd{};
  TwoShare<R> mask_abc{};
  TwoShare<R> mask_abd{};
  TwoShare<R> mask_acd{};
  TwoShare<R> mask_bcd{};
  TwoShare<R> mask_ab{};
  TwoShare<R> mask_ac{};
  TwoShare<R> mask_ad{};
  TwoShare<R> mask_bc{};
  TwoShare<R> mask_bd{};
  TwoShare<R> mask_cd{};


  PreprocMult4Gate() = default;
  PreprocMult4Gate(const TwoShare<R>& mask,
                  const TwoShare<R>& mask_ab,
                  const TwoShare<R>& mask_ac,
                  const TwoShare<R>& mask_ad,
                  const TwoShare<R>& mask_bc,
                  const TwoShare<R>& mask_bd,
                  const TwoShare<R>& mask_cd,
                  const TwoShare<R>& mask_abc,
                  const TwoShare<R>& mask_abd,
                  const TwoShare<R>& mask_acd,
                  const TwoShare<R>& mask_bcd,
                  const TwoShare<R>& mask_abcd)
      : PreprocGate<R>(mask), mask_ab(mask_ab),
                              mask_ac(mask_ac),
                              mask_ad(mask_ad),
                              mask_bc(mask_bc),
                              mask_bd(mask_bd),
                              mask_cd(mask_cd), 
                              mask_abc(mask_abc),
                              mask_abd(mask_abd),
                              mask_acd(mask_acd),
                              mask_bcd(mask_bcd),
                              mask_abcd(mask_abcd) {}
};


};  // namespace dmAsyncAsterisksGOD
