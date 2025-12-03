#pragma once

#include "helpers.h"

using namespace utils;

namespace dmAsyncAsteriskGOD {
    template <class R>
    class TwoShare {
        R value_;
        R mac_component_;

        public:
        TwoShare() = default;
        explicit TwoShare(R value) : value_{value} {}

        R getValue() const { return value_; }
        R getMACComponent() { return mac_component_; }

        void setValue(R value) { value_ = value; }
        void setMACComponent(R mac_component) { mac_component_ = mac_component; }

        TwoShare<R>& add(R val, int id) {
            if (id==0) {
                value_ += val;
            }
            return *this;
        }

        TwoShare<R>& operator+=(const TwoShare<R>& rhs) {
            value_ += rhs.value_;
            mac_component_ += rhs.mac_component_;
            return *this;
        }

        TwoShare<R>& operator-=(const TwoShare<R>& rhs) {
            (*this) += (rhs * R(-1));
            return *this;
        }

        friend TwoShare<R> operator+(TwoShare<R> lhs, const TwoShare<R>& rhs) {
            lhs += rhs;
            return lhs;
        }

        friend TwoShare<R> operator-(TwoShare<R> lhs, const TwoShare<R>& rhs) {
            lhs -= rhs;
            return lhs;
        }

        friend TwoShare<R> operator+(TwoShare<R> lhs, const R& rhs) {
            lhs.value_ += rhs;
            lhs.mac_component_ += rhs;
            return lhs;
        }

        // friend TwoShare<R> operator-(TwoShare<R> lhs, const R& rhs) {
        //     lhs.value_ -= rhs;
        //     lhs.mac_component_ -= rhs;
        //     return lhs;
        // }

        // friend TwoShare<R> operator-(const R& lhs, TwoShare<R> rhs) {
        //     return rhs - lhs;
        // }

        friend TwoShare<R> operator*(TwoShare<R> lhs, const R& rhs) {
            lhs.value_ *= rhs;
            lhs.mac_component_ *= rhs;
            return lhs;
        }

    };

    template <typename R>
    void computeMACs(R key, std::vector<R>& tp_mac_components, const std::vector<R>& p_shares) {
        for (int i = 0; i < tp_mac_components.size(); i++) {
            tp_mac_components[i] = key * p_shares[i] - tp_mac_components[i];
        }
    }
}; // namespace dmAsyncAsteriskGOD