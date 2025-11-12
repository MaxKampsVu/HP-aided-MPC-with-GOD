#pragma once

#include "helpers.h"

using namespace utils;

namespace dmAsyncAsteriskGOD {
    template <class R>
    class TwoShare {
        R value_;
        R tag_;

        public:
        TwoShare() = default;
        explicit TwoShare(R value, R tag) : value_{value}, tag_{tag} {}

        R getValue() const { return value_; }
        R getTag() { return tag_; }

        void setValue(R value) { value_ = value; }
        void setTag(R tag) { tag_ = tag; }

        TwoShare<R>& operator+=(const TwoShare<R>& rhs) {
            value_ += rhs.value_;
            tag_ += rhs.tag_;
            return *this;
        }

        friend TwoShare<R> operator+(TwoShare<R> lhs, const TwoShare<R>& rhs) {
            lhs += rhs;
            return lhs;
        }

        TwoShare<R>& operator-=(const TwoShare<R>& rhs) {
            (*this) += (rhs * R(-1));
            return *this;
        }

        friend TwoShare<R> operator-(TwoShare<R> lhs, const TwoShare<R>& rhs) {
            lhs -= rhs;
            return lhs;
        }

        TwoShare<R>& operator*=(const TwoShare<R>& rhs) {
            value_ += rhs.value_;
            tag_ += rhs.tag_;
            return *this;
        }

        friend TwoShare<R> operator*(TwoShare<R> lhs, const R& rhs) {
            lhs.value_ *= rhs;
            return lhs;
        }
 
        friend R computeTag(const TwoShare<R>& shareP, const TwoShare<R>& shareTP, R key) {
            return shareP.getValue() * key - shareTP.getValue();
        }

        friend R computeTag(const R& sharePvalue, const TwoShare<R>& shareTP, R key) {
            return sharePvalue * key - shareTP.getValue();
        }
    };

}; // namespace dmAsyncAsteriskGOD