#pragma once

#include "helpers.h"

using namespace utils;

namespace hmAsyncAsterisk {
    template <class R>
    class AuthShamirShare {
        R key_sh_;
        R value_;
        R tag_;

        public:
        AuthShamirShare() = default;
        explicit AuthShamirShare(R key_sh, R value, R tag) : key_sh_{key_sh}, value_{value}, tag_{tag} {}

        R& valueAt() { return value_; }
        R& tagAt() { return tag_; }
        R& keySh() { return key_sh_; }

        void pushValue(R val) { value_ = val; }
        void pushTag(R tag) { tag_ = tag; }
        void setKey(R key) { key_sh_ = key; }
  
        R valueAt() const { return value_; }
        R tagAt() const { return tag_; }
        R keySh() const { return key_sh_; }

        AuthShamirShare<R>& operator+=(const AuthShamirShare<R>& rhs) {
            value_ += rhs.value_;
            tag_ += rhs.tag_;
            key_sh_ = rhs.key_sh_;
            return *this;
        }

        friend AuthShamirShare<R> operator+(AuthShamirShare<R> lhs, const AuthShamirShare<R>& rhs) {
            lhs += rhs;
            return lhs;
        }

        AuthShamirShare<R>& operator-=(const AuthShamirShare<R>& rhs) {
            (*this) += (rhs * R(-1));
            return *this;
        }

        friend AuthShamirShare<R> operator-(AuthShamirShare<R> lhs, const AuthShamirShare<R>& rhs) {
            lhs -= rhs;
            return lhs;
        }

        AuthShamirShare<R>& operator*=(const R& rhs) {
            value_ *= rhs;
            tag_ *= rhs;
            return *this;
        }

        friend AuthShamirShare<R> operator*(AuthShamirShare<R> lhs, const R& rhs) {
            lhs *= rhs;
            return lhs;
        }

        AuthShamirShare<R>& add(R val) {
            value_ += val;
            tag_ += key_sh_*val;
            return *this;
        }
    };

    template <class R>
    class TPShamirShare {
        R key_;
        std::vector<R> key_sh_;
        std::vector<R> values_;
        std::vector<R> tags_;

        public:
        TPShamirShare() = default;
        explicit TPShamirShare(R key, std::vector<R> key_sh, std::vector<R> values, std::vector<R> tags) : key_{key}, key_sh_{key_sh}, values_{values}, tags_{tags} {}

        R& operator[](size_t idx) { return values_.at(idx); }
        R operator[](size_t idx) const { return values_.at(idx); }
        R& macKey() { return key_; }
        R& commonValueWithParty(int pid) { return values_.at(pid); }
        R& commonTagWithParty(int pid)  { return tags_.at(pid); }
        R& commonKeyWithParty(int pid) { return key_sh_.at(pid); }
        [[nodiscard]] R commonValueWithParty(int pid) const { return values_.at(pid);}
        [[nodiscard]] R commonTagWithParty(int pid) const { return tags_.at(pid); }
        [[nodiscard]] R commonKeyWithParty(int pid) const { return key_sh_.at(pid);}
        [[nodiscard]] R secret() const { 
            int t = (values_.size() - 2) / 2;
            std::vector<Field> x(t + 1);
            std::vector<Field> y(t + 1);
            for (int i = 0; i <= t; ++i) {
                x[i] = Field(i + 1); 
                y[i] = values_[i + 1]; 
            }
            fieldPoly poly = reconstructPolynomial(x, y);
            return coeff(poly, 0);
        }

        void setKey(R key) { key_ = key;}
        void pushValues(R val) { values_.push_back(val); }
        void pushTags(R tag) { tags_.push_back(tag);}
        void setKeySh(R keysh) { key_sh_.push_back(keysh); }

        TPShamirShare<R>& operator+=(const TPShamirShare<R>& rhs) {
            for (size_t i = 1; i < values_.size(); i++) {
                values_[i] += rhs.values_[i];
                tags_[i] += rhs.tags_[i];
                key_sh_[i] = rhs.key_sh_[i];
            }
            key_ = rhs.key_;
            return *this;
        }

        friend TPShamirShare<R> operator+(TPShamirShare<R> lhs, const TPShamirShare<R>& rhs) {
            lhs += rhs;
            return lhs;
        }

        TPShamirShare<R>& operator-=(const TPShamirShare<R>& rhs) {
            (*this) += (rhs * R(-1));
            return *this;
        }

        friend TPShamirShare<R> operator-(TPShamirShare<R> lhs, const TPShamirShare<R>& rhs) {
            lhs -= rhs;
            return lhs;
        }

        TPShamirShare<R>& operator*=(const R& rhs) {
            for(size_t i = 1; i < values_.size(); i++) {
                values_[i] *= rhs;
                tags_[i] *= rhs;
            }
            return *this;
        }

        friend TPShamirShare<R> operator*(TPShamirShare<R> lhs, const R& rhs) {
            lhs *= rhs;
            return lhs;
        }

        AuthShamirShare<R> getAAS(size_t pid) {
            return AuthShamirShare<R>({key_sh_.at(pid), values_.at(pid), tags_.at(pid)});
        }
    };

    template <class R>
    class DummyShare {
        std::vector<R> key_sh_;
        std::vector<R> values_;
        std::vector<R> tags_;
        R key_;
        R value_;
        R tag_;

        public:
        DummyShare() = default;
        explicit DummyShare(int nP) : key_sh_(nP), values_(nP), tags_(nP) {}

        std::vector<R> getKeySh() const { return key_sh_; }
        std::vector<R> getValues() const { return values_; }
        std::vector<R> getTags() { return tags_; }

        void setKeySh(std::vector<R> key_sh) { key_sh_ = key_sh; }
        void setValues(std::vector<R> values) { values_ = values; }
        void setTags(std::vector<R> tags) { tags_ = tags; }

        void pushValues(R val) { values_.push_back(val); }
        void pushTags(R tag) { tags_.push_back(tag);}
        void setKeySh(R keysh) { key_sh_.push_back(keysh); }

        void setKey(R key) { key_ = key; }
        void setValue(R value) { value_ = value; }
        void setTag(R tag) { tag_ = tag; }

        R secretKey() const {
            return key_;
        }

        R secretVal() const {
            return value_;
        }

        R secretTag() const {
            return tag_;
        }

        DummyShare<R>& operator+=(const DummyShare<R>& rhs) {
            for (size_t i = 1; i < values_.size(); i++) {
                values_[i] += rhs.values_[i];
                tags_[i] += rhs.tags_[i];
                key_sh_[i] = rhs.key_sh_[i];
            }
            key_ = rhs.key_;
            value_ += rhs.value_;
            tag_ += rhs.tag_;
            return *this;
        }

        friend DummyShare<R> operator+(DummyShare<R> lhs, const DummyShare<R>& rhs) {
            lhs += rhs;
            return lhs;
        }

        DummyShare<R>& operator-=(const DummyShare<R>& rhs) {
            (*this) += (rhs * R(-1));
            return *this;
        }

        friend DummyShare<R> operator-(DummyShare<R> lhs, const DummyShare<R>& rhs) {
            lhs -= rhs;
            return lhs;
        }

        DummyShare<R>& operator*=(const R& rhs) {
            for(size_t i = 1; i < values_.size(); i++) {
                values_[i] *= rhs;
                tags_[i] *= rhs;
            }
            value_ *= rhs;
            tag_ *= rhs;
            return *this;
        }

        friend DummyShare<R> operator*(DummyShare<R> lhs, const R& rhs) {
            lhs *= rhs;
            return lhs;
        }

        AuthShamirShare<R> getASS(int id) {
            if (id == 0) {
                return AuthShamirShare<R>(secretKey(), Field(0), Field(0));
            }
            else {
                return AuthShamirShare<R>(key_sh_[id-1], values_[id-1], tags_[id-1]);
            }            
        }

        TPShamirShare<R> getTSS(int id) {
            if (id == 0) {
                return TPShamirShare<R>(secretKey(), key_sh_, values_, tags_);
            }
            else {
                TPShamirShare<R> tss;
                return tss;
            }            
        }
    };
}; // namespace hmAsyncAsterisk