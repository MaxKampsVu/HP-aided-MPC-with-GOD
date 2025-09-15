#pragma once

#include "helpers.h"

using namespace utils;

namespace dmAsyncAsterisk {
    template <class R>
    class RepShare {
        std::vector<R> key_sh_;
        std::vector<R> values_;
        std::vector<R> tags_;

        public:
        RepShare() = default;
        explicit RepShare(std::vector<R> key_sh, std::vector<R> values, std::vector<R> tags) : key_sh_{key_sh}, values_{values}, tags_{tags} {}

        std::vector<R> getKeySh() const { return key_sh_; }
        std::vector<R> getValues() const { return values_; }
        std::vector<R> getTags() { return tags_; }

        void setKeySh(std::vector<R> key_sh) { key_sh_ = key_sh; }
        void setValues(std::vector<R> values) { values_ = values; }
        void setTags(std::vector<R> tags) { tags_ = tags; }

        void pushValues(R val) { values_.push_back(val); }
        void pushTags(R tag) { tags_.push_back(tag);}
        void setKeySh(R keysh) { key_sh_.push_back(keysh); }

        RepShare<R>& operator+=(const RepShare<R>& rhs) {
            for (size_t i = 0; i < values_.size(); i++) {
                values_[i] += rhs.values_[i];
                key_sh_[i] = rhs.key_sh_[i];
            }
            for (size_t i = 0; i < tags_.size(); i++) {
                tags_[i] += rhs.tags_[i];
            }
            return *this;
        }

        friend RepShare<R> operator+(RepShare<R> lhs, const RepShare<R>& rhs) {
            lhs += rhs;
            return lhs;
        }

        RepShare<R>& operator-=(const RepShare<R>& rhs) {
            (*this) += (rhs * R(-1));
            return *this;
        }

        friend RepShare<R> operator-(RepShare<R> lhs, const RepShare<R>& rhs) {
            lhs -= rhs;
            return lhs;
        }

        RepShare<R>& operator*=(const R& rhs) {
            for(size_t i = 0; i < values_.size(); i++) {
                values_[i] *= rhs;
            }
            for(size_t i = 0; i < tags_.size(); i++) {
                tags_[i] *= rhs;
            }
            return *this;
        }

        friend RepShare<R> operator*(RepShare<R> lhs, const R& rhs) {
            lhs *= rhs;
            return lhs;
        }

        RepShare<R>& add(R val, int id) {
            if (id!=0) {
                values_[1] += val;
                for (size_t i = 0; i < tags_.size(); i++) {
                    tags_[i] += key_sh_[i]*val;
                }
            }
            else {
                for (size_t i = 0; i < tags_.size(); i++) {
                    tags_[i] += key_sh_[i]*val;
                }
            }
            return *this;
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
        explicit DummyShare(int nP) : key_sh_(nP+1), values_(nP+1), tags_(nP+1) {}

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
            for (size_t i = 0; i < values_.size(); i++) {
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
            for(size_t i = 0; i < values_.size(); i++) {
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

        RepShare<R> getRSS(int id) {
            return RepShare<R>(getKeyShare(id), getValShare(id), getTagShare(id));
        }

        std::vector<R> getKeyShare(int id) {
            std::vector<R> key_sh;
            if (id == 0) {
                for (size_t i = 0; i < key_sh_.size()-1; i++) {
                    key_sh.push_back(key_sh_[i]);
                }
            }
            else {
                key_sh.push_back(key_sh_[id-1]);
                key_sh.push_back(key_sh_[key_sh_.size()-1]);
            }
            return key_sh;
        }

        std::vector<R> getValShare(int id) {
            std::vector<R> val_sh;
            if (id == 0) {
                for (size_t i = 0; i < values_.size()-1; i++) {
                    val_sh.push_back(values_[i]);
                }
            }
            else {
                val_sh.push_back(values_[id-1]);
                val_sh.push_back(values_[values_.size()-1]);
            }
            return val_sh;
        }

        std::vector<R> getTagShare(int id) {
            std::vector<R> tag_sh;
            if (id == 0) {
                for (size_t i = 0; i < tags_.size()-1; i++) {
                    tag_sh.push_back(tags_[i]);
                }
            }
            else {
                tag_sh.push_back(tags_[id-1]);
                tag_sh.push_back(tags_[tags_.size()-1]);
            }
            return tag_sh;
        }
    };
}; // namespace dmAsyncAsterisk