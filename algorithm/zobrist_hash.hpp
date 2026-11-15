#pragma once

#include <bit>
#include <cassert>
#include <chrono>
#include <concepts>
#include <random>
#include <ranges>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <vector>

static_assert(sizeof(int) == 4, "32-bit int is required");
static_assert(sizeof(long) == 8, "64-bit long is required");
static_assert(sizeof(unsigned long) == 8, "64-bit unsigned long is required");

namespace zobrist {

inline unsigned long splitmix64(unsigned long x) {
    x += 0x9e3779b97f4a7c15UL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9UL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebUL;
    return x ^ (x >> 31);
}

inline unsigned long default_seed() {
    static const unsigned long seed = []() {
        std::random_device rd;
        unsigned long t = static_cast<unsigned long>(std::chrono::steady_clock::now().time_since_epoch().count());
        unsigned long r = (static_cast<unsigned long>(rd()) << 32) ^ static_cast<unsigned long>(rd());
        return splitmix64(t ^ r ^ 0x123456789abcdef0UL);
    }();
    return seed;
}

inline unsigned long hash_combine(unsigned long h1, unsigned long h2) {
    return splitmix64(h1 ^ (h2 + 0x9e3779b97f4a7c15UL + (h1 << 6) + (h1 >> 2)));
}

template <typename T>
inline unsigned long hash_value(const T& x, unsigned long seed = default_seed()) {
    using U = std::decay_t<T>;
    if constexpr (std::integral<U>) {
        return splitmix64(static_cast<unsigned long>(x) ^ seed);
    } else if constexpr (std::floating_point<U>) {
        double d = static_cast<double>(x);
        if (d == 0.0) d = 0.0;  // -0.0 と +0.0 を同一視
        return splitmix64(std::bit_cast<unsigned long>(d) ^ seed);
    } else if constexpr (std::convertible_to<U, std::string_view>) {
        std::string_view sv(x);
        unsigned long h = seed;
        for (char c : sv) {
            h = hash_combine(h, splitmix64(static_cast<unsigned long>(c) ^ seed));
        }
        return hash_combine(h, static_cast<unsigned long>(sv.size()));
    } else if constexpr (requires { std::tuple_size<U>::value; } && !std::ranges::input_range<U>) {
        unsigned long h = seed;
        std::apply([&](const auto&... args) { ((h = hash_combine(h, hash_value(args, seed))), ...); }, x);
        return hash_combine(h, static_cast<unsigned long>(std::tuple_size_v<U>));
    } else if constexpr (std::ranges::input_range<U>) {
        unsigned long h = seed;
        unsigned long cnt = 0;
        for (const auto& elem : x) {
            h = hash_combine(h, hash_value(elem, seed));
            ++cnt;
        }
        return hash_combine(h, cnt);
    } else {
        static_assert(sizeof(U) == 0, "Unsupported type for zobrist::hash_value");
        return 0;
    }
}

struct ZobristHasher {
    unsigned long seed;

    explicit ZobristHasher(unsigned long s = default_seed()) : seed(s) {}

    template <typename T>
    unsigned long operator()(const T& x) const {
        return hash_value(x, seed);
    }
};

// 集合のトグル（出現回数の偶奇）を管理するクラス
class ZobristSet {
   public:
    explicit ZobristSet(unsigned long seed = default_seed()) : seed_(seed) {}

    template <typename T>
    void toggle(const T& x) {
        val_ ^= hash_value(x, seed_);
    }

    void clear() { val_ = 0; }
    unsigned long hash() const { return val_; }

    bool operator==(const ZobristSet& other) const {
        assert(seed_ == other.seed_);
        return val_ == other.val_;
    }

   private:
    unsigned long seed_;
    unsigned long val_ = 0;
};

// 多重集合（要素の重複を考慮）を動的に管理するクラス
class ZobristMultiset {
   public:
    explicit ZobristMultiset(unsigned long seed = default_seed()) : seed_(seed) {}

    template <typename T>
    void add(const T& x, long cnt = 1) {
        val_ += hash_value(x, seed_) * static_cast<unsigned long>(cnt);
    }

    template <typename T>
    void remove(const T& x, long cnt = 1) {
        val_ -= hash_value(x, seed_) * static_cast<unsigned long>(cnt);
    }

    void clear() { val_ = 0; }
    unsigned long hash() const { return val_; }

    bool operator==(const ZobristMultiset& other) const {
        assert(seed_ == other.seed_);
        return val_ == other.val_;
    }

   private:
    unsigned long seed_;
    unsigned long val_ = 0;
};

// 静的配列の区間 [l, r) に対する XOR ハッシュ（出現回数の偶奇判定用）
class RangeZobristXor {
   public:
    RangeZobristXor() = default;

    template <std::ranges::input_range R>
    explicit RangeZobristXor(const R& a, unsigned long seed = default_seed()) : seed_(seed) {
        if constexpr (std::ranges::sized_range<R>) {
            pref_.reserve(std::ranges::size(a) + 1);
        }
        pref_.push_back(0);
        for (const auto& x : a) {
            pref_.push_back(pref_.back() ^ hash_value(x, seed_));
        }
        n_ = static_cast<int>(pref_.size()) - 1;
    }

    int size() const { return n_; }

    unsigned long query(int l, int r) const {
        assert(0 <= l && l <= r && r <= n_);
        return pref_[r] ^ pref_[l];
    }

    bool same(int l1, int r1, int l2, int r2) const { return query(l1, r1) == query(l2, r2); }

   private:
    unsigned long seed_ = default_seed();
    int n_ = 0;
    std::vector<unsigned long> pref_{0};
};

// 静的配列の区間 [l, r) に対する多重集合ハッシュ（並び替え一致判定用）
class RangeZobristMultiset {
   public:
    RangeZobristMultiset() = default;

    template <std::ranges::input_range R>
    explicit RangeZobristMultiset(const R& a, unsigned long seed = default_seed()) : seed_(seed) {
        if constexpr (std::ranges::sized_range<R>) {
            pref_.reserve(std::ranges::size(a) + 1);
        }
        pref_.push_back(0);
        for (const auto& x : a) {
            pref_.push_back(pref_.back() + hash_value(x, seed_));
        }
        n_ = static_cast<int>(pref_.size()) - 1;
    }

    int size() const { return n_; }

    unsigned long query(int l, int r) const {
        assert(0 <= l && l <= r && r <= n_);
        return pref_[r] - pref_[l];
    }

    bool same(int l1, int r1, int l2, int r2) const { return query(l1, r1) == query(l2, r2); }

   private:
    unsigned long seed_ = default_seed();
    int n_ = 0;
    std::vector<unsigned long> pref_{0};
};

}  // namespace zobrist
