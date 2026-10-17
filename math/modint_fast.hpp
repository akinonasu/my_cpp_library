#pragma once

#include <cassert>
#include <istream>
#include <ostream>
#include <type_traits>

#if defined(__GNUC__) && !defined(__clang__) && (defined(__x86_64__) || defined(__i386__))
#pragma GCC push_options
#pragma GCC optimize("O3")
#endif

namespace fps_detail {
template <class I>
constexpr bool mint_integral = std::is_integral_v<I>
#if defined(__SIZEOF_INT128__)
                               || std::is_same_v<I, __int128> || std::is_same_v<I, unsigned __int128>
#endif
    ;
struct Mint {
    using mint = Mint;
    unsigned v = 0;
    static constexpr unsigned mod() { return 998244353; }
    static unsigned reduce(unsigned long long x) { return (x + (unsigned long long)(unsigned(x) * 998244351u) * mod()) >> 32; }
    Mint() = default;
    template <class I, std::enable_if_t<mint_integral<I>, int> = 0>
    Mint(I x) {
        long long r = x % (long long)mod();
        if (r < 0) r += mod();
        v = reduce((unsigned long long)r * 932051910);
    }
    static Mint raw(unsigned x) {
        Mint a;
        a.v = reduce((unsigned long long)x * 932051910);
        return a;
    }
    static Mint word(unsigned x) {
        Mint a;
        a.v = x;
        return a;
    }
    unsigned val() const {
        unsigned r = reduce(v);
        return r < mod() ? r : r - mod();
    }
    friend std::istream& operator>>(std::istream& is, Mint& x) {
        long long value;
        if (is >> value) x = Mint(value);
        return is;
    }
    friend std::ostream& operator<<(std::ostream& os, const Mint& x) { return os << x.val(); }
    Mint& operator+=(Mint b) {
        v += b.v;
        if (v >= 2 * mod()) v -= 2 * mod();
        return *this;
    }
    Mint& operator-=(Mint b) {
        v -= b.v;
        if (v >> 31) v += 2 * mod();
        return *this;
    }
    Mint& operator*=(Mint b) {
        v = reduce((unsigned long long)v * b.v);
        return *this;
    }
    Mint& operator/=(Mint b) { return *this *= b.inv(); }
    Mint& operator++() { return *this += Mint(1); }
    Mint& operator--() { return *this -= Mint(1); }
    Mint operator++(int) {
        Mint a = *this;
        ++*this;
        return a;
    }
    Mint operator--(int) {
        Mint a = *this;
        --*this;
        return a;
    }
    Mint operator+() const { return *this; }
    Mint operator-() const { return word(v ? 2 * mod() - v : 0); }
    Mint operator+(Mint b) const { return Mint(*this) += b; }
    Mint operator-(Mint b) const { return Mint(*this) -= b; }
    Mint operator*(Mint b) const { return Mint(*this) *= b; }
    Mint operator/(Mint b) const { return Mint(*this) /= b; }
    bool operator==(Mint b) const { return (v < mod() ? v : v - mod()) == (b.v < mod() ? b.v : b.v - mod()); }
    bool operator!=(Mint b) const { return !(*this == b); }
    template <class I, std::enable_if_t<mint_integral<I>, int> = 0>
    friend Mint operator+(I a, Mint b) {
        return Mint(a) + b;
    }
    template <class I, std::enable_if_t<mint_integral<I>, int> = 0>
    friend Mint operator-(I a, Mint b) {
        return Mint(a) - b;
    }
    template <class I, std::enable_if_t<mint_integral<I>, int> = 0>
    friend Mint operator*(I a, Mint b) {
        return Mint(a) * b;
    }
    template <class I, std::enable_if_t<mint_integral<I>, int> = 0>
    friend Mint operator/(I a, Mint b) {
        return Mint(a) / b;
    }
    template <class I, std::enable_if_t<mint_integral<I>, int> = 0>
    friend bool operator==(I a, Mint b) {
        return Mint(a) == b;
    }
    template <class I, std::enable_if_t<mint_integral<I>, int> = 0>
    friend bool operator!=(I a, Mint b) {
        return Mint(a) != b;
    }
    Mint pow(long long n) const {
        assert(n >= 0);
        if (!n) return Mint(1);
        Mint r(1), a = *this;
        for (; n > 1; n >>= 1, a *= a)
            if (n & 1) r *= a;
        return r * a;
    }
    [[gnu::noinline]] Mint inv() const {
        assert(*this != Mint(0));
        Mint a2 = *this * *this, a3 = a2 * *this;
        Mint a7 = a3 * a3 * *this;
        Mint a31 = a7;
        for (int i = 0; i < 2; ++i)
            a31 *= a31;
        a31 *= a3;
        Mint r = a7;
        for (int i = 0; i < 3; ++i)
            r *= r;
        r *= a3;
        r *= r;
        for (int j = 0; j < 4; ++j) {
            for (int i = 0; i < 5; ++i)
                r *= r;
            r *= a31;
        }
        for (int i = 0; i < 3; ++i)
            r *= r;
        return r * a7;
    }
};
}  // namespace fps_detail

using fps_mint = fps_detail::Mint;
using mint = fps_mint;

#if defined(__GNUC__) && !defined(__clang__) && (defined(__x86_64__) || defined(__i386__))
#pragma GCC pop_options
#endif
