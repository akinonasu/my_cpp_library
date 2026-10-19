#include <algorithm>
#include <bit>
#include <cassert>
#include <concepts>
#include <iostream>
#include <limits>
#include <numeric>
#include <print>
#include <vector>

// --- 基本セグ木実装 ---
template <class S, auto Op, auto E>
    requires std::invocable<decltype(Op), S, S> && std::same_as<std::invoke_result_t<decltype(Op), S, S>, S> &&
             std::invocable<decltype(E)> && std::same_as<std::invoke_result_t<decltype(E)>, S>
class SegTree {
   private:
    int _n;
    int size;
    std::vector<S> d;

    void update(int k) { d[k] = Op(d[2 * k], d[2 * k + 1]); }

   public:
    SegTree() : SegTree(0) {}
    explicit SegTree(int n) : SegTree(std::vector<S>(n, E())) {}
    explicit SegTree(const std::vector<S>& v) : _n(static_cast<int>(v.size())) {
        size = std::bit_ceil(static_cast<unsigned int>(_n));
        d.assign(2 * size, E());
        for (int i = 0; i < _n; i++)
            d[size + i] = v[i];
        for (int i = size - 1; i >= 1; i--)
            update(i);
    }

    void set(int p, S x) {
        assert(0 <= p && p < _n);
        p += size;
        d[p] = x;
        for (p >>= 1; p >= 1; p >>= 1)
            update(p);
    }

    S get(int p) const {
        assert(0 <= p && p < _n);
        return d[p + size];
    }

    S prod(int l, int r) const {
        assert(0 <= l && l <= r && r <= _n);
        S sml = E(), smr = E();
        l += size;
        r += size;
        while (l < r) {
            if (l & 1) sml = Op(sml, d[l++]);
            if (r & 1) smr = Op(d[--r], smr);
            l >>= 1;
            r >>= 1;
        }
        return Op(sml, smr);
    }

    S all_prod() const { return d[1]; }

    template <class F>
        requires std::predicate<F, S>
    int max_right(int l, F f) const {
        assert(0 <= l && l <= _n);
        assert(f(E()));
        if (l == _n) return _n;
        l += size;
        S sm = E();
        do {
            while ((l & 1) == 0)
                l >>= 1;
            if (!f(Op(sm, d[l]))) {
                while (l < size) {
                    l = (2 * l);
                    if (f(Op(sm, d[l]))) {
                        sm = Op(sm, d[l]);
                        l++;
                    }
                }
                return l - size;
            }
            sm = Op(sm, d[l]);
            l++;
        } while ((l & -l) != l);
        return _n;
    }

    template <class F>
        requires std::predicate<F, S>
    int min_left(int r, F f) const {
        assert(0 <= r && r <= _n);
        assert(f(E()));
        if (r == 0) return 0;
        r += size;
        S sm = E();
        do {
            r--;
            while (r > 1 && (r & 1))
                r >>= 1;
            if (!f(Op(d[r], sm))) {
                while (r < size) {
                    r = (2 * r + 1);
                    if (f(Op(d[r], sm))) {
                        sm = Op(d[r], sm);
                        r--;
                    }
                }
                return r + 1 - size;
            }
            sm = Op(d[r], sm);
        } while ((r & -r) != r);
        return 0;
    }
};

// --- 頻出演算の定義用ヘルパー ---
namespace seg_ops {
template <typename T>
struct Min {
    static constexpr T op(T a, T b) { return std::min(a, b); }
    static constexpr T e() { return std::numeric_limits<T>::max(); }
};

template <typename T>
struct Max {
    static constexpr T op(T a, T b) { return std::max(a, b); }
    static constexpr T e() { return std::numeric_limits<T>::lowest(); }
};

template <typename T>
struct Add {
    static constexpr T op(T a, T b) { return a + b; }
    static constexpr T e() { return T{0}; }
};

template <typename T>
struct Gcd {
    static constexpr T op(T a, T b) { return std::gcd(a, b); }
    static constexpr T e() { return T{0}; }
};

template <typename T>
struct BitXor {
    static constexpr T op(T a, T b) { return a ^ b; }
    static constexpr T e() { return T{0}; }
};

template <typename T>
struct BitAnd {
    static constexpr T op(T a, T b) { return a & b; }
    static constexpr T e() { return ~T{0}; }
};

template <typename T>
struct BitOr {
    static constexpr T op(T a, T b) { return a | b; }
    static constexpr T e() { return T{0}; }
};
}  // namespace seg_ops

// --- 頻出セグ木の型エイリアス (デフォルトは default_type = long) ---

using default_type = long;
template <typename T = default_type>
using RangeMin = SegTree<T, seg_ops::Min<T>::op, seg_ops::Min<T>::e>;

template <typename T = default_type>
using RangeMax = SegTree<T, seg_ops::Max<T>::op, seg_ops::Max<T>::e>;

template <typename T = default_type>
using RangeSum = SegTree<T, seg_ops::Add<T>::op, seg_ops::Add<T>::e>;

template <typename T = default_type>
using RangeGcd = SegTree<T, seg_ops::Gcd<T>::op, seg_ops::Gcd<T>::e>;

template <typename T = default_type>
using RangeXor = SegTree<T, seg_ops::BitXor<T>::op, seg_ops::BitXor<T>::e>;

template <typename T = default_type>
using RangeAnd = SegTree<T, seg_ops::BitAnd<T>::op, seg_ops::BitAnd<T>::e>;

template <typename T = default_type>
using RangeOr = SegTree<T, seg_ops::BitOr<T>::op, seg_ops::BitOr<T>::e>;