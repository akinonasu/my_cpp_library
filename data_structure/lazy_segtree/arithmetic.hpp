#pragma once

#include <cassert>
#include <vector>

#include "../lazy_segtree.hpp"

template <typename T = long>
struct ArithNode {
    T sum_a = 0;  // sum(A[i])
    T sum_i = 0;  // sum(i)
    T len = 0;    // 区間長

    ArithNode() = default;
    ArithNode(T idx, T val) : sum_a(val), sum_i(idx), len(1) {}
    ArithNode(T sa, T si, T l) : sum_a(sa), sum_i(si), len(l) {}
};

// A[i] に c0 + c1 * i を加算する作用素
template <typename T = long>
struct ArithAdd {
    T c0 = 0;
    T c1 = 0;
};

namespace arith_ops {
template <typename T>
constexpr ArithNode<T> op(const ArithNode<T>& l, const ArithNode<T>& r) {
    return {l.sum_a + r.sum_a, l.sum_i + r.sum_i, l.len + r.len};
}

template <typename T>
constexpr ArithNode<T> e() {
    return {0, 0, 0};
}

template <typename T>
constexpr ArithNode<T> mapping(const ArithAdd<T>& f, const ArithNode<T>& s) {
    return {s.sum_a + f.c0 * s.len + f.c1 * s.sum_i, s.sum_i, s.len};
}

template <typename T>
constexpr ArithAdd<T> composition(const ArithAdd<T>& f, const ArithAdd<T>& g) {
    return {f.c0 + g.c0, f.c1 + g.c1};
}

template <typename T>
constexpr ArithAdd<T> id() {
    return {0, 0};
}
}  // namespace arith_ops

// --- ラッパークラス ---
template <typename T = long>
class RangeArithmetic {
   private:
    using Node = ArithNode<T>;
    using F = ArithAdd<T>;
    using Tree =
        LazySegTree<Node, arith_ops::op<T>, arith_ops::e<T>, F, arith_ops::mapping<T>, arith_ops::composition<T>, arith_ops::id<T>>;
    Tree seg;
    int _n;

   public:
    RangeArithmetic(int n, T val = 0) : RangeArithmetic(std::vector<T>(n, val)) {}
    RangeArithmetic(const std::vector<T>& a) : _n(static_cast<int>(a.size())) {
        std::vector<Node> init(_n);
        for (int i = 0; i < _n; i++)
            init[i] = Node(static_cast<T>(i), a[i]);
        seg = Tree(init);
    }

    // 通常の区間加算 (+val)
    void apply_add(int l, int r, T val) { seg.apply(l, r, F{.c0 = val, .c1 = 0}); }

    // 【メイン機能】等差数列加算:
    // A[l] += a, A[l+1] += a+d, A[l+2] += a+2d, ... を加算
    void apply_arithmetic_add(int l, int r, T a, T d) {
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return;
        // a + (i - l)*d = (a - l*d) + d * i
        T c0 = a - static_cast<T>(l) * d;
        T c1 = d;
        seg.apply(l, r, F{.c0 = c0, .c1 = c1});
    }

    // 区間和 sum(A[i])
    T prod_sum(int l, int r) { return seg.prod(l, r).sum_a; }
};
