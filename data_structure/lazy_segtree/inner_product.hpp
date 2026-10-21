#pragma once

#include <cassert>
#include <vector>

#include "../lazy_segtree.hpp"

// ノード: sum(A*B), sum(A), sum(B), len
template <typename T = long>
struct InnerProductNode {
    T sum_ab = 0;
    T sum_a = 0;
    T sum_b = 0;
    T len = 0;

    InnerProductNode() = default;
    InnerProductNode(T a, T b) : sum_ab(a * b), sum_a(a), sum_b(b), len(1) {}
    InnerProductNode(T ab, T a, T b, T l) : sum_ab(ab), sum_a(a), sum_b(b), len(l) {}
};

// 作用素: (Aへの加算値, Bへの加算値)
template <typename T = long>
struct AddPair {
    T add_a = 0;
    T add_b = 0;
};

namespace inner_product_ops {
template <typename T>
constexpr InnerProductNode<T> op(const InnerProductNode<T>& l, const InnerProductNode<T>& r) {
    return {l.sum_ab + r.sum_ab, l.sum_a + r.sum_a, l.sum_b + r.sum_b, l.len + r.len};
}

template <typename T>
constexpr InnerProductNode<T> e() {
    return {0, 0, 0, 0};
}

template <typename T>
constexpr InnerProductNode<T> mapping(const AddPair<T>& f, const InnerProductNode<T>& s) {
    return {s.sum_ab + f.add_a * s.sum_b + f.add_b * s.sum_a + f.add_a * f.add_b * s.len, s.sum_a + f.add_a * s.len,
            s.sum_b + f.add_b * s.len, s.len};
}

template <typename T>
constexpr AddPair<T> composition(const AddPair<T>& f, const AddPair<T>& g) {
    return {f.add_a + g.add_a, f.add_b + g.add_b};
}

template <typename T>
constexpr AddPair<T> id() {
    return {0, 0};
}
}  // namespace inner_product_ops

// --- ラッパークラス ---
template <typename T = long>
class RangeInnerProduct {
   private:
    using Node = InnerProductNode<T>;
    using F = AddPair<T>;
    using Tree = LazySegTree<Node, inner_product_ops::op<T>, inner_product_ops::e<T>, F, inner_product_ops::mapping<T>,
                             inner_product_ops::composition<T>, inner_product_ops::id<T>>;
    Tree seg;
    int _n;

   public:
    RangeInnerProduct(const std::vector<T>& a, const std::vector<T>& b) : _n(static_cast<int>(a.size())) {
        assert(a.size() == b.size());
        std::vector<Node> init(_n);
        for (int i = 0; i < _n; i++)
            init[i] = Node(a[i], b[i]);
        seg = Tree(init);
    }

    // A[l, r) に val を加算
    void apply_add_a(int l, int r, T val) { seg.apply(l, r, F{.add_a = val, .add_b = 0}); }

    // B[l, r) に val を加算
    void apply_add_b(int l, int r, T val) { seg.apply(l, r, F{.add_a = 0, .add_b = val}); }

    // 両方に同時に加算
    void apply_add_both(int l, int r, T val_a, T val_b) { seg.apply(l, r, F{.add_a = val_a, .add_b = val_b}); }

    // sum(A[i] * B[i]) を取得
    T prod_inner_product(int l, int r) { return seg.prod(l, r).sum_ab; }

    // sum(A[i]) を取得
    T prod_sum_a(int l, int r) { return seg.prod(l, r).sum_a; }

    // sum(B[i]) を取得
    T prod_sum_b(int l, int r) { return seg.prod(l, r).sum_b; }
};
