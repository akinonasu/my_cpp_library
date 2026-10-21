#pragma once

#include <cassert>
#include <format>
#include <vector>

#include "../lazy_segtree.hpp"

// --- ノード定義 ---
template <typename T = long>
struct IndexSumNode {
    T sum_a = 0;   // sum(A[i])
    T sum_ia = 0;  // sum(i * A[i])
    T sum_i = 0;   // sum(i)
    T len = 0;     // 区間幅

    IndexSumNode() = default;
    IndexSumNode(T idx, T val) : sum_a(val), sum_ia(idx * val), sum_i(idx), len(1) {}
    IndexSumNode(T sa, T sia, T si, T l) : sum_a(sa), sum_ia(sia), sum_i(si), len(l) {}
};

// std::print 用フォーマッタ (sum_ia と sum_a を表示)
template <typename T>
struct std::formatter<IndexSumNode<T>> {
    constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }
    auto format(const IndexSumNode<T>& node, auto& ctx) const {
        return std::format_to(ctx.out(), "{{sum_ia: {}, sum_a: {}}}", node.sum_ia, node.sum_a);
    }
};

// 作用素 (一次関数 f(x) = ax + b)
template <typename T = long>
struct IndexSumAffine {
    T a = 1;
    T b = 0;
};

// 演算の定義
namespace index_sum_ops {
template <typename T>
constexpr IndexSumNode<T> op(const IndexSumNode<T>& l, const IndexSumNode<T>& r) {
    return IndexSumNode<T>{l.sum_a + r.sum_a, l.sum_ia + r.sum_ia, l.sum_i + r.sum_i, l.len + r.len};
}

template <typename T>
constexpr IndexSumNode<T> e() {
    return IndexSumNode<T>{0, 0, 0, 0};
}

template <typename T>
constexpr IndexSumNode<T> mapping(const IndexSumAffine<T>& f, const IndexSumNode<T>& s) {
    return IndexSumNode<T>{s.sum_a * f.a + f.b * s.len, s.sum_ia * f.a + f.b * s.sum_i, s.sum_i, s.len};
}

template <typename T>
constexpr IndexSumAffine<T> composition(const IndexSumAffine<T>& f, const IndexSumAffine<T>& g) {
    return IndexSumAffine<T>{f.a * g.a, f.a * g.b + f.b};
}

template <typename T>
constexpr IndexSumAffine<T> id() {
    return IndexSumAffine<T>{1, 0};
}
}  // namespace index_sum_ops

// --- ラッパークラス (型引数デフォルト: long) ---
template <typename T = long>
class RangeIndexSum {
   public:
    using Node = IndexSumNode<T>;
    using F = IndexSumAffine<T>;

   private:
    using Tree = LazySegTree<Node, index_sum_ops::op<T>, index_sum_ops::e<T>, F, index_sum_ops::mapping<T>, index_sum_ops::composition<T>,
                             index_sum_ops::id<T>>;
    Tree seg;
    int _n;

   public:
    RangeIndexSum() : RangeIndexSum(0) {}
    explicit RangeIndexSum(int n, T val = 0) : RangeIndexSum(std::vector<T>(n, val)) {}
    explicit RangeIndexSum(const std::vector<T>& a) : _n(static_cast<int>(a.size())) {
        std::vector<Node> init(_n);
        for (int i = 0; i < _n; i++) {
            init[i] = Node(static_cast<T>(i), a[i]);
        }
        seg = Tree(init);
    }

    // [l, r) に +val を加算
    void apply_add(int l, int r, T val) { seg.apply(l, r, F{.a = 1, .b = val}); }

    // [l, r) を val に一括代入
    void apply_update(int l, int r, T val) { seg.apply(l, r, F{.a = 0, .b = val}); }

    // [l, r) に一次関数 x -> a*x + b を適用
    void apply_affine(int l, int r, T a, T b) { seg.apply(l, r, F{.a = a, .b = b}); }

    // a[p] = val に 1点更新
    void set(int p, T val) {
        assert(0 <= p && p < _n);
        seg.set(p, Node(static_cast<T>(p), val));
    }

    // a[p] の値を取得
    T get(int p) {
        assert(0 <= p && p < _n);
        return seg.get(p).sum_a;
    }

    // 区間のノード情報そのまま [l, r)
    Node prod(int l, int r) { return seg.prod(l, r); }

    // 区間和: sum_{i=l}^{r-1} A[i]
    T prod_sum(int l, int r) { return seg.prod(l, r).sum_a; }

    // 重み付き区間和: sum_{i=l}^{r-1} (i * A[i])
    T prod_weighted_sum(int l, int r) { return seg.prod(l, r).sum_ia; }

    // 【典型クエリ】左端を 1 とした相対重み付き和:
    // sum_{k=0}^{r-l-1} (k + 1) * A[l + k]
    // = 1*A[l] + 2*A[l+1] + 3*A[l+2] + ...
    T prod_relative_weighted_sum(int l, int r) {
        if (l >= r) return 0;
        Node res = seg.prod(l, r);
        return res.sum_ia - static_cast<T>(l - 1) * res.sum_a;
    }
};
