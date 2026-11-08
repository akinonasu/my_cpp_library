#pragma once

#include <algorithm>
#include <limits>
#include <vector>

#include "data_structure/segtree.hpp"

namespace applied_seg {

// =========================================================================
// 1. ArgMin / ArgMax (値とそのインデックスを同時に管理)
// =========================================================================
template <typename T = long>
struct ArgMin {
    struct S {
        T val;
        int idx;
    };
    static constexpr S op(S a, S b) {
        if (a.val <= b.val) return a;
        return b;
    }
    static constexpr S e() { return S{std::numeric_limits<T>::max(), -1}; }
};

template <typename T = long>
struct ArgMax {
    struct S {
        T val;
        int idx;
    };
    static constexpr S op(S a, S b) {
        if (a.val >= b.val) return a;
        return b;
    }
    static constexpr S e() { return S{std::numeric_limits<T>::lowest(), -1}; }
};

// =========================================================================
// 2. 最大部分配列和 (Maximum Subarray Sum)
// 区間内の空でない連続する部分列の総和の最大値を O(log N) で取得
// =========================================================================
template <typename T = long>
struct MaxSubarraySum {
    static constexpr T INF = std::numeric_limits<T>::max() / 4;

    struct S {
        T sum;      // 区間全体の総和
        T lmax;     // 左端から始まる最大プレフィックス和
        T rmax;     // 右端で終わる最大サフィックス和
        T max_sum;  // 区間内の最大部分配列和

        // 単位元判定
        constexpr bool is_empty() const { return max_sum <= -INF / 2; }
    };

    static constexpr S make(T val) { return S{val, val, val, val}; }

    static constexpr S op(S l, S r) {
        if (l.is_empty()) return r;
        if (r.is_empty()) return l;

        return S{l.sum + r.sum, std::max(l.lmax, l.sum + r.lmax), std::max(r.rmax, r.sum + l.rmax), std::max({l.max_sum, r.max_sum, l.rmax + r.lmax})};
    }

    static constexpr S e() { return S{0, -INF, -INF, -INF}; }
};

// =========================================================================
// 3. アフィン変換 / 一次関数の合成 (Composite Functions)
// f(x) = a * x + b の合成: (g ∘ f)(x) = g(f(x)) = c*(a*x + b) + d
// =========================================================================
template <typename T = long, T Mod = 998244353>
struct Affine {
    struct S {
        T a, b;  // f(x) = a * x + b
        constexpr T eval(T x) const { return (a * x + b) % Mod; }
    };

    // f の後に g を適用: g(f(x))
    static constexpr S op(S f, S g) { return S{(g.a * f.a) % Mod, (g.a * f.b + g.b) % Mod}; }

    // 恒等写像: id(x) = 1 * x + 0
    static constexpr S e() { return S{1, 0}; }
};

// =========================================================================
// 4. 括弧列クエリ (Bracket Sequence)
// 正しい括弧列のペア数、マッチしていない '(' と ')' の数を管理
// =========================================================================
struct BracketSequence {
    struct S {
        int open;     // 未マッチの '('
        int close;    // 未マッチの ')'
        int matched;  // マッチしたペア数
    };

    static constexpr S make(char c) {
        if (c == '(') return S{1, 0, 0};
        if (c == ')') return S{0, 1, 0};
        return S{0, 0, 0};
    }

    static constexpr S op(S l, S r) {
        int match = std::min(l.open, r.close);
        return S{l.open + r.open - match, l.close + r.close - match, l.matched + r.matched + match};
    }

    static constexpr S e() { return S{0, 0, 0}; }
};

}  // namespace applied_seg

// --- エイリアス定義 ---
template <typename T = long>
using RangeArgMin = SegTree<typename applied_seg::ArgMin<T>::S, applied_seg::ArgMin<T>::op, applied_seg::ArgMin<T>::e>;

template <typename T = long>
using RangeArgMax = SegTree<typename applied_seg::ArgMax<T>::S, applied_seg::ArgMax<T>::op, applied_seg::ArgMax<T>::e>;

template <typename T = long>
using RangeMaxSubarraySum = SegTree<typename applied_seg::MaxSubarraySum<T>::S, applied_seg::MaxSubarraySum<T>::op, applied_seg::MaxSubarraySum<T>::e>;

template <typename T = long, T Mod = 998244353>
using RangeAffine = SegTree<typename applied_seg::Affine<T, Mod>::S, applied_seg::Affine<T, Mod>::op, applied_seg::Affine<T, Mod>::e>;

using RangeBracket = SegTree<applied_seg::BracketSequence::S, applied_seg::BracketSequence::op, applied_seg::BracketSequence::e>;