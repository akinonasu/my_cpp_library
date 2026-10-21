#include <algorithm>
#include <bit>
#include <cassert>
#include <concepts>
#include <iostream>
#include <limits>
#include <optional>
#include <print>
#include <vector>

/*
// 区間加算・区間最小値 RangeAddMin
// 区間加算・区間最大値 RangeAddMax
// 区間変更・区間最小値 RangeUpdateMin
// 区間変更・区間最大値 RangeUpdateMax
// 区間変更・区間和 RangeUpdateSum
// 区間アフィン変換 RangeAffineSum
*/

// --- 遅延セグメントツリー本体 ---
// S: 要素の型
// Op: 要素同士のマージ (S, S) -> S
// E: 要素の単位元 () -> S
// F: 作用素の型
// Mapping: 作用素を要素に適用 (F, S) -> S
// Composition: 作用素同士の合成 (F 新, F 旧) -> F
// Id: 作用素の恒等写像 () -> F
template <class S, auto Op, auto E, class F, auto Mapping, auto Composition, auto Id>
    requires std::invocable<decltype(Op), S, S> && std::same_as<std::invoke_result_t<decltype(Op), S, S>, S> &&
             std::invocable<decltype(E)> && std::same_as<std::invoke_result_t<decltype(E)>, S> && std::invocable<decltype(Mapping), F, S> &&
             std::same_as<std::invoke_result_t<decltype(Mapping), F, S>, S> && std::invocable<decltype(Composition), F, F> &&
             std::same_as<std::invoke_result_t<decltype(Composition), F, F>, F> && std::invocable<decltype(Id)> &&
             std::same_as<std::invoke_result_t<decltype(Id)>, F>
class LazySegTree {
   private:
    int _n;
    int size;
    int log;
    std::vector<S> d;
    std::vector<F> lz;

    void update(int k) { d[k] = Op(d[2 * k], d[2 * k + 1]); }

    void all_apply(int k, F f) {
        d[k] = Mapping(f, d[k]);
        if (k < size) lz[k] = Composition(f, lz[k]);
    }

    void push(int k) {
        all_apply(2 * k, lz[k]);
        all_apply(2 * k + 1, lz[k]);
        lz[k] = Id();
    }

   public:
    LazySegTree() : LazySegTree(0) {}
    explicit LazySegTree(int n) : LazySegTree(std::vector<S>(n, E())) {}

    explicit LazySegTree(const std::vector<S>& v) : _n(static_cast<int>(v.size())) {
        size = std::bit_ceil(static_cast<unsigned int>(_n));
        log = std::countr_zero(static_cast<unsigned int>(size));
        d.assign(2 * size, E());
        lz.assign(size, Id());
        for (int i = 0; i < _n; i++)
            d[size + i] = v[i];
        for (int i = size - 1; i >= 1; i--)
            update(i);
    }

    // 別型 (例: long) の配列から S への暗黙/明示変換コンストラクタ
    template <typename U>
        requires(!std::same_as<U, S> && std::constructible_from<S, U>)
    explicit LazySegTree(const std::vector<U>& v) : _n(static_cast<int>(v.size())) {
        size = std::bit_ceil(static_cast<unsigned int>(_n));
        log = std::countr_zero(static_cast<unsigned int>(size));
        d.assign(2 * size, E());
        lz.assign(size, Id());
        for (int i = 0; i < _n; i++)
            d[size + i] = S(v[i]);
        for (int i = size - 1; i >= 1; i--)
            update(i);
    }

    // 1点更新
    void set(int p, S x) {
        assert(0 <= p && p < _n);
        p += size;
        for (int i = log; i >= 1; i--)
            push(p >> i);
        d[p] = x;
        for (int i = 1; i <= log; i++)
            update(p >> i);
    }

    // 1点取得
    S get(int p) {
        assert(0 <= p && p < _n);
        p += size;
        for (int i = log; i >= 1; i--)
            push(p >> i);
        return d[p];
    }

    // 区間取得 [l, r)
    S prod(int l, int r) {
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return E();
        l += size;
        r += size;
        for (int i = log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
        S sml = E(), smr = E();
        while (l < r) {
            if (l & 1) sml = Op(sml, d[l++]);
            if (r & 1) smr = Op(d[--r], smr);
            l >>= 1;
            r >>= 1;
        }
        return Op(sml, smr);
    }

    S all_prod() const { return d[1]; }

    // 1点作用
    void apply(int p, F f) {
        assert(0 <= p && p < _n);
        p += size;
        for (int i = log; i >= 1; i--)
            push(p >> i);
        d[p] = Mapping(f, d[p]);
        for (int i = 1; i <= log; i++)
            update(p >> i);
    }

    // 区間作用 [l, r) に f を適用
    void apply(int l, int r, F f) {
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return;
        l += size;
        r += size;
        for (int i = log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
        {
            int l2 = l, r2 = r;
            while (l < r) {
                if (l & 1) all_apply(l++, f);
                if (r & 1) all_apply(--r, f);
                l >>= 1;
                r >>= 1;
            }
            l = l2;
            r = r2;
        }
        for (int i = 1; i <= log; i++) {
            if (((l >> i) << i) != l) update(l >> i);
            if (((r >> i) << i) != r) update((r - 1) >> i);
        }
    }

    // セグ木上の二分探索 (右方向)
    template <class Pred>
        requires std::predicate<Pred, S>
    int max_right(int l, Pred g) {
        assert(0 <= l && l <= _n);
        assert(g(E()));
        if (l == _n) return _n;
        l += size;
        for (int i = log; i >= 1; i--)
            push(l >> i);
        S sm = E();
        do {
            while (l % 2 == 0)
                l >>= 1;
            if (!g(Op(sm, d[l]))) {
                while (l < size) {
                    push(l);
                    l = (2 * l);
                    if (g(Op(sm, d[l]))) {
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

    // セグ木上の二分探索 (左方向)
    template <class Pred>
        requires std::predicate<Pred, S>
    int min_left(int r, Pred g) {
        assert(0 <= r && r <= _n);
        assert(g(E()));
        if (r == 0) return 0;
        r += size;
        for (int i = log; i >= 1; i--)
            push((r - 1) >> i);
        S sm = E();
        do {
            r--;
            while (r > 1 && (r % 2))
                r >>= 1;
            if (!g(Op(d[r], sm))) {
                while (r < size) {
                    push(r);
                    r = (2 * r + 1);
                    if (g(Op(d[r], sm))) {
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

// --- 区間和クエリ用のノード定義 ---
template <typename T>
struct SumNode {
    T sum = 0;
    T size = 1;
    SumNode() : sum(0), size(0) {}
    SumNode(T val) : sum(val), size(1) {}
    SumNode(T sum, T size) : sum(sum), size(size) {}
    // 自動で T (long など) にキャストできるようにする
    operator T() const { return sum; }
};

// std::print 用のフォーマッタ (SumNode をそのまま print 可能にする)
template <typename T>
struct std::formatter<SumNode<T>> : std::formatter<T> {
    auto format(const SumNode<T>& node, auto& ctx) const { return std::formatter<T>::format(node.sum, ctx); }
};

// アフィン変換用構造体 f(x) = a * x + b
template <typename T>
struct Affine {
    T a = 1;
    T b = 0;
};

// --- 頻出演算ヘルパー ---
namespace lazy_ops {
// 1. Range Add & Range Min
template <typename T>
struct AddMin {
    static constexpr T INF = std::numeric_limits<T>::max() / 2;
    static constexpr T op(T a, T b) { return std::min(a, b); }
    static constexpr T e() { return INF; }
    static constexpr T mapping(T f, T s) { return (s >= INF ? INF : s + f); }
    static constexpr T composition(T f, T g) { return f + g; }
    static constexpr T id() { return 0; }
};

// 2. Range Add & Range Max
template <typename T>
struct AddMax {
    static constexpr T INF = std::numeric_limits<T>::max() / 2;
    static constexpr T op(T a, T b) { return std::max(a, b); }
    static constexpr T e() { return -INF; }
    static constexpr T mapping(T f, T s) { return (s <= -INF ? -INF : s + f); }
    static constexpr T composition(T f, T g) { return f + g; }
    static constexpr T id() { return 0; }
};

// 3. Range Add & Range Sum
template <typename T>
struct AddSum {
    using S = SumNode<T>;
    static constexpr S op(S a, S b) { return {a.sum + b.sum, a.size + b.size}; }
    static constexpr S e() { return {0, 0}; }
    static constexpr S mapping(T f, S s) { return {s.sum + f * s.size, s.size}; }
    static constexpr T composition(T f, T g) { return f + g; }
    static constexpr T id() { return 0; }
};

// 4. Range Update & Range Min
template <typename T>
struct UpdateMin {
    static constexpr T INF = std::numeric_limits<T>::max() / 2;
    using F = std::optional<T>;
    static constexpr T op(T a, T b) { return std::min(a, b); }
    static constexpr T e() { return INF; }
    static constexpr T mapping(F f, T s) { return f.value_or(s); }
    static constexpr F composition(F f, F g) { return f.has_value() ? f : g; }
    static constexpr F id() { return std::nullopt; }
};

// 5. Range Update & Range Max
template <typename T>
struct UpdateMax {
    static constexpr T INF = std::numeric_limits<T>::max() / 2;
    using F = std::optional<T>;
    static constexpr T op(T a, T b) { return std::max(a, b); }
    static constexpr T e() { return -INF; }
    static constexpr T mapping(F f, T s) { return f.value_or(s); }
    static constexpr F composition(F f, F g) { return f.has_value() ? f : g; }
    static constexpr F id() { return std::nullopt; }
};

// 6. Range Update & Range Sum
template <typename T>
struct UpdateSum {
    using S = SumNode<T>;
    using F = std::optional<T>;
    static constexpr S op(S a, S b) { return {a.sum + b.sum, a.size + b.size}; }
    static constexpr S e() { return {0, 0}; }
    static constexpr S mapping(F f, S s) { return f.has_value() ? S{*f * s.size, s.size} : s; }
    static constexpr F composition(F f, F g) { return f.has_value() ? f : g; }
    static constexpr F id() { return std::nullopt; }
};

// 7. Range Affine & Range Sum (x -> a*x + b)
template <typename T>
struct AffineSum {
    using S = SumNode<T>;
    using F = Affine<T>;
    static constexpr S op(S a, S b) { return {a.sum + b.sum, a.size + b.size}; }
    static constexpr S e() { return {0, 0}; }
    static constexpr S mapping(F f, S s) { return {s.sum * f.a + f.b * s.size, s.size}; }
    static constexpr F composition(F f, F g) { return {f.a * g.a, f.a * g.b + f.b}; }
    static constexpr F id() { return {1, 0}; }
};
}  // namespace lazy_ops

// --- 頻出遅延セグ木の型エイリアス (デフォルト: long) ---

// 区間加算・区間最小値
template <typename T = long>
using RangeAddMin = LazySegTree<T, lazy_ops::AddMin<T>::op, lazy_ops::AddMin<T>::e, T, lazy_ops::AddMin<T>::mapping,
                                lazy_ops::AddMin<T>::composition, lazy_ops::AddMin<T>::id>;

// 区間加算・区間最大値
template <typename T = long>
using RangeAddMax = LazySegTree<T, lazy_ops::AddMax<T>::op, lazy_ops::AddMax<T>::e, T, lazy_ops::AddMax<T>::mapping,
                                lazy_ops::AddMax<T>::composition, lazy_ops::AddMax<T>::id>;

// 区間加算・区間和
template <typename T = long>
using RangeAddSum = LazySegTree<SumNode<T>, lazy_ops::AddSum<T>::op, lazy_ops::AddSum<T>::e, T, lazy_ops::AddSum<T>::mapping,
                                lazy_ops::AddSum<T>::composition, lazy_ops::AddSum<T>::id>;

// 区間変更・区間最小値
template <typename T = long>
using RangeUpdateMin = LazySegTree<T, lazy_ops::UpdateMin<T>::op, lazy_ops::UpdateMin<T>::e, std::optional<T>,
                                   lazy_ops::UpdateMin<T>::mapping, lazy_ops::UpdateMin<T>::composition, lazy_ops::UpdateMin<T>::id>;

// 区間変更・区間最大値
template <typename T = long>
using RangeUpdateMax = LazySegTree<T, lazy_ops::UpdateMax<T>::op, lazy_ops::UpdateMax<T>::e, std::optional<T>,
                                   lazy_ops::UpdateMax<T>::mapping, lazy_ops::UpdateMax<T>::composition, lazy_ops::UpdateMax<T>::id>;

// 区間変更・区間和
template <typename T = long>
using RangeUpdateSum = LazySegTree<SumNode<T>, lazy_ops::UpdateSum<T>::op, lazy_ops::UpdateSum<T>::e, std::optional<T>,
                                   lazy_ops::UpdateSum<T>::mapping, lazy_ops::UpdateSum<T>::composition, lazy_ops::UpdateSum<T>::id>;

// 区間アフィン変換 (a*x + b)・区間和
template <typename T = long>
using RangeAffineSum = LazySegTree<SumNode<T>, lazy_ops::AffineSum<T>::op, lazy_ops::AffineSum<T>::e, Affine<T>,
                                   lazy_ops::AffineSum<T>::mapping, lazy_ops::AffineSum<T>::composition, lazy_ops::AffineSum<T>::id>;