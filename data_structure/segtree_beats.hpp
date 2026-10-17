#pragma once

#include <algorithm>
#include <cassert>
#include <limits>
#include <type_traits>
#include <vector>

#include "components/default_init_allocator.hpp"

template <class T = long long, class Sum = T>
class SegtreeBeats {
    static_assert(std::is_integral_v<T> && std::is_signed_v<T>);
    static constexpr T inf = std::numeric_limits<T>::max();
    static constexpr T neg_inf = std::numeric_limits<T>::lowest();
    struct Node {
        Sum sum;
        T lo, hi, lo2, hi2, add;
        int nlo, nhi;
    };
    int n = 0;
    mutable std::vector<Node, DefaultInitAllocator<Node>> tree;
    static int left(int k) { return k + 1; }
    static int right(int k, int l, int m) { return k + 2 * (m - l); }
    static void fill(Node& a, T x, int len) { a = {Sum(x) * len, x, x, inf, neg_inf, 0, len, len}; }
    static void add_node(Node& a, T x, int len) {
        a.sum += Sum(x) * len;
        a.lo += x;
        a.hi += x;
        if (a.lo2 != inf) a.lo2 += x;
        if (a.hi2 != neg_inf) a.hi2 += x;
        a.add += x;
    }
    static void lower(Node& a, T x) {
        a.sum += (Sum(x) - Sum(a.hi)) * a.nhi;
        if (a.lo == a.hi)
            a.lo = x;
        else if (a.lo2 == a.hi)
            a.lo2 = x;
        a.hi = x;
    }
    static void raise(Node& a, T x) {
        a.sum += (Sum(x) - Sum(a.lo)) * a.nlo;
        if (a.hi == a.lo)
            a.hi = x;
        else if (a.hi2 == a.lo)
            a.hi2 = x;
        a.lo = x;
    }
    void pull(int k, int lc, int rc) {
        auto& a = tree[k];
        const auto &b = tree[lc], &c = tree[rc];
        a.sum = b.sum + c.sum;
        if (b.lo < c.lo)
            a.lo = b.lo, a.lo2 = std::min(b.lo2, c.lo), a.nlo = b.nlo;
        else if (b.lo > c.lo)
            a.lo = c.lo, a.lo2 = std::min(b.lo, c.lo2), a.nlo = c.nlo;
        else
            a.lo = b.lo, a.lo2 = std::min(b.lo2, c.lo2), a.nlo = b.nlo + c.nlo;
        if (b.hi > c.hi)
            a.hi = b.hi, a.hi2 = std::max(b.hi2, c.hi), a.nhi = b.nhi;
        else if (b.hi < c.hi)
            a.hi = c.hi, a.hi2 = std::max(b.hi, c.hi2), a.nhi = c.nhi;
        else
            a.hi = b.hi, a.hi2 = std::max(b.hi2, c.hi2), a.nhi = b.nhi + c.nhi;
    }
    void push(int k, int l, int m, int r) const {
        auto &a = tree[k], &b = tree[left(k)], &c = tree[right(k, l, m)];
        if (a.lo == a.hi) {
            fill(b, a.lo, m - l);
            fill(c, a.lo, r - m);
            a.add = 0;
            return;
        }
        if (a.add) {
            add_node(b, a.add, m - l);
            add_node(c, a.add, r - m);
            a.add = 0;
        }
        if (b.hi > a.hi) lower(b, a.hi);
        if (b.lo < a.lo) raise(b, a.lo);
        if (c.hi > a.hi) lower(c, a.hi);
        if (c.lo < a.lo) raise(c, a.lo);
    }
    template <class F>
    void build(int k, int l, int r, const F& f) {
        if (r - l == 1) return fill(tree[k], f(l), 1);
        int m = (l + r) / 2, lc = left(k), rc = right(k, l, m);
        build(lc, l, m, f);
        build(rc, m, r, f);
        tree[k].add = 0;
        pull(k, lc, rc);
    }
    template <int op>
    void update(int k, int l, int r, int ql, int qr, T x) {
        auto& a = tree[k];
        if constexpr (op == 0) {
            if (a.hi <= x) return;
        }
        if constexpr (op == 1) {
            if (a.lo >= x) return;
        }
        if constexpr (op == 3) {
            if (a.lo == x && a.hi == x) return;
        }
        if (ql <= l && r <= qr) {
            if constexpr (op == 0) {
                if (x <= a.lo) return fill(a, x, r - l);
                if (a.hi2 < x) return lower(a, x);
            }
            if constexpr (op == 1) {
                if (x >= a.hi) return fill(a, x, r - l);
                if (a.lo2 > x) return raise(a, x);
            }
            if constexpr (op == 2) return add_node(a, x, r - l);
            if constexpr (op == 3) return fill(a, x, r - l);
        }
        int m = (l + r) / 2, lc = left(k), rc = right(k, l, m);
        push(k, l, m, r);
        if (ql < m) update<op>(lc, l, m, ql, qr, x);
        if (m < qr) update<op>(rc, m, r, ql, qr, x);
        pull(k, lc, rc);
    }
    template <int op>
    using Result = std::conditional_t<op == 0, Sum, T>;
    template <int op>
    Result<op> query(int k, int l, int r, int ql, int qr) const {
        const auto& a = tree[k];
        if (ql <= l && r <= qr) {
            if constexpr (op == 0) return a.sum;
            if constexpr (op == 1) return a.lo;
            if constexpr (op == 2) return a.hi;
        }
        if (a.lo == a.hi) {
            if constexpr (op == 0)
                return Sum(a.lo) * (std::min(r, qr) - std::max(l, ql));
            else
                return a.lo;
        }
        int m = (l + r) / 2;
        push(k, l, m, r);
        if (qr <= m) return query<op>(left(k), l, m, ql, qr);
        if (m <= ql) return query<op>(right(k, l, m), m, r, ql, qr);
        auto b = query<op>(left(k), l, m, ql, qr), c = query<op>(right(k, l, m), m, r, ql, qr);
        if constexpr (op == 0) return b + c;
        if constexpr (op == 1) return std::min(b, c);
        if constexpr (op == 2) return std::max(b, c);
    }
    void collect(int k, int l, int r, std::vector<T>& a) const {
        if (tree[k].lo == tree[k].hi) {
            std::fill(a.begin() + l, a.begin() + r, tree[k].lo);
            return;
        }
        int m = (l + r) / 2;
        push(k, l, m, r);
        collect(left(k), l, m, a);
        collect(right(k, l, m), m, r, a);
    }
    void check(int l, int r) const {
        assert(0 <= l && l <= r && r <= n);
        (void)l;
        (void)r;
    }

   public:
    SegtreeBeats() = default;
    explicit SegtreeBeats(int n, T x = 0) : n(n), tree(n > 0 ? 2 * std::size_t(n) - 1 : 0) {
        assert(n >= 0);
        if (n) build(0, 0, n, [x](int) { return x; });
    }
    explicit SegtreeBeats(const std::vector<T>& a) : n(int(a.size())), tree(n ? 2 * std::size_t(n) - 1 : 0) {
        if (n) build(0, 0, n, [&](int i) { return a[i]; });
    }
    int size() const { return n; }
    void chmin(int l, int r, T x) {
        check(l, r);
        if (l < r) update<0>(0, 0, n, l, r, x);
    }
    void chmax(int l, int r, T x) {
        check(l, r);
        if (l < r) update<1>(0, 0, n, l, r, x);
    }
    void add(int l, int r, T x) {
        check(l, r);
        if (l < r && x) update<2>(0, 0, n, l, r, x);
    }
    void assign(int l, int r, T x) {
        check(l, r);
        if (l < r) update<3>(0, 0, n, l, r, x);
    }
    Sum sum(int l, int r) const {
        check(l, r);
        return l < r ? query<0>(0, 0, n, l, r) : Sum(0);
    }
    T min(int l, int r) const {
        check(l, r);
        return l < r ? query<1>(0, 0, n, l, r) : inf;
    }
    T max(int l, int r) const {
        check(l, r);
        return l < r ? query<2>(0, 0, n, l, r) : neg_inf;
    }
    Sum all_sum() const { return n ? tree[0].sum : Sum(0); }
    T all_min() const { return n ? tree[0].lo : inf; }
    T all_max() const { return n ? tree[0].hi : neg_inf; }
    T get(int p) const {
        assert(0 <= p && p < n);
        return min(p, p + 1);
    }
    void set(int p, T x) {
        assert(0 <= p && p < n);
        assign(p, p + 1, x);
    }
    std::vector<T> values() const {
        std::vector<T> a(n);
        if (n) collect(0, 0, n, a);
        return a;
    }
};
