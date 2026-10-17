#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>

#include "components/splitmix64.hpp"

template <class S, auto op, auto e, class F, auto mapping, auto composition, auto id, bool commutative = false>
class LazyReversibleTreap {
    struct Empty {};
    using Back = std::conditional_t<commutative, Empty, S>;
    struct Node {
        int left = 0, right = 0;
        std::uint32_t count : 30 = 1;
        std::uint32_t pending : 1 = 0;
        std::uint32_t reversed : 1 = 0;
        std::uint32_t priority;
        S value, product;
        [[no_unique_address]] Back back;
        F lazy = id();
        static Back make_back(const S& value) {
            if constexpr (commutative)
                return {};
            else
                return value;
        }
        Node(S value, std::uint32_t priority)
            : priority(priority), value(std::move(value)), product(this->value), back(make_back(this->value)) {}
    };
    std::vector<Node> t;
    SplitMix64 rng;
    int root = 0, free = 0;
    int count(int x) const { return t[x].count; }
    void pull(int x) {
        auto& v = t[x];
        int l = v.left, r = v.right;
        v.count = count(l) + count(r) + 1;
        v.product = l ? op(t[l].product, v.value) : v.value;
        if (r) v.product = op(std::move(v.product), t[r].product);
        if constexpr (!commutative) {
            v.back = r ? op(t[r].back, v.value) : v.value;
            if (l) v.back = op(std::move(v.back), t[l].back);
        }
    }
    void flip(int x) {
        if (!x) return;
        auto& v = t[x];
        std::swap(v.left, v.right);
        if constexpr (!commutative) std::swap(v.product, v.back);
        v.reversed ^= 1;
    }
    void act(int x, const F& f) {
        if (!x) return;
        auto& v = t[x];
        v.value = mapping(f, v.value, 1);
        v.product = mapping(f, v.product, count(x));
        if constexpr (!commutative) v.back = mapping(f, v.back, count(x));
        v.lazy = v.pending ? composition(f, v.lazy) : f;
        v.pending = true;
    }
    void push(int x) {
        auto& v = t[x];
        if (v.reversed) {
            flip(v.left);
            flip(v.right);
            v.reversed = 0;
        }
        if (v.pending) {
            act(v.left, v.lazy);
            act(v.right, v.lazy);
            v.pending = false;
        }
    }
    int make(S value) {
        if (free) {
            int x = free;
            free = t[x].left;
            t[x] = Node(std::move(value), std::uint32_t(rng()));
            return x;
        }
        assert(t.size() <= (std::size_t(1) << 30) - 1);
        t.emplace_back(std::move(value), std::uint32_t(rng()));
        return int(t.size()) - 1;
    }
    std::pair<int, int> split(int x, int k) {
        if (!x) return {0, 0};
        if (!k) return {0, x};
        if (k == count(x)) return {x, 0};
        push(x);
        int left = count(t[x].left);
        if (k <= left) {
            auto [a, b] = split(t[x].left, k);
            t[x].left = b;
            pull(x);
            return {a, x};
        }
        auto [a, b] = split(t[x].right, k - left - 1);
        t[x].right = a;
        pull(x);
        return {x, b};
    }
    int merge(int a, int b) {
        if (!a || !b) return a | b;
        if (t[a].priority > t[b].priority) {
            push(a);
            t[a].right = merge(t[a].right, b);
            pull(a);
            return a;
        }
        push(b);
        t[b].left = merge(a, t[b].left);
        pull(b);
        return b;
    }
    int insert_node(int x, int p, int added) {
        if (!x) return added;
        if (t[added].priority >= t[x].priority) {
            auto [a, b] = split(x, p);
            t[added].left = a;
            t[added].right = b;
            pull(added);
            return added;
        }
        push(x);
        int left = count(t[x].left);
        if (p <= left)
            t[x].left = insert_node(t[x].left, p, added);
        else
            t[x].right = insert_node(t[x].right, p - left - 1, added);
        pull(x);
        return x;
    }
    int erase_node(int x, int p, S& value) {
        push(x);
        int left = count(t[x].left);
        if (p == left) {
            value = std::move(t[x].value);
            int result = merge(t[x].left, t[x].right);
            t[x].left = free;
            free = x;
            return result;
        }
        if (p < left)
            t[x].left = erase_node(t[x].left, p, value);
        else
            t[x].right = erase_node(t[x].right, p - left - 1, value);
        pull(x);
        return x;
    }
    void set_node(int x, int p, S value) {
        push(x);
        int left = count(t[x].left);
        if (p == left)
            t[x].value = std::move(value);
        else if (p < left)
            set_node(t[x].left, p, std::move(value));
        else
            set_node(t[x].right, p - left - 1, std::move(value));
        pull(x);
    }
    S fold(int x, int l, int r) {
        if (!l && r == count(x)) return t[x].product;
        push(x);
        int left = count(t[x].left);
        if (r <= left) return fold(t[x].left, l, r);
        if (left < l) return fold(t[x].right, l - left - 1, r - left - 1);
        S result = t[x].value;
        if (l < left) result = op(fold(t[x].left, l, left), std::move(result));
        if (left + 1 < r) result = op(std::move(result), fold(t[x].right, 0, r - left - 1));
        return result;
    }
    void apply_node(int x, int l, int r, const F& f) {
        if (!l && r == count(x)) {
            act(x, f);
            return;
        }
        push(x);
        int left = count(t[x].left);
        if (l < left) apply_node(t[x].left, l, std::min(r, left), f);
        if (l <= left && left < r) t[x].value = mapping(f, t[x].value, 1);
        if (left + 1 < r) apply_node(t[x].right, std::max(l - left - 1, 0), r - left - 1, f);
        pull(x);
    }
    int reverse_node(int x, int l, int r) {
        if (!l && r == count(x)) {
            flip(x);
            return x;
        }
        push(x);
        int left = count(t[x].left);
        if (r <= left)
            t[x].left = reverse_node(t[x].left, l, r);
        else if (left < l)
            t[x].right = reverse_node(t[x].right, l - left - 1, r - left - 1);
        else {
            auto [a, b] = split(t[x].left, l);
            auto [c, d] = split(t[x].right, r - left - 1);
            t[x].left = b;
            t[x].right = c;
            pull(x);
            flip(x);
            return merge(a, merge(x, d));
        }
        pull(x);
        return x;
    }
    void collect(int x, std::vector<S>& values) {
        if (!x) return;
        push(x);
        collect(t[x].left, values);
        values.push_back(t[x].value);
        collect(t[x].right, values);
    }
    void bounds(int l, int r) const {
        assert(0 <= l && l <= r && r <= size());
        (void)l;
        (void)r;
    }

   public:
    using value_type = S;
    using action_type = F;
    explicit LazyReversibleTreap(std::uint64_t seed = 0) : rng(seed) {
        t.emplace_back(e(), 0);
        t[0].count = 0;
    }
    explicit LazyReversibleTreap(const std::vector<S>& values, std::uint64_t seed = 0) : LazyReversibleTreap(seed) { assign(values); }
    LazyReversibleTreap(const LazyReversibleTreap&) = default;
    LazyReversibleTreap& operator=(const LazyReversibleTreap&) = default;
    LazyReversibleTreap(LazyReversibleTreap&& other) : LazyReversibleTreap() { swap(other); }
    LazyReversibleTreap& operator=(LazyReversibleTreap&& other) {
        if (this != &other) {
            LazyReversibleTreap moved(std::move(other));
            swap(moved);
        }
        return *this;
    }
    void swap(LazyReversibleTreap& other) noexcept {
        t.swap(other.t);
        std::swap(rng, other.rng);
        std::swap(root, other.root);
        std::swap(free, other.free);
    }
    int size() const { return count(root); }
    bool empty() const { return !root; }
    void reserve(std::size_t n) {
        assert(n < (std::size_t(1) << 30));
        t.reserve(n + 1);
    }
    std::size_t memory_usage() const { return t.capacity() * sizeof(Node); }
    void clear() {
        t.erase(t.begin() + 1, t.end());
        root = free = 0;
    }
    void assign(const std::vector<S>& values) {
        clear();
        reserve(values.size());
        std::vector<int> path;
        for (const S& value : values) {
            int x = make(value), last = 0;
            while (!path.empty() && t[path.back()].priority <= t[x].priority) {
                last = path.back();
                path.pop_back();
                pull(last);
            }
            t[x].left = last;
            if (!path.empty()) t[path.back()].right = x;
            path.push_back(x);
        }
        if (!path.empty()) root = path[0];
        while (!path.empty()) {
            pull(path.back());
            path.pop_back();
        }
    }
    void insert(int p, S value) {
        bounds(p, p);
        root = insert_node(root, p, make(std::move(value)));
    }
    S erase(int p) {
        bounds(p, p + 1);
        S value = e();
        root = erase_node(root, p, value);
        return value;
    }
    S get(int p) {
        bounds(p, p + 1);
        int x = root;
        for (;;) {
            push(x);
            int left = count(t[x].left);
            if (p == left) return t[x].value;
            if (p < left)
                x = t[x].left;
            else
                p -= left + 1, x = t[x].right;
        }
    }
    void set(int p, S value) {
        bounds(p, p + 1);
        set_node(root, p, std::move(value));
    }
    S prod(int l, int r) {
        bounds(l, r);
        return l == r ? e() : fold(root, l, r);
    }
    S all_prod() const { return t[root].product; }
    void apply(int p, const F& f) { apply(p, p + 1, f); }
    void apply(int l, int r, const F& f) {
        bounds(l, r);
        if (l != r) apply_node(root, l, r, f);
    }
    void reverse(int l, int r) {
        bounds(l, r);
        if (r - l > 1) root = reverse_node(root, l, r);
    }
    void rotate(int l, int m, int r) {
        bounds(l, m);
        bounds(m, r);
        if (l == m || m == r) return;
        auto [a, b] = split(root, l);
        auto [c, d] = split(b, m - l);
        auto [f, g] = split(d, r - m);
        root = merge(a, merge(f, merge(c, g)));
    }
    void move(int l, int r, int p) {
        bounds(l, r);
        assert(0 <= p && p <= size() - (r - l));
        if (l == r || l == p) return;
        auto [a, b] = split(root, l);
        auto [c, d] = split(b, r - l);
        auto [f, g] = split(merge(a, d), p);
        root = merge(f, merge(c, g));
    }
    std::vector<S> to_vector() {
        std::vector<S> values;
        values.reserve(size());
        collect(root, values);
        return values;
    }
};
