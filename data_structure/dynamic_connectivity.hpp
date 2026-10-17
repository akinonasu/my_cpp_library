#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "components/splitmix64.hpp"

template <class S = void, auto op = nullptr, auto e = nullptr>
class DynamicConnectivity {
    static constexpr bool aggregate = !std::is_void_v<S>;
    struct Empty {};
    using Value = std::conditional_t<aggregate, S, Empty>;
    struct alignas(32) Node {
        int ch[2] = {}, parent = 0, size = 0, key = 0, aux = 0;
        std::uint8_t own = 0, flags = 0;
    };
    struct Edge {
        std::uint64_t key = 0;
        int top = 0, count = 0, prev[2] = {}, next[2] = {};
    };
    int n = 0, groups = 0, active = 0, free_edge = 0, free_arc = 0;
    std::int64_t m = 0;
    std::uint64_t seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::vector<Node> t;
    std::vector<Edge> edges;
    std::vector<int> table;
    std::vector<std::vector<int>> vertices;
    std::vector<Value> values, sums;
    static std::uint64_t key(int u, int v) {
        if (u > v) std::swap(u, v);
        return (std::uint64_t(unsigned(u)) << 32) | unsigned(v);
    }
    std::uint64_t hash(std::uint64_t x) const { return splitmix64(x + seed); }
    int level_of(int id) const {
        int x = edges[id].top;
        return x > 0 ? (t[x].own >> 3) : -x - 1;
    }
    int slot(std::uint64_t k) const {
        int mask = static_cast<int>(table.size()) - 1, i = hash(k) & mask;
        while (table[i] && edges[table[i]].key != k)
            i = (i + 1) & mask;
        return i;
    }
    void rehash(int capacity) {
        table.assign(capacity, 0);
        for (int i = 1; i < static_cast<int>(edges.size()); ++i)
            if (edges[i].count) table[slot(edges[i].key)] = i;
    }
    void erase_slot(int i) {
        int mask = static_cast<int>(table.size()) - 1;
        for (int j = (i + 1) & mask; table[j]; j = (j + 1) & mask) {
            int k = hash(edges[table[j]].key) & mask;
            if (((j - k) & mask) >= ((j - i) & mask)) table[i] = table[j], i = j;
        }
        table[i] = 0;
    }
    void check(int v) const {
        assert(0 <= v && v < n);
        (void)v;
    }
    [[gnu::always_inline]] void pull(int x) {
        auto& v = t[x];
        int l = v.ch[0], r = v.ch[1];
        v.size = t[l].size + t[r].size + (v.key > 0);
        v.flags = t[l].flags | t[r].flags | (v.own & 3);
        if constexpr (aggregate) {
            if (v.own & 4) sums[x] = op(op(sums[l], v.key > 0 ? values[v.key - 1] : sums[0]), sums[r]);
        }
    }
    [[gnu::always_inline]] void rotate(int x) {
        int p = t[x].parent, g = t[p].parent;
        int d = t[p].ch[1] == x, c = t[x].ch[d ^ 1];
        if (g) t[g].ch[t[g].ch[1] == p] = x;
        t[x].parent = g;
        t[p].ch[d] = c;
        if (c) t[c].parent = p;
        t[x].ch[d ^ 1] = p;
        t[p].parent = x;
        pull(p);
    }
    template <bool update = true>
    [[gnu::always_inline]] void splay(int x) {
        if (!t[x].parent) return;
        while (t[x].parent) {
            int p = t[x].parent, g = t[p].parent;
            if (!g) {
                rotate(x);
                break;
            }
            int gg = t[g].parent, d = t[p].ch[1] == x, side = t[g].ch[1] == p;
            if (d == side) {
                int a = t[p].ch[d ^ 1], b = t[x].ch[d ^ 1];
                t[g].ch[d] = a;
                if (a) t[a].parent = g;
                t[p].ch[d] = b;
                if (b) t[b].parent = p;
                t[p].ch[d ^ 1] = g;
                t[g].parent = p;
                t[x].ch[d ^ 1] = p;
                t[p].parent = x;
                pull(g);
                pull(p);
            } else {
                int a = t[x].ch[d], b = t[x].ch[d ^ 1];
                t[p].ch[d] = b;
                if (b) t[b].parent = p;
                t[g].ch[side] = a;
                if (a) t[a].parent = g;
                t[x].ch[d] = g;
                t[x].ch[d ^ 1] = p;
                t[p].parent = t[g].parent = x;
                pull(p);
                pull(g);
            }
            t[x].parent = gg;
            if (gg) t[gg].ch[t[gg].ch[1] == g] = x;
        }
        if constexpr (update) pull(x);
    }
    bool same_nodes(int u, int v) {
        if (u == v) return true;
        splay(u);
        splay(v);
        return t[u].parent != 0;
    }
    int join(int a, int b) {
        if (!a || !b) return a ? a : b;
        int x = a;
        while (t[x].ch[1])
            x = t[x].ch[1];
        splay<false>(x);
        t[x].ch[1] = b;
        t[b].parent = x;
        pull(x);
        return x;
    }
    int reroot(int x) {
        splay<false>(x);
        int l = t[x].ch[0];
        if (!l) {
            pull(x);
            return x;
        }
        t[l].parent = t[x].ch[0] = 0;
        pull(x);
        return join(x, l);
    }
    int vertex(int level, int v) {
        if (level == static_cast<int>(vertices.size())) vertices.emplace_back(n, 0);
        int& x = vertices[level][v];
        if (!x) {
            assert(t.size() < std::size_t(std::numeric_limits<int>::max()));
            x = static_cast<int>(t.size());
            t.emplace_back();
            t[x].key = v + 1;
            t[x].size = 1;
            if constexpr (aggregate) sums.push_back(e());
        }
        return x;
    }
    int allocate_arcs() {
        int x = free_arc;
        if (x) {
            free_arc = t[x].ch[0];
            t[x] = t[x + 1] = Node{};
            if constexpr (aggregate) sums[x] = sums[x + 1] = e();
        } else {
            assert(t.size() <= std::size_t(std::numeric_limits<int>::max() - 2));
            x = static_cast<int>(t.size());
            t.emplace_back();
            t.emplace_back();
            if constexpr (aggregate) sums.insert(sums.end(), 2, e());
        }
        return x;
    }
    void link_tree(int id, int level, bool exact) {
        int u = vertex(level, edges[id].key >> 32), v = vertex(level, unsigned(edges[id].key));
        u = reroot(u), v = reroot(v);
        int a = allocate_arcs(), b = a + 1;
        t[a].key = -id;
        t[a].aux = edges[id].top;
        t[a].own = (level << 3) | (aggregate && level == 0 ? 4 : 0) | exact;
        t[b].own = aggregate && level == 0 ? 4 : 0;
        t[a].ch[0] = u, t[a].ch[1] = v;
        t[u].parent = t[v].parent = a;
        t[a].size = t[u].size + t[v].size;
        t[a].flags = t[u].flags | t[v].flags | exact;
        if constexpr (aggregate)
            if (level == 0) sums[a] = op(sums[u], sums[v]);
        t[b].ch[0] = a;
        t[a].parent = b;
        t[b].size = t[a].size;
        t[b].flags = t[a].flags;
        if constexpr (aggregate)
            if (level == 0) sums[b] = sums[a];
        edges[id].top = a;
    }
    std::pair<int, int> cut_tree(int a) {
        int b = a + 1;
        splay<false>(a);
        int l = t[a].ch[0], r = t[a].ch[1];
        if (l) t[l].parent = 0;
        if (r) t[r].parent = 0;
        join(r, l);
        splay<false>(b);
        l = t[b].ch[0], r = t[b].ch[1];
        if (l) t[l].parent = 0;
        if (r) t[r].parent = 0;
        t[a].ch[0] = free_arc;
        free_arc = a;
        return {l, r};
    }
    int find_flag(int x, int flag) {
        splay(x);
        if (!(t[x].flags & flag)) return 0;
        while (!(t[x].own & flag)) {
            int l = t[x].ch[0];
            x = t[l].flags & flag ? l : t[x].ch[1];
        }
        splay(x);
        return x;
    }
    void attach(int id, int level) {
        for (int side = 0; side < 2; ++side) {
            int x = vertex(level, side ? unsigned(edges[id].key) : edges[id].key >> 32);
            int old = t[x].aux, token = id * 2 + side;
            edges[id].prev[side] = 0;
            edges[id].next[side] = old;
            if (old) edges[old >> 1].prev[old & 1] = token;
            t[x].aux = token;
            if (!old) splay<false>(x), t[x].own |= 2, pull(x);
        }
        edges[id].top = -level - 1;
    }
    void detach(int id) {
        int level = level_of(id);
        for (int side = 0; side < 2; ++side) {
            int x = vertices[level][side ? unsigned(edges[id].key) : edges[id].key >> 32];
            int p = edges[id].prev[side], q = edges[id].next[side];
            if (p)
                edges[p >> 1].next[p & 1] = q;
            else
                t[x].aux = q;
            if (q) edges[q >> 1].prev[q & 1] = p;
            if (!t[x].aux) splay<false>(x), t[x].own &= ~2, pull(x);
        }
    }
    bool reconnect(const std::pair<int, int>* roots, int level) {
        for (; level >= 0; --level) {
            auto [a, b] = roots[level];
            if (t[a].size > t[b].size) std::swap(a, b);
            if (!(t[a].flags & 2)) continue;
            bool promoted = false;
            while (int x = find_flag(a, 2)) {
                while (t[x].aux) {
                    int id = t[x].aux >> 1;
                    int s = edges[id].key >> 32, d = unsigned(edges[id].key);
                    if (!same_nodes(vertices[level][s], vertices[level][d])) {
                        detach(id);
                        edges[id].top = 0;
                        for (int i = 0; i < level + 1; ++i)
                            link_tree(id, i, i == level);
                        return true;
                    }
                    if (!promoted) {
                        int token = edges[id].next[t[x].aux & 1];
                        if (token) {
                            int candidate = token >> 1;
                            int c = edges[candidate].key >> 32, d = unsigned(edges[candidate].key);
                            if (!same_nodes(vertices[level][c], vertices[level][d])) {
                                detach(candidate);
                                edges[candidate].top = 0;
                                for (int i = 0; i < level + 1; ++i)
                                    link_tree(candidate, i, i == level);
                                return true;
                            }
                        }

                        while (int y = find_flag(a, 1)) {
                            int edge = -t[y].key;
                            assert(edges[edge].top == y && level_of(edge) == level);
                            t[y].own &= ~1;
                            pull(y);
                            link_tree(edge, level + 1, true);
                        }
                        promoted = true;
                    }
                    detach(id);
                    attach(id, level + 1);
                }
            }
        }
        return false;
    }
    void init(int count) {
        assert(0 <= count && count < std::numeric_limits<int>::max() / 3);
        n = groups = count;
        t.reserve(std::size_t(n) * 3 + 1);
        t.resize(std::size_t(n) + 1);
        vertices.emplace_back(n);
        for (int v = 0; v < n; ++v) {
            vertices[0][v] = v + 1;
            t[v + 1].key = v + 1;
            t[v + 1].size = 1;
            t[v + 1].own = aggregate ? 4 : 0;
        }
        edges.emplace_back();
        table.assign(8, 0);
        if constexpr (aggregate) {
            values.assign(n, e());
            sums.reserve(std::size_t(n) * 3 + 1);
            sums.assign(std::size_t(n) + 1, e());
        }
    }

   public:
    explicit DynamicConnectivity(int n = 0) { init(n); }
    explicit DynamicConnectivity(const std::vector<Value>& a)
        requires aggregate
    {
        assert(a.size() < std::size_t(std::numeric_limits<int>::max() / 3));
        init(static_cast<int>(a.size()));
        values = a;
        for (int i = 0; i < n; ++i)
            sums[i + 1] = a[i];
    }
    void reserve(int capacity) {
        assert(0 <= capacity && capacity < std::numeric_limits<int>::max() / 4);
        edges.reserve(std::size_t(capacity) + 1);
        int buckets = std::bit_ceil(unsigned(std::max(8, capacity * 2)));
        if (buckets > static_cast<int>(table.size())) rehash(buckets);
    }
    int size() const { return n; }
    int size(int v) {
        check(v);
        splay(v + 1);
        return t[v + 1].size;
    }
    int component_count() const { return groups; }
    std::int64_t edge_count() const { return m; }
    int count(int u, int v) const {
        check(u);
        check(v);
        return edges[table[slot(key(u, v))]].count;
    }
    bool same(int u, int v) {
        check(u);
        check(v);
        return u == v || groups == 1 || (groups != n && same_nodes(u + 1, v + 1));
    }
    bool link(int u, int v) {
        check(u);
        check(v);
        std::uint64_t k = key(u, v);
        int pos = slot(k), id = table[pos];
        ++m;
        if (id) {
            assert(edges[id].count < std::numeric_limits<int>::max());
            ++edges[id].count;
            return false;
        }
        if ((active + 1) * 2 > static_cast<int>(table.size())) rehash(static_cast<int>(table.size()) * 2), pos = slot(k);
        id = free_edge;
        if (id)
            free_edge = edges[id].next[0], edges[id] = Edge{};
        else {
            assert(edges.size() < std::size_t(std::numeric_limits<int>::max() / 2));
            id = static_cast<int>(edges.size());
            edges.emplace_back();
        }
        table[pos] = id;
        edges[id].key = k;
        edges[id].count = 1;
        ++active;
        if (u == v) return false;
        if (groups != 1 && !same_nodes(u + 1, v + 1)) {
            link_tree(id, 0, true);
            --groups;
            return true;
        }
        attach(id, 0);
        return false;
    }
    bool cut(int u, int v) {
        check(u);
        check(v);
        int pos = slot(key(u, v)), id = table[pos];
        if (!id) return false;
        --m;
        if (--edges[id].count) return false;
        erase_slot(pos);
        --active;
        bool split = false;
        if (edges[id].top > 0) {
            int a = edges[id].top;
            int level = level_of(id), i = level;
            std::array<std::pair<int, int>, 32> roots;
            while (a) {
                int down = t[a].aux;
                roots[i--] = cut_tree(a);
                a = down;
            }
            split = !reconnect(roots.data(), level);
            groups += split;
        } else if (u != v)
            detach(id);
        edges[id].next[0] = free_edge;
        free_edge = id;
        return split;
    }
    Value get(int v) const
        requires aggregate
    {
        check(v);
        return values[v];
    }
    void set(int v, Value value)
        requires aggregate
    {
        check(v);
        splay<false>(v + 1);
        values[v] = std::move(value);
        pull(v + 1);
    }
    Value prod(int v)
        requires aggregate
    {
        check(v);
        splay(v + 1);
        return sums[v + 1];
    }
};
