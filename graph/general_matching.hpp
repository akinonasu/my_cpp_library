
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <utility>
#include <vector>

class GeneralMatching {
    int n, cardinality = 0, tail = 0;
    std::uint32_t stamp = 0;
    std::vector<std::pair<int, int>> input;
    std::vector<int> start, adj, base, parent, queue;
    std::vector<std::uint32_t> seen;
    std::vector<signed char> state;
    int find(int v) {
        while (base[v] != v)
            v = base[v] = base[base[v]];
        return v;
    }
    int lca(int u, int v) {
        if (!++stamp) std::fill(seen.begin(), seen.end(), 0), ++stamp;
        for (;; std::swap(u, v))
            if (u != -1) {
                u = find(u);
                if (seen[u] == stamp) return u;
                seen[u] = stamp;
                u = mate[u] == -1 ? -1 : parent[mate[u]];
            }
    }
    void contract(int u, int v, int b) {
        while (find(u) != b) {
            parent[u] = v;
            v = mate[u];
            if (state[v] == 1) state[v] = 0, queue[tail++] = v;
            if (base[u] == u) base[u] = b;
            if (base[v] == v) base[v] = b;
            u = parent[v];
        }
    }
    void flip(int v) {
        while (v != -1) {
            int u = parent[v], next = mate[u];
            mate[u] = v;
            mate[v] = u;
            v = next;
        }
    }
    bool augment(int root) {
        std::iota(base.begin(), base.end(), 0);
        std::fill(state.begin(), state.end(), -1);
        tail = 1;
        queue[0] = root;
        state[root] = 0;
        for (int head = 0; head < tail; ++head) {
            int u = queue[head];
            for (int e = start[u]; e < start[u + 1]; ++e) {
                int v = adj[e];
                if (state[v] == -1) {
                    parent[v] = u;
                    state[v] = 1;
                    if (mate[v] == -1) {
                        flip(v);
                        return true;
                    }
                    state[mate[v]] = 0;
                    queue[tail++] = mate[v];
                } else if (state[v] == 0 && find(u) != find(v)) {
                    int b = lca(u, v);
                    contract(u, v, b);
                    contract(v, u, b);
                }
            }
        }
        return false;
    }
    void short_paths() {
        std::vector<int> first(n, -1), second(n, -1);
        auto save = [&](int u, int v) {
            if (first[u] == -1)
                first[u] = v;
            else if (first[u] != v)
                second[u] = v;
        };
        for (auto [u, v] : input) {
            if (mate[u] == -1) save(v, u);
            if (mate[v] == -1) save(u, v);
        }
        for (int u = 0; u < n; ++u)
            if (u < mate[u]) {
                int v = mate[u];
                for (int a : {first[u], second[u]})
                    if (a != -1 && mate[a] == -1) {
                        int b = first[v];
                        if (b == -1 || b == a || mate[b] != -1) b = second[v];
                        if (b == -1 || b == a || mate[b] != -1) continue;
                        mate[u] = a;
                        mate[a] = u;
                        mate[v] = b;
                        mate[b] = v;
                        ++cardinality;
                        break;
                    }
            }
    }
    int component_bound(std::vector<int>& roots) {
        std::fill(parent.begin(), parent.end(), -1);
        int bound = 0;
        for (int root = 0; root < n; ++root)
            if (parent[root] == -1) {
                int end = 1, count[2] = {1, 0};
                bool bipartite = true;
                queue[0] = root;
                parent[root] = 0;
                for (int head = 0; head < end; ++head) {
                    int u = queue[head];
                    for (int e = start[u]; e < start[u + 1]; ++e) {
                        int v = adj[e];
                        if (parent[v] == -1) {
                            parent[v] = parent[u] ^ 1;
                            ++count[parent[v]];
                            queue[end++] = v;
                        } else if (parent[v] == parent[u])
                            bipartite = false;
                    }
                }
                bound += bipartite ? std::min(count[0], count[1]) : end / 2;
                int side = count[0] > count[1];
                for (int i = 0; i < end; ++i) {
                    int u = queue[i];
                    if (mate[u] == -1 && (!bipartite || parent[u] == side)) roots.push_back(u);
                }
            }
        return bound;
    }

   public:
    std::vector<int> mate;
    explicit GeneralMatching(int n) : n(n), mate(n, -1) { assert(n >= 0); }
    void reserve(std::size_t m) { input.reserve(m); }
    void add_edge(int u, int v) {
        assert(0 <= u && u < n && 0 <= v && v < n);
        if (u == v || cardinality == n / 2) return;
        input.emplace_back(u, v);
        if (mate[u] == -1 && mate[v] == -1) {
            mate[u] = v;
            mate[v] = u;
            ++cardinality;
        }
    }
    int build() {
        if (input.empty() || cardinality == n / 2) return cardinality;
        short_paths();
        if (cardinality == n / 2) return cardinality;
        start.assign(n + 1, 0);
        for (auto [u, v] : input)
            ++start[u + 1], ++start[v + 1];
        std::partial_sum(start.begin(), start.end(), start.begin());
        adj.resize(input.size() * 2);
        parent = start;
        for (auto [u, v] : input)
            adj[parent[u]++] = v, adj[parent[v]++] = u;
        parent.resize(n);
        queue.resize(n);
        std::vector<int> roots;
        int limit = component_bound(roots);
        if (cardinality == limit) return cardinality;
        base.resize(n);
        seen.resize(n);
        state.resize(n);
        for (int u : roots)
            if (mate[u] == -1 && augment(u)) {
                if (++cardinality == limit) break;
            }
        return cardinality;
    }
    std::vector<std::pair<int, int>> edges() const {
        std::vector<std::pair<int, int>> result;
        result.reserve(cardinality);
        for (int u = 0; u < n; ++u)
            if (u < mate[u]) result.emplace_back(u, mate[u]);
        return result;
    }
};
