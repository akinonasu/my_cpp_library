#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

struct LowLink {
    struct Arc {
        int to, id;
    };
    struct Frame {
        int v, parent, next;
    };
    int n;
    std::vector<std::pair<int, int>> edges;
    std::vector<int> ord, low, articulation, bridges;
    explicit LowLink(int n) : n(n) { assert(n >= 0); }
    int add_edge(int u, int v) {
        assert(0 <= u && u < n && 0 <= v && v < n);
        assert(edges.size() < std::size_t(std::numeric_limits<int>::max() / 2));
        int id = static_cast<int>(edges.size());
        edges.emplace_back(u, v);
        return id;
    }
    void build() {
        assert(n >= 0 && edges.size() <= std::size_t(std::numeric_limits<int>::max() / 2));
        std::vector<int> start(std::size_t(n) + 1);
        for (auto [u, v] : edges) {
            assert(0 <= u && u < n && 0 <= v && v < n);
            ++start[u + 1], ++start[v + 1];
        }
        for (int v = 0; v < n; ++v)
            start[v + 1] += start[v];
        std::unique_ptr<Arc[]> adj(new Arc[2 * edges.size()]);
        low.assign(start.begin(), start.end() - 1);
        for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
            auto [u, v] = edges[i];
            adj[low[u]++] = {v, i};
            adj[low[v]++] = {u, i};
        }
        ord.assign(n, -1);
        low.assign(n, -1);
        articulation.clear();
        bridges.clear();
        std::vector<std::uint8_t> cut(n);
        std::vector<Frame> stack;
        stack.reserve(n);
        int timer = 0;
        for (int root = 0; root < n; ++root) {
            if (ord[root] != -1) continue;
            int children = 0;
            ord[root] = low[root] = timer++;
            stack.push_back({root, -1, start[root]});
            while (!stack.empty()) {
                auto& f = stack.back();
                int v = f.v;
                if (f.next == start[v + 1]) {
                    int id = f.parent;
                    stack.pop_back();
                    if (stack.empty())
                        cut[v] = children > 1;
                    else {
                        int p = stack.back().v;
                        low[p] = std::min(low[p], low[v]);
                        if (ord[p] < low[v]) bridges.push_back(id);
                        if (p != root && ord[p] <= low[v]) cut[p] = 1;
                    }
                    continue;
                }
                auto [to, id] = adj[f.next++];
                if (id == f.parent) continue;
                if (ord[to] != -1) {
                    low[v] = std::min(low[v], ord[to]);
                    continue;
                }
                children += v == root;
                ord[to] = low[to] = timer++;
                stack.push_back({to, id, start[to]});
            }
        }
        for (int v = 0; v < n; ++v)
            if (cut[v]) articulation.push_back(v);
    }
    bool is_bridge(int id) const {
        assert(0 <= id && id < static_cast<int>(edges.size()) && static_cast<int>(ord.size()) == n);
        auto [u, v] = edges[id];
        if (ord[u] > ord[v]) std::swap(u, v);
        return ord[u] < low[v];
    }
};
