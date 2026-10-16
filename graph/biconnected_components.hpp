#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

/** @brief Vertex-biconnected blocks of a loop-free undirected multigraph. */
class BiconnectedComponents {
   public:
    BiconnectedComponents(int n, const std::vector<std::pair<int, int>>& edges) : n_(n) {
        assert(n >= 0);
        assert(edges.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()));
        std::vector<std::size_t> offset(static_cast<std::size_t>(n) + 1);
        for (const auto& [u, v] : edges) {
            assert(0 <= u && u < n && 0 <= v && v < n && u != v);
            ++offset[u + 1];
            ++offset[v + 1];
        }
        for (int v = 0; v < n; ++v)
            offset[v + 1] += offset[v];
        std::vector<int> adjacency(offset.back());
        auto next = offset;
        for (int id = 0; id < static_cast<int>(edges.size()); ++id) {
            const auto [u, v] = edges[id];
            adjacency[next[u]++] = id;
            adjacency[next[v]++] = id;
        }
        next = offset;
        std::vector<int> order(n, -1), low(n), parent(n, -1), parent_edge(n, -1), stamp(n, -1);
        std::vector<int> stack, edge_stack;
        stack.reserve(n);
        edge_stack.reserve(edges.size());
        int timer = 0;
        for (int root = 0; root < n; ++root) {
            if (order[root] != -1) continue;
            order[root] = low[root] = timer++;
            if (offset[root] == offset[root + 1]) {
                groups_.push_back({root});
                continue;
            }
            stack.push_back(root);
            while (!stack.empty()) {
                const int v = stack.back();
                if (next[v] < offset[v + 1]) {
                    const int id = adjacency[next[v]++];
                    if (id == parent_edge[v]) continue;
                    const auto [a, b] = edges[id];
                    const int to = a ^ b ^ v;
                    if (order[to] == -1) {
                        parent[to] = v;
                        parent_edge[to] = id;
                        edge_stack.push_back(id);
                        order[to] = low[to] = timer++;
                        stack.push_back(to);
                    } else if (order[to] < order[v]) {
                        edge_stack.push_back(id);
                        low[v] = std::min(low[v], order[to]);
                    }
                } else {
                    stack.pop_back();
                    const int p = parent[v];
                    if (p == -1) continue;
                    low[p] = std::min(low[p], low[v]);
                    if (low[v] >= order[p]) {
                        const int block = static_cast<int>(groups_.size());
                        groups_.emplace_back();
                        int id;
                        do {
                            id = edge_stack.back();
                            edge_stack.pop_back();
                            const auto [a, b] = edges[id];
                            for (int vertex : {a, b})
                                if (stamp[vertex] != block) {
                                    stamp[vertex] = block;
                                    groups_.back().push_back(vertex);
                                }
                        } while (id != parent_edge[v]);
                    }
                }
            }
        }
    }

    const std::vector<std::vector<int>>& groups() const { return groups_; }

    std::vector<std::vector<int>> block_cut_tree() const {
        assert(groups_.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max() - n_));
        std::vector<std::vector<int>> result(n_ + groups_.size());
        for (std::size_t block = 0; block < groups_.size(); ++block) {
            const int node = n_ + static_cast<int>(block);
            result[node] = groups_[block];
            for (int v : groups_[block])
                result[v].push_back(node);
        }
        return result;
    }

   private:
    int n_;
    std::vector<std::vector<int>> groups_;
};
