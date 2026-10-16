#pragma once
#include <cassert>
#include <cstddef>
#include <limits>
#include <optional>
#include <vector>

// Returns nullopt for a cyclic graph, including a self-loop.
inline std::optional<std::vector<int>> topological_sort(const std::vector<std::vector<int>>& graph) {
    assert(graph.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()));
    const int n = static_cast<int>(graph.size());
    std::vector<std::size_t> indegree(n);
    for (const auto& row : graph)
        for (int v : row) {
            assert(0 <= v && v < n);
            ++indegree[v];
        }
    std::vector<int> order;
    order.reserve(n);
    for (int v = 0; v < n; ++v)
        if (!indegree[v]) order.push_back(v);
    for (std::size_t i = 0; i < order.size(); ++i)
        for (int v : graph[order[i]])
            if (--indegree[v] == 0) order.push_back(v);
    if (order.size() != graph.size()) return std::nullopt;
    return order;
}