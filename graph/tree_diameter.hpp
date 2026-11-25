#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>
#include <utility>
#include <vector>

namespace internal {

template <class T>
long argmax_index(const std::vector<T>& v) {
    return static_cast<long>(std::max_element(v.begin(), v.end()) - v.begin());
}

inline std::pair<std::vector<long>, std::vector<long>> bfs_tree(const std::vector<std::vector<long>>& g, long s) {
    long n = std::ssize(g);
    assert(0 <= s && s < n);
    std::vector<long> dist(n, -1);
    std::vector<long> parent(n, -1);
    std::vector<long> que;
    que.reserve(n);
    dist[s] = 0;
    que.push_back(s);
    for (long head = 0; head < std::ssize(que); ++head) {
        long v = que[head];
        for (long to : g[v]) {
            if (dist[to] != -1) continue;
            dist[to] = dist[v] + 1;
            parent[to] = v;
            que.push_back(to);
        }
    }
    assert(std::ssize(que) == n);
    return {dist, parent};
}

inline std::pair<std::vector<long>, std::vector<long>> bfs_tree(const std::vector<std::vector<std::pair<long, long>>>& g, long s) {
    long n = std::ssize(g);
    assert(0 <= s && s < n);
    std::vector<long> dist(n, -1);
    std::vector<long> parent(n, -1);
    std::vector<long> que;
    que.reserve(n);
    dist[s] = 0;
    que.push_back(s);
    for (long head = 0; head < std::ssize(que); ++head) {
        long v = que[head];
        for (auto [to, cost] : g[v]) {
            assert(cost >= 0);
            if (dist[to] != -1) continue;
            dist[to] = dist[v] + cost;
            parent[to] = v;
            que.push_back(to);
        }
    }
    assert(std::ssize(que) == n);
    return {dist, parent};
}

}  // namespace internal

// 木の直径の結果
struct TreeDiameter {
    long length;
    long u;
    long v;
    std::vector<long> path;
};

namespace internal {

template <class Graph>
TreeDiameter tree_diameter_impl(const Graph& g) {
    auto dist0 = bfs_tree(g, 0).first;
    long u = argmax_index(dist0);
    auto [du, parent] = bfs_tree(g, u);
    long v = argmax_index(du);

    std::vector<long> path;
    for (long cur = v; cur != -1; cur = parent[cur]) {
        path.push_back(cur);
    }
    std::reverse(path.begin(), path.end());
    return TreeDiameter{du[v], u, v, std::move(path)};
}

template <class Graph>
std::vector<long> tree_eccentricities_impl(const Graph& g) {
    auto dist0 = bfs_tree(g, 0).first;
    long u = argmax_index(dist0);
    auto du = bfs_tree(g, u).first;
    long v = argmax_index(du);
    auto dv = bfs_tree(g, v).first;

    long n = std::ssize(g);
    std::vector<long> ecc(n);
    for (long i = 0; i < n; ++i) {
        ecc[i] = std::max(du[i], dv[i]);
    }
    return ecc;
}

}  // namespace internal

// 頂点 s から各頂点への距離（重みなし木）
// 計算量: 時間 O(N), 空間 O(N)
// 事前条件: g は頂点数 N >= 1 の木, 0 <= s < N
inline std::vector<long> tree_distances(const std::vector<std::vector<long>>& g, long s) {
    return internal::bfs_tree(g, s).first;
}

// 頂点 s から各頂点への距離（重みあり木）
// 計算量: 時間 O(N), 空間 O(N)
// 事前条件: g は頂点数 N >= 1 の木（各辺の重みは非負）, 0 <= s < N
inline std::vector<long> tree_distances(const std::vector<std::vector<std::pair<long, long>>>& g, long s) {
    return internal::bfs_tree(g, s).first;
}

// 木の直径（長さ・両端点・単純パス）（重みなし木）
// 計算量: 時間 O(N), 空間 O(N)
// 事前条件: g は頂点数 N >= 1 の木
inline TreeDiameter tree_diameter(const std::vector<std::vector<long>>& g) {
    return internal::tree_diameter_impl(g);
}

// 木の直径（長さ・両端点・単純パス）（重みあり木）
// 計算量: 時間 O(N), 空間 O(N)
// 事前条件: g は頂点数 N >= 1 の木（各辺の重みは非負）
inline TreeDiameter tree_diameter(const std::vector<std::vector<std::pair<long, long>>>& g) {
    return internal::tree_diameter_impl(g);
}

// 木の中心（直径上の中心頂点：辺数が偶数なら 1 頂点、奇数なら 2 頂点）
// 計算量: 時間 O(N), 空間 O(N)
// 事前条件: g は頂点数 N >= 1 の木
inline std::vector<long> tree_centers(const std::vector<std::vector<long>>& g) {
    auto dia = tree_diameter(g);
    long k = std::ssize(dia.path);
    if (k % 2 == 1) return {dia.path[k / 2]};
    return {dia.path[k / 2 - 1], dia.path[k / 2]};
}

// 各頂点の偏心度（各頂点からの最遠頂点までの距離）（重みなし木）
// 計算量: 時間 O(N), 空間 O(N)
// 事前条件: g は頂点数 N >= 1 の木
inline std::vector<long> tree_eccentricities(const std::vector<std::vector<long>>& g) {
    return internal::tree_eccentricities_impl(g);
}

// 各頂点の偏心度（各頂点からの最遠頂点までの距離）（重みあり木）
// 計算量: 時間 O(N), 空間 O(N)
// 事前条件: g は頂点数 N >= 1 の木（各辺の重みは非負）
inline std::vector<long> tree_eccentricities(const std::vector<std::vector<std::pair<long, long>>>& g) {
    return internal::tree_eccentricities_impl(g);
}
