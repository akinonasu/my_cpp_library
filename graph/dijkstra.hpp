#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <queue>
#include <type_traits>
#include <utility>
#include <vector>

#include "components/radixheap.hpp"

inline constexpr std::int64_t DIJKSTRA_INF = 1'000'000'000'000'000'000;  // 1e18

template <class T = std::int64_t>
class Dijkstra {
   public:
    struct Edge {
        int to;
        T cost;
        int id;  // 辺番号（必要に応じて経路復元等で使用）
    };

    int n;
    T inf;
    std::vector<std::vector<Edge>> g;
    std::vector<T> dist;
    std::vector<int> prev_vertex;  // 経路復元用の直前頂点
    std::vector<int> prev_edge;    // 経路復元用の直前辺ID

    static constexpr T default_inf() {
        if constexpr (std::is_floating_point_v<T>) {
            return std::numeric_limits<T>::infinity();
        } else {
            return static_cast<T>(std::min<std::int64_t>(DIJKSTRA_INF, std::numeric_limits<T>::max()));
        }
    }

    explicit Dijkstra(int n, T inf = default_inf()) : n(n), inf(inf), g(n), dist(n, inf), prev_vertex(n, -1), prev_edge(n, -1) {
        assert(n >= 0);
        assert(inf > 0);
    }

    // 有向辺の追加 (id: 辺番号を持たせたい場合に指定)
    void add_edge(int from, int to, T cost, int id = -1) {
        assert(0 <= from && from < n);
        assert(0 <= to && to < n);
        assert(cost >= 0);
        g[from].push_back(Edge{to, cost, id});
    }

    void add_directed_edge(int from, int to, T cost, int id = -1) { add_edge(from, to, cost, id); };

    // 無向辺の追加
    void add_undirected_edge(int u, int v, T cost, int id = -1) {
        add_edge(u, v, cost, id);
        add_edge(v, u, cost, id);
    }

    // 単一始点ダイクストラの実行
    void build(int st) { build(std::vector<int>{st}); }

    // 複数始点ダイクストラの実行
    void build(const std::vector<int>& starts) {
        dist.assign(n, inf);
        prev_vertex.assign(n, -1);
        prev_edge.assign(n, -1);

        using State = std::pair<T, int>;
        constexpr bool radix = std::is_integral_v<T> && sizeof(T) <= 8;
        std::conditional_t<radix, RadixHeap<T>, std::priority_queue<State, std::vector<State>, std::greater<State>>> pq;

        for (int st : starts) {
            assert(0 <= st && st < n);
            dist[st] = 0;
            pq.push({0, st});
        }

        while (!pq.empty()) {
            const auto [cost, v] = [&] {
                if constexpr (radix) {
                    return pq.pop_value();
                } else {
                    State s = pq.top();
                    pq.pop();
                    return s;
                }
            }();

            if (dist[v] < cost) continue;

            for (const auto& [to, c, id] : g[v]) {
                if (c >= inf - cost) continue;
                T next = cost + c;
                if (next < dist[to]) {
                    dist[to] = next;
                    prev_vertex[to] = v;
                    prev_edge[to] = id;
                    pq.push({next, to});
                }
            }
        }
    }

    // 距離の取得（直接アクセス）
    T operator[](int v) const {
        assert(0 <= v && v < n);
        return dist[v];
    }

    // 到達判定
    bool reachable(int v) const {
        assert(0 <= v && v < n);
        return dist[v] != inf;
    }

    // 到達不能時に unreachable_val (デフォルト -1) を返す取得メソッド
    T get_dist(int v, T unreachable_val = -1) const {
        assert(0 <= v && v < n);
        return reachable(v) ? dist[v] : unreachable_val;
    }

    // 全頂点の距離ベクトルを取得（到達不能を -1 に置換する旧互換用）
    std::vector<T> get_distances(T unreachable_val = -1) const {
        std::vector<T> res = dist;
        for (auto& val : res) {
            if (val == inf) val = unreachable_val;
        }
        return res;
    }

    // 頂点の経路復元 (始点 -> ... -> goal)
    std::vector<int> restore_path(int goal) const {
        assert(0 <= goal && goal < n);
        if (!reachable(goal)) return {};
        std::vector<int> path;
        for (int v = goal; v != -1; v = prev_vertex[v]) {
            path.push_back(v);
        }
        std::ranges::reverse(path);
        return path;
    }

    // 辺IDの経路復元 (使用した辺のID列)
    std::vector<int> restore_edge_path(int goal) const {
        assert(0 <= goal && goal < n);
        if (!reachable(goal)) return {};
        std::vector<int> path;
        for (int v = goal; prev_edge[v] != -1; v = prev_vertex[v]) {
            path.push_back(prev_edge[v]);
        }
        std::ranges::reverse(path);
        return path;
    }
};

// 従来の関数インターフェース（後方互換用）
template <class T>
std::vector<T> dijkstra(const std::vector<std::vector<std::pair<int, T>>>& g, int st, T inf = Dijkstra<T>::default_inf()) {
    Dijkstra<T> dij(static_cast<int>(g.size()), inf);
    for (int u = 0; u < static_cast<int>(g.size()); ++u) {
        for (const auto& [v, c] : g[u]) {
            dij.add_edge(u, v, c);
        }
    }
    dij.build(st);
    return dij.get_distances(-1);
}