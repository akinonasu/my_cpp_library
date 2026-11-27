#pragma once

#include <cassert>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

namespace internal {

inline long edge_to(long to) {
    return to;
}

inline long edge_to(const std::pair<long, long>& e) {
    return e.first;
}

template <class PutEdge, class T>
inline T apply_edge(PutEdge put_edge, const T& val, long p, long to) {
    return put_edge(val, p, to);
}

template <class PutEdge, class T>
inline T apply_edge(PutEdge put_edge, const T& val, long p, const std::pair<long, long>& e) {
    return put_edge(val, p, e.first, e.second);
}

template <class Edge>
inline std::pair<std::vector<long>, std::vector<long>> tree_order_impl(const std::vector<std::vector<Edge>>& g, long root) {
    long n = std::ssize(g);
    assert(0 <= root && root < n);
    std::vector<long> order;
    order.reserve(n);
    std::vector<long> parent(n, -1);
    order.push_back(root);
    for (long head = 0; head < std::ssize(order); ++head) {
        long u = order[head];
        for (const auto& e : g[u]) {
            long v = edge_to(e);
            if (v == parent[u]) continue;
            assert(v != root && parent[v] == -1);
            parent[v] = u;
            order.push_back(v);
        }
    }
    assert(std::ssize(order) == n);
    return {order, parent};
}

template <class Edge, class E, class Merge, class PutEdge, class PutVertex>
inline auto tree_dp_impl(const std::vector<std::vector<Edge>>& g, long root, E e, Merge merge, PutEdge put_edge, PutVertex put_vertex) {
    using T = std::decay_t<decltype(e())>;
    auto [order, parent] = tree_order_impl(g, root);
    long n = std::ssize(g);
    std::vector<T> dp(n, e());
    for (long i = n - 1; i >= 0; --i) {
        long u = order[i];
        T acc = e();
        for (const auto& edge : g[u]) {
            long v = edge_to(edge);
            if (v == parent[u]) continue;
            acc = merge(acc, apply_edge(put_edge, dp[v], u, edge));
        }
        dp[u] = put_vertex(acc, u);
    }
    return dp;
}

template <class Edge, class E, class Merge, class PutEdge, class PutVertex>
inline auto rerooting_impl(const std::vector<std::vector<Edge>>& g, E e, Merge merge, PutEdge put_edge, PutVertex put_vertex) {
    using T = std::decay_t<decltype(e())>;
    long n = std::ssize(g);
    if (n == 0) return std::vector<T>{};

    auto [order, parent] = tree_order_impl(g, 0);
    std::vector<T> dp_sub(n, e());
    for (long i = n - 1; i >= 0; --i) {
        long u = order[i];
        T acc = e();
        for (const auto& edge : g[u]) {
            long v = edge_to(edge);
            if (v == parent[u]) continue;
            acc = merge(acc, apply_edge(put_edge, dp_sub[v], u, edge));
        }
        dp_sub[u] = put_vertex(acc, u);
    }

    std::vector<T> from_par(n, e());
    std::vector<T> ans(n, e());
    std::vector<T> contrib, pref, suff;

    for (long u : order) {
        long deg = std::ssize(g[u]);
        contrib.clear();
        contrib.reserve(deg);
        for (const auto& edge : g[u]) {
            long v = edge_to(edge);
            const T& sub_val = (v == parent[u]) ? from_par[u] : dp_sub[v];
            contrib.push_back(apply_edge(put_edge, sub_val, u, edge));
        }

        pref.assign(deg + 1, e());
        suff.assign(deg + 1, e());
        for (long i = 0; i < deg; ++i) {
            pref[i + 1] = merge(pref[i], contrib[i]);
        }
        for (long i = deg - 1; i >= 0; --i) {
            suff[i] = merge(contrib[i], suff[i + 1]);
        }

        ans[u] = put_vertex(pref[deg], u);

        for (long i = 0; i < deg; ++i) {
            long v = edge_to(g[u][i]);
            if (v == parent[u]) continue;
            from_par[v] = put_vertex(merge(pref[i], suff[i + 1]), u);
        }
    }
    return ans;
}

}  // namespace internal

// 根 root からの BFS 順 order と親配列 parent (root の親は -1) を返す
// 計算量: 時間 O(N), 空間 O(N)
// 事前条件: g は N 頂点の連結な無向木（0 <= root < N）
inline std::pair<std::vector<long>, std::vector<long>> tree_order(const std::vector<std::vector<long>>& g, long root = 0) {
    return internal::tree_order_impl(g, root);
}

inline std::pair<std::vector<long>, std::vector<long>> tree_order(const std::vector<std::vector<std::pair<long, long>>>& g, long root = 0) {
    return internal::tree_order_impl(g, root);
}

// 通常のボトムアップ木DP（各頂点 v を根とする部分木の DP 値を返す）
// 計算量: 時間 O(N), 空間 O(N)（各関数オブジェクトの呼び出しを O(1) と仮定）
// 事前条件: g は N 頂点の連結な無向木（0 <= root < N）
template <class E, class Merge, class PutEdge, class PutVertex>
inline auto tree_dp(const std::vector<std::vector<long>>& g, long root, E e, Merge merge, PutEdge put_edge, PutVertex put_vertex) {
    return internal::tree_dp_impl(g, root, e, merge, put_edge, put_vertex);
}

template <class E, class Merge, class PutVertex>
inline auto tree_dp(const std::vector<std::vector<long>>& g, long root, E e, Merge merge, PutVertex put_vertex) {
    using T = std::decay_t<decltype(e())>;
    auto id = [](const T& val, long, long) { return val; };
    return internal::tree_dp_impl(g, root, e, merge, id, put_vertex);
}

template <class E, class Merge, class PutEdge, class PutVertex>
inline auto tree_dp(const std::vector<std::vector<std::pair<long, long>>>& g, long root, E e, Merge merge, PutEdge put_edge,
                    PutVertex put_vertex) {
    return internal::tree_dp_impl(g, root, e, merge, put_edge, put_vertex);
}

// 全方位木DP（各頂点 v を木全体の根としたときの DP 値を返す）
// 計算量: 時間 O(N), 空間 O(N)（各関数オブジェクトの呼び出しを O(1) と仮定）
// 事前条件: g は N 頂点の連結な無向木、かつ (T, merge, e()) は可換モノイドをなすこと
template <class E, class Merge, class PutEdge, class PutVertex>
inline auto rerooting(const std::vector<std::vector<long>>& g, E e, Merge merge, PutEdge put_edge, PutVertex put_vertex) {
    return internal::rerooting_impl(g, e, merge, put_edge, put_vertex);
}

template <class E, class Merge, class PutVertex>
inline auto rerooting(const std::vector<std::vector<long>>& g, E e, Merge merge, PutVertex put_vertex) {
    using T = std::decay_t<decltype(e())>;
    auto id = [](const T& val, long, long) { return val; };
    return internal::rerooting_impl(g, e, merge, id, put_vertex);
}

template <class E, class Merge, class PutEdge, class PutVertex>
inline auto rerooting(const std::vector<std::vector<std::pair<long, long>>>& g, E e, Merge merge, PutEdge put_edge, PutVertex put_vertex) {
    return internal::rerooting_impl(g, e, merge, put_edge, put_vertex);
}
