#include <algorithm>
#include <concepts>
#include <iostream>
#include <limits>
#include <print>
#include <vector>

/**
 * @brief 動的更新対応 ワーシャルフロイド (Floyd-Warshall)
 * @tparam T コストの型 (デフォルト: long long)
 */
template <typename T = long long>
    requires std::is_arithmetic_v<T>
class DynamicFloydWarshall {
   public:
    static constexpr T INF = std::numeric_limits<T>::max() / 3;

   private:
    int n;
    std::vector<std::vector<T>> d;
    std::vector<std::vector<int>> nxt;

    // 有向辺 (u -> v) の更新を O(V^2) で行う内部関数
    bool update_directed(int u, int v, T cost) {
        // 既存の距離より短くならなければ改善の見込みなし
        if (cost >= d[u][v]) return false;

        d[u][v] = cost;
        nxt[u][v] = v;
        bool updated = false;

        // 全てのペア (i, j) に対し、i -> u -> v -> j のパスで短縮できるか検証
        for (int i = 0; i < n; ++i) {
            if (d[i][u] >= INF) continue;
            for (int j = 0; j < n; ++j) {
                if (d[v][j] >= INF) continue;

                T new_dist = d[i][u] + cost + d[v][j];
                if (new_dist < d[i][j]) {
                    d[i][j] = new_dist;
                    nxt[i][j] = (i == u ? v : nxt[i][u]);
                    updated = true;
                }
            }
        }
        return updated;
    }

   public:
    explicit DynamicFloydWarshall(int n) : n(n), d(n, std::vector<T>(n, INF)), nxt(n, std::vector<int>(n, -1)) {
        for (int i = 0; i < n; ++i) {
            d[i][i] = 0;
            nxt[i][i] = i;
        }
    }

    /**
     * @brief 初期グラフへの辺の追加 (build() の前に呼ぶ)
     */
    void add_edge(int from, int to, T cost, bool directed = false) {
        if (cost < d[from][to]) {
            d[from][to] = cost;
            nxt[from][to] = to;
        }
        if (!directed) {
            if (cost < d[to][from]) {
                d[to][from] = cost;
                nxt[to][from] = from;
            }
        }
    }

    /**
     * @brief 全点対間最短経路の初期構築 O(V^3)
     */
    void build() {
        for (int k = 0; k < n; ++k) {
            for (int i = 0; i < n; ++i) {
                if (d[i][k] >= INF) continue;
                for (int j = 0; j < n; ++j) {
                    if (d[k][j] >= INF) continue;
                    if (d[i][k] + d[k][j] < d[i][j]) {
                        d[i][j] = d[i][k] + d[k][j];
                        nxt[i][j] = nxt[i][k];
                    }
                }
            }
        }
    }

    /**
     * @brief 辺の追加 / コスト短縮による動的更新 O(V^2)
     * @return 最短距離の更新が発生したかどうか
     */
    bool update(int from, int to, T cost, bool directed = false) {
        if (directed) {
            return update_directed(from, to, cost);
        } else {
            bool r1 = update_directed(from, to, cost);
            bool r2 = update_directed(to, from, cost);
            return r1 || r2;
        }
    }

    /**
     * @brief C++23 多次元添字演算子: fw[u, v] で距離取得 O(1)
     */
    [[nodiscard]] constexpr T operator[](int from, int to) const noexcept { return d[from][to]; }

    [[nodiscard]] constexpr T dist(int from, int to) const noexcept { return d[from][to]; }

    /**
     * @brief 到達可能か判定
     */
    [[nodiscard]] bool is_reachable(int from, int to) const noexcept { return d[from][to] < INF / 2; }

    /**
     * @brief 負閉路が存在するか判定
     */
    [[nodiscard]] bool has_negative_cycle() const noexcept {
        for (int i = 0; i < n; ++i) {
            if (d[i][i] < 0) return true;
        }
        return false;
    }

    /**
     * @brief 最短経路の頂点列を取得 (到達不能時は空配列)
     */
    [[nodiscard]] std::vector<int> get_path(int from, int to) const {
        if (!is_reachable(from, to)) return {};
        std::vector<int> path;
        int curr = from;
        path.push_back(curr);
        while (curr != to) {
            curr = nxt[curr][to];
            if (curr == -1) return {};
            path.push_back(curr);
            if (static_cast<int>(path.size()) > n) {
                // 負閉路による無限ループの検出
                return {};
            }
        }
        return path;
    }
};
