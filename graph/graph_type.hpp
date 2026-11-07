#include <algorithm>
#include <iostream>
#include <numeric>
#include <print>
#include <queue>
#include <vector>

// 無向グラフの形状判定クラス
class UndirectedGraphClassifier {
    int n;
    int m = 0;
    std::vector<std::vector<int>> adj;
    std::vector<int> deg;
    bool simple = true;
    bool self_loops = false;
    bool multi_edges = false;

   public:
    explicit UndirectedGraphClassifier(int n) : n(n), adj(n), deg(n, 0) {}

    // 0-indexedで辺を追加
    void add_edge(int u, int v) {
        m++;
        if (u == v) self_loops = true;
        adj[u].push_back(v);
        adj[v].push_back(u);
        deg[u]++;
        deg[v]++;
    }

    // グラフの初期化・多重辺のチェック（辺を追加し終えた後に呼ぶ）
    void build() {
        for (int i = 0; i < n; ++i) {
            std::ranges::sort(adj[i]);
            for (int j = 1; j < (int)adj[i].size(); ++j) {
                if (adj[i][j] == adj[i][j - 1]) {
                    multi_edges = true;
                }
            }
        }
        simple = (!self_loops && !multi_edges);
    }

    int vertex_count() const { return n; }
    int edge_count() const { return m; }
    bool is_simple() const { return simple; }

    // 連結成分数を数える
    int count_connected_components() const {
        std::vector<bool> visited(n, false);
        int components = 0;
        for (int i = 0; i < n; ++i) {
            if (!visited[i]) {
                components++;
                std::queue<int> q;
                q.push(i);
                visited[i] = true;
                while (!q.empty()) {
                    int u = q.front();
                    q.pop();
                    for (int v : adj[u]) {
                        if (!visited[v]) {
                            visited[v] = true;
                            q.push(v);
                        }
                    }
                }
            }
        }
        return components;
    }

    // 連結グラフか判定
    bool is_connected() const {
        if (n == 0) return true;
        return count_connected_components() == 1;
    }

    // 二部グラフ判定
    bool is_bipartite() const {
        std::vector<int> color(n, -1);
        for (int i = 0; i < n; ++i) {
            if (color[i] != -1) continue;
            std::queue<int> q;
            q.push(i);
            color[i] = 0;
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                for (int v : adj[u]) {
                    if (color[v] == -1) {
                        color[v] = 1 - color[u];
                        q.push(v);
                    } else if (color[v] == color[u]) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    // 木 (Tree): 連結かつ辺数が N - 1 かつ単純グラフ
    bool is_tree() const {
        if (!simple) return false;
        return m == n - 1 && is_connected();
    }

    // 森 (Forest): サイクルが存在しない単純グラフ
    bool is_forest() const {
        if (!simple) return false;
        return m == n - count_connected_components();
    }

    // パスグラフ (Path): 1本の一連のパス
    bool is_path() const {
        if (!is_tree()) return false;
        if (n <= 1) return true;
        int deg1 = 0, deg2 = 0;
        for (int d : deg) {
            if (d == 1)
                deg1++;
            else if (d == 2)
                deg2++;
            else
                return false;
        }
        return (deg1 == 2 && deg2 == n - 2);
    }

    // サイクルグラフ (Cycle): 単純閉路 (N >= 3)
    bool is_cycle() const {
        if (!simple || n < 3) return false;
        if (m != n || !is_connected()) return false;
        for (int d : deg) {
            if (d != 2) return false;
        }
        return true;
    }

    // スターグラフ (Star): 1つの頂点に他のすべての頂点が接続 (N <= 3 も網羅)
    bool is_star() const {
        if (!is_tree()) return false;
        if (n <= 2) return true;
        int center_count = 0;
        int leaf_count = 0;
        for (int d : deg) {
            if (d == n - 1)
                center_count++;
            else if (d == 1)
                leaf_count++;
        }
        return (center_count == 1 && leaf_count == n - 1);
    }

    // 完全グラフ (Complete Graph): すべての頂点対間に辺が存在する
    bool is_complete() const {
        if (!simple) return false;
        long max_edges = (long)n * (n - 1) / 2;
        return m == max_edges;
    }

    // なもりグラフ (Pseudo-tree / Namori graph): 連結かつ M = N (サイクルをちょうど1つ含む)
    bool is_namori() const {
        if (!simple) return false;
        return m == n && is_connected();
    }

    // 正則グラフ (Regular Graph): 全ての頂点の次数が k
    bool is_regular(int k = -1) const {
        if (n == 0) return true;
        int target = (k == -1) ? deg[0] : k;
        for (int d : deg) {
            if (d != target) return false;
        }
        return true;
    }

    // オイラーグラフ判定 (Eulerian circuit を持つか: 全ての次数が偶数かつ孤立点を除いて連結)
    bool is_eulerian() const {
        for (int d : deg) {
            if (d % 2 != 0) return false;
        }
        // 辺が存在する頂点同士が同一連結成分にあるか判定
        int non_zero_v = -1;
        for (int i = 0; i < n; ++i) {
            if (deg[i] > 0) {
                non_zero_v = i;
                break;
            }
        }
        if (non_zero_v == -1) return true;  // 辺がない場合は真

        std::vector<bool> visited(n, false);
        std::queue<int> q;
        q.push(non_zero_v);
        visited[non_zero_v] = true;
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (int v : adj[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    q.push(v);
                }
            }
        }
        for (int i = 0; i < n; ++i) {
            if (deg[i] > 0 && !visited[i]) return false;
        }
        return true;
    }
};
