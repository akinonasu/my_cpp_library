#include <algorithm>
#include <cassert>
#include <iostream>
#include <ranges>
#include <vector>

class StronglyConnectedComponents {
   private:
    int n_;
    std::vector<std::vector<int>> g_;      // 順方向
    std::vector<std::vector<int>> rev_g_;  // 逆方向

    std::vector<int> comp_id_;              // 各頂点が属する成分ID (0-indexed, トポロジカル順)
    std::vector<std::vector<int>> groups_;  // 成分ごとの頂点リスト
    bool is_built_ = false;

    void dfs(int v, std::vector<bool>& visited, std::vector<int>& order) {
        visited[v] = true;
        for (int to : g_[v]) {
            if (!visited[to]) dfs(to, visited, order);
        }
        order.push_back(v);
    }

    void rdfs(int v, int cur_id, std::vector<int>& group) {
        comp_id_[v] = cur_id;
        group.push_back(v);
        for (int to : rev_g_[v]) {
            if (comp_id_[to] == -1) rdfs(to, cur_id, group);
        }
    }

    void build() {
        if (is_built_) return;
        is_built_ = true;

        std::vector<bool> visited(n_, false);
        std::vector<int> order;
        order.reserve(n_);

        for (int i = 0; i < n_; ++i) {
            if (!visited[i]) dfs(i, visited, order);
        }

        comp_id_.assign(n_, -1);
        int cur_id = 0;

        // 帰りがけ順の逆から探索することで、トポロジカル順に成分が構築される
        for (int v : order | std::views::reverse) {
            if (comp_id_[v] == -1) {
                groups_.emplace_back();
                rdfs(v, cur_id, groups_.back());
                cur_id++;
            }
        }
    }

   public:
    explicit StronglyConnectedComponents(int n) : n_(n), g_(n), rev_g_(n) {}

    // 有向辺 (from -> to) を追加
    void add_edge(int from, int to) {
        assert(0 <= from && from < n_);
        assert(0 <= to && to < n_);
        g_[from].push_back(to);
        rev_g_[to].push_back(from);
        is_built_ = false;  // 辺が追加されたらキャッシュを無効化
    }

    // 強連結成分に分解した頂点群のリストを返す (トポロジカル順)
    // 戻り値: std::vector<std::vector<int>>
    const std::vector<std::vector<int>>& scc() {
        build();
        return groups_;
    }

    // 各頂点がどの強連結成分に属するか（ID）を返す
    // 戻り値: 各頂点の成分ID (0 <= id < 成分数)
    const std::vector<int>& ids() {
        build();
        return comp_id_;
    }

    // 縮約後のDAG（有向非巡回グラフ）を構築して返す
    // 重複した辺は自動的に除去されます
    std::vector<std::vector<int>> build_dag() {
        build();
        int num_groups = static_cast<int>(groups_.size());
        std::vector<std::vector<int>> dag(num_groups);

        for (int from = 0; from < n_; ++from) {
            for (int to : g_[from]) {
                int u = comp_id_[from];
                int v = comp_id_[to];
                if (u != v) {
                    dag[u].push_back(v);
                }
            }
        }

        // 多重辺の削除 (C++20/23 ranges)
        for (auto& adj : dag) {
            std::ranges::sort(adj);
            auto ret = std::ranges::unique(adj);
            adj.erase(ret.begin(), ret.end());
        }

        return dag;
    }

    // 成分の総数を返す
    int size() {
        build();
        return static_cast<int>(groups_.size());
    }

    // 特定の頂点が属する成分IDを返す
    int operator[](int v) {
        build();
        assert(0 <= v && v < n_);
        return comp_id_[v];
    }
};
