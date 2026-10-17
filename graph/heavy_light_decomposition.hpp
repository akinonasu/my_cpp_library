#include <algorithm>
#include <cassert>
#include <iostream>
#include <utility>
#include <vector>

class HeavyLightDecomposition {
    int n;
    std::vector<std::vector<int>> g;
    std::vector<int> parent_, depth_, sz_, heavy_, head_, in_, out_, rev_;

    // 1回目のDFS: 部分木サイズ、深さ、親、heavy edgeの決定
    void dfs_sz(int u, int p, int d) {
        parent_[u] = p;
        depth_[u] = d;
        sz_[u] = 1;
        int max_c_sz = 0;
        for (int v : g[u]) {
            if (v == p) continue;
            dfs_sz(v, u, d + 1);
            sz_[u] += sz_[v];
            if (sz_[v] > max_c_sz) {
                max_c_sz = sz_[v];
                heavy_[u] = v;
            }
        }
    }

    // 2回目のDFS: 各パスの先頭(head)とオイラーツアー順(in, out)の決定
    void dfs_hld(int u, int h, int& timer) {
        head_[u] = h;
        in_[u] = timer++;
        rev_[in_[u]] = u;
        if (heavy_[u] != -1) {
            dfs_hld(heavy_[u], h, timer);
        }
        for (int v : g[u]) {
            if (v == parent_[u] || v == heavy_[u]) continue;
            dfs_hld(v, v, timer);
        }
        out_[u] = timer;
    }

   public:
    explicit HeavyLightDecomposition(int n)
        : n(n), g(n), parent_(n, -1), depth_(n, 0), sz_(n, 0), heavy_(n, -1), head_(n, 0), in_(n, 0), out_(n, 0), rev_(n, 0) {}

    // 無向辺を追加
    void add_edge(int u, int v) {
        assert(0 <= u && u < n && 0 <= v && v < n);
        g[u].push_back(v);
        g[v].push_back(u);
    }

    // HLDを構築 (すべての辺を追加し終えた後に必ず1回呼ぶ)
    void build(int root = 0) {
        assert(0 <= root && root < n);
        dfs_sz(root, -1, 0);
        int timer = 0;
        dfs_hld(root, root, timer);
    }

    // -------------------------------------------------------------
    // 各種情報へのアクセサ (ゲッター)
    // -------------------------------------------------------------

    // 頂点 u の親
    [[nodiscard]] int parent(int u) const { return parent_[u]; }

    // 頂点 u の深さ (根は 0)
    [[nodiscard]] int depth(int u) const { return depth_[u]; }

    // 部分木 u のサイズ
    [[nodiscard]] int subtree_size(int u) const { return sz_[u]; }

    // 頂点 u に対応する配列上のインデックス (get_pos と同じ)
    [[nodiscard]] int index(int u) const { return in_[u]; }
    [[nodiscard]] int get_pos(int u) const { return in_[u]; }

    // 部分木 u に対応する半開区間 [subtree_l, subtree_r)
    [[nodiscard]] int subtree_l(int u) const { return in_[u]; }
    [[nodiscard]] int subtree_r(int u) const { return out_[u]; }
    [[nodiscard]] std::pair<int, int> subtree_range(int u) const { return {in_[u], out_[u]}; }

    // 辺 (u, v) のうち、子側（根から遠い方）の頂点を返す
    [[nodiscard]] int edge_to_child(int u, int v) const { return depth_[u] > depth_[v] ? u : v; }

    // -------------------------------------------------------------
    // 木の基本クエリ
    // -------------------------------------------------------------

    // LCA (最小共通祖先)
    [[nodiscard]] int lca(int u, int v) const {
        while (head_[u] != head_[v]) {
            if (depth_[head_[u]] > depth_[head_[v]]) std::swap(u, v);
            v = parent_[head_[v]];
        }
        return depth_[u] < depth_[v] ? u : v;
    }

    // 頂点間の距離 (辺の数)
    [[nodiscard]] int dist(int u, int v) const { return depth_[u] + depth_[v] - 2 * depth_[lca(u, v)]; }

    // 頂点 u から根方向に k 個さかのぼった祖先 (存在しない場合は -1)
    [[nodiscard]] int kth_ancestor(int u, int k) const {
        if (k < 0 || depth_[u] < k) return -1;
        while (u != -1) {
            int d = depth_[u] - depth_[head_[u]];
            if (k <= d) return rev_[in_[u] - k];
            k -= d + 1;
            u = parent_[head_[u]];
        }
        return -1;
    }

    // 頂点 u から v への最短パス上で k 番目の頂点 (0-indexed: k=0 は u)
    [[nodiscard]] int jump(int u, int v, int k) const {
        int l = lca(u, v);
        int d1 = depth_[u] - depth_[l];
        int d2 = depth_[v] - depth_[l];
        if (k < 0 || k > d1 + d2) return -1;
        if (k <= d1) return kth_ancestor(u, k);
        return kth_ancestor(v, d1 + d2 - k);
    }

    // -------------------------------------------------------------
    // セグメント木連携用クエリ
    // -------------------------------------------------------------

    // 可換演算用パスクエリ
    // edge = false : 頂点クエリ
    // edge = true  : 辺クエリ (LCA が除外される)
    template <typename F>
    void path_query(int u, int v, bool edge, F&& f) const {
        while (head_[u] != head_[v]) {
            if (depth_[head_[u]] > depth_[head_[v]]) std::swap(u, v);
            f(in_[head_[v]], in_[v] + 1);
            v = parent_[head_[v]];
        }
        if (depth_[u] > depth_[v]) std::swap(u, v);
        int l = in_[u] + (edge ? 1 : 0);
        int r = in_[v] + 1;
        if (l < r) {
            f(l, r);
        }
    }

    // 非可換演算用パスクエリ (u -> v の向きを考慮)
    // f_up(l, r):   u -> lca 方向 (逆順・右から左に合成)
    // f_down(l, r): lca -> v 方向 (正順・左から右に合成)
    template <typename F1, typename F2>
    void path_non_commutative_query(int u, int v, bool edge, F1&& f_up, F2&& f_down) const {
        int l = lca(u, v);
        int curr_u = u;
        while (head_[curr_u] != head_[l]) {
            f_up(in_[head_[curr_u]], in_[curr_u] + 1);
            curr_u = parent_[head_[curr_u]];
        }
        if (curr_u != l) {
            f_up(in_[l] + 1, in_[curr_u] + 1);
        }
        if (!edge) {
            f_up(in_[l], in_[l] + 1);
        }

        std::vector<std::pair<int, int>> down_intervals;
        int curr_v = v;
        while (head_[curr_v] != head_[l]) {
            down_intervals.emplace_back(in_[head_[curr_v]], in_[curr_v] + 1);
            curr_v = parent_[head_[curr_v]];
        }
        if (curr_v != l) {
            down_intervals.emplace_back(in_[l] + 1, in_[curr_v] + 1);
        }
        for (auto it = down_intervals.rbegin(); it != down_intervals.rend(); ++it) {
            f_down(it->first, it->second);
        }
    }

    // 部分木クエリ
    // edge = false : 頂点クエリ
    // edge = true  : 辺クエリ (根 u に入る辺が除外される)
    template <typename F>
    void subtree_query(int u, bool edge, F&& f) const {
        int l = in_[u] + (edge ? 1 : 0);
        int r = out_[u];
        if (l < r) {
            f(l, r);
        }
    }
};
