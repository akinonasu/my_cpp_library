#include <algorithm>
#include <cassert>
#include <print>
#include <tuple>
#include <utility>
#include <vector>

class Tree {
   public:
    struct PathSegment {
        int l;     // 半開区間 [l, r) の左端 (0 <= l < r <= n)
        int r;     // 半開区間 [l, r) の右端
        bool rev;  // true なら r-1 から l の向き（u->LCA の上り方向）
    };

    Tree() = default;

    // parent[i]: 頂点 i の親 (-1 は根)
    // self_is_root == true の場合、parent[i] == i も根として扱う
    explicit Tree(const std::vector<int>& parent, bool self_is_root = true) : n_(static_cast<int>(parent.size())), par_(n_, -1), children_(n_), dep_(n_, 0), sz_(n_, 1), root_(n_, -1), head_(n_, -1), tin_(n_, -1), tout_(n_, -1), inv_(n_, -1) {
        if (n_ == 0) return;

        for (int i = 0; i < n_; ++i) {
            if (parent[i] == -1 || (self_is_root && parent[i] == i)) {
                par_[i] = -1;
                roots_.push_back(i);
            } else {
                assert(0 <= parent[i] && parent[i] < n_);
                par_[i] = parent[i];
                children_[parent[i]].push_back(i);
            }
        }

        // BFS でトポロジカル順（親が先）と深さ・根を計算
        std::vector<int> order;
        order.reserve(n_);
        for (int r : roots_) {
            root_[r] = r;
            dep_[r] = 0;
            order.push_back(r);
        }
        for (int idx = 0; idx < static_cast<int>(order.size()); ++idx) {
            int v = order[idx];
            for (int c : children_[v]) {
                root_[c] = root_[v];
                dep_[c] = dep_[v] + 1;
                order.push_back(c);
            }
        }
        // サイクルが存在する場合、到達できない頂点が残りここで検出される
        assert(static_cast<int>(order.size()) == n_);

        // 逆順で部分木サイズを計算し、最大の重さの子 (heavy child) を先頭に移動
        for (int idx = n_ - 1; idx >= 0; --idx) {
            int v = order[idx];
            for (int& c : children_[v]) {
                sz_[v] += sz_[c];
                if (sz_[c] > sz_[children_[v][0]]) {
                    std::swap(c, children_[v][0]);
                }
            }
        }

        // 非再帰 DFS で HLD (head) と Euler Tour (tin, tout, inv) を構築
        int timer = 0;
        std::vector<int> st;
        st.reserve(n_);
        for (int r : roots_) {
            head_[r] = r;
            st.push_back(r);
            while (!st.empty()) {
                int v = st.back();
                st.pop_back();
                if (v >= 0) {
                    tin_[v] = timer++;
                    inv_[tin_[v]] = v;
                    st.push_back(~v);
                    for (int i = static_cast<int>(children_[v].size()) - 1; i >= 0; --i) {
                        int c = children_[v][i];
                        head_[c] = (i == 0 ? head_[v] : c);
                        st.push_back(c);
                    }
                } else {
                    tout_[~v] = timer;
                }
            }
        }
    }

    int size() const { return n_; }
    const std::vector<int>& roots() const { return roots_; }

    int root(int v) const {
        assert(0 <= v && v < n_);
        return root_[v];
    }

    int parent(int v) const {
        assert(0 <= v && v < n_);
        return par_[v];
    }

    const std::vector<int>& children(int v) const {
        assert(0 <= v && v < n_);
        return children_[v];
    }

    int depth(int v) const {
        assert(0 <= v && v < n_);
        return dep_[v];
    }

    int subtree_size(int v) const {
        assert(0 <= v && v < n_);
        return sz_[v];
    }

    // u が v の祖先か (u == v も true)
    bool is_ancestor(int u, int v) const {
        assert(0 <= u && u < n_);
        assert(0 <= v && v < n_);
        return tin_[u] <= tin_[v] && tout_[v] <= tout_[u];
    }

    // 最小共通祖先 (別の木なら -1)
    int lca(int u, int v) const {
        assert(0 <= u && u < n_);
        assert(0 <= v && v < n_);
        if (root_[u] != root_[v]) return -1;
        while (head_[u] != head_[v]) {
            if (dep_[head_[u]] < dep_[head_[v]]) std::swap(u, v);
            u = par_[head_[u]];
        }
        return dep_[u] < dep_[v] ? u : v;
    }

    // 2頂点間の距離 (別の木なら -1)
    int dist(int u, int v) const {
        int w = lca(u, v);
        if (w == -1) return -1;
        return dep_[u] + dep_[v] - 2 * dep_[w];
    }

    // v から k 個親を辿った祖先 (k == 0 なら v, 存在しなければ -1)
    int kth_ancestor(int v, int k) const {
        assert(0 <= v && v < n_);
        assert(k >= 0);
        if (k > dep_[v]) return -1;
        int target_dep = dep_[v] - k;
        while (dep_[head_[v]] > target_dep) {
            v = par_[head_[v]];
        }
        return inv_[tin_[v] - (dep_[v] - target_dep)];
    }

    // u から v へ k 歩進んだ頂点 (到達不可・別の木なら -1)
    int jump(int u, int v, int k) const {
        assert(0 <= u && u < n_);
        assert(0 <= v && v < n_);
        assert(k >= 0);
        int w = lca(u, v);
        if (w == -1) return -1;
        int d1 = dep_[u] - dep_[w];
        int d2 = dep_[v] - dep_[w];
        if (k > d1 + d2) return -1;
        if (k <= d1) return kth_ancestor(u, k);
        return kth_ancestor(v, d1 + d2 - k);
    }

    // u から v への最短パス上の頂点リスト (別の木なら空)
    std::vector<int> path(int u, int v) const {
        auto segs = path_segments(u, v, false);
        std::vector<int> res;
        for (const auto& [l, r, rev] : segs) {
            if (!rev) {
                for (int t = l; t < r; ++t)
                    res.push_back(inv_[t]);
            } else {
                for (int t = r - 1; t >= l; --t)
                    res.push_back(inv_[t]);
            }
        }
        return res;
    }

    // --- セグ木・BIT連携用 (Euler Tour / HLD) ---
    int tin(int v) const {
        assert(0 <= v && v < n_);
        return tin_[v];
    }

    int tout(int v) const {
        assert(0 <= v && v < n_);
        return tout_[v];
    }

    int vertex_at(int t) const {
        assert(0 <= t && t < n_);
        return inv_[t];
    }

    // 部分木 v に対応する半開区間 [l, r)
    std::pair<int, int> subtree_range(int v) const {
        assert(0 <= v && v < n_);
        return {tin_[v], tout_[v]};
    }

    // u -> v パス上の HLD 区間列 (edge == true なら辺属性＝LCA を除く)
    // 非可換モノイドでもそのまま畳み込めるよう順序と rev フラグを保持して返す
    std::vector<PathSegment> path_segments(int u, int v, bool edge = false) const {
        assert(0 <= u && u < n_);
        assert(0 <= v && v < n_);
        if (root_[u] != root_[v]) return {};

        std::vector<PathSegment> up, down;
        while (head_[u] != head_[v]) {
            if (dep_[head_[u]] > dep_[head_[v]]) {
                up.push_back({tin_[head_[u]], tin_[u] + 1, true});
                u = par_[head_[u]];
            } else {
                down.push_back({tin_[head_[v]], tin_[v] + 1, false});
                v = par_[head_[v]];
            }
        }
        if (dep_[u] > dep_[v]) {
            int l = tin_[v] + (edge ? 1 : 0);
            int r = tin_[u] + 1;
            if (l < r) up.push_back({l, r, true});
        } else {
            int l = tin_[u] + (edge ? 1 : 0);
            int r = tin_[v] + 1;
            if (l < r) down.push_back({l, r, false});
        }
        std::ranges::reverse(down);
        up.insert(up.end(), down.begin(), down.end());
        return up;
    }

    // v が属する木の直径: {直径の長さ, 端点1, 端点2}
    std::tuple<int, int, int> diameter(int v = 0) const {
        assert(0 <= v && v < n_);
        auto farthest = [&](int start) -> std::pair<int, int> {
            std::vector<std::tuple<int, int, int>> st;
            st.emplace_back(start, -1, 0);
            int best_d = -1, best_v = start;
            while (!st.empty()) {
                auto [cur, prev, d] = st.back();
                st.pop_back();
                if (d > best_d) {
                    best_d = d;
                    best_v = cur;
                }
                if (par_[cur] != -1 && par_[cur] != prev) {
                    st.emplace_back(par_[cur], cur, d + 1);
                }
                for (int c : children_[cur]) {
                    if (c != prev) {
                        st.emplace_back(c, cur, d + 1);
                    }
                }
            }
            return {best_d, best_v};
        };

        auto [_, u] = farthest(v);
        auto [len, w] = farthest(u);
        return {len, u, w};
    }

    // v が属する木の重心 (1個または2個)
    std::vector<int> centroids(int v = 0) const {
        assert(0 <= v && v < n_);
        int c = root_[v];
        int total = sz_[c];
        while (!children_[c].empty() && sz_[children_[c][0]] * 2 > total) {
            c = children_[c][0];
        }
        std::vector<int> res = {c};
        if (!children_[c].empty() && sz_[children_[c][0]] * 2 == total) {
            res.push_back(children_[c][0]);
        }
        return res;
    }

   private:
    int n_ = 0;
    std::vector<int> roots_;
    std::vector<int> par_;
    std::vector<std::vector<int>> children_;
    std::vector<int> dep_;
    std::vector<int> sz_;
    std::vector<int> root_;
    std::vector<int> head_;
    std::vector<int> tin_;
    std::vector<int> tout_;
    std::vector<int> inv_;
};
