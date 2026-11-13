#include <algorithm>
#include <cassert>
#include <print>
#include <utility>
#include <vector>

struct CycleResult {
    bool ok;
    std::vector<int> vertices;  // 閉路の頂点列 [v_0, v_1, ..., v_{k-1}]
    std::vector<int> edges;     // 閉路の辺ID列 [e_0, e_1, ..., e_{k-1}]
};

class CycleDetector {
    int n;
    bool directed;
    struct Edge {
        int to;
        int id;
    };
    std::vector<std::vector<Edge>> adj;
    int edge_cnt = 0;

    // 無向グラフの単純奇閉路 (BFS: O(N + M))
    CycleResult find_undirected_odd_cycle() const {
        std::vector<int> color(n, -1);
        std::vector<int> depth(n, 0);
        std::vector<int> parent_v(n, -1);
        std::vector<int> parent_e(n, -1);
        std::vector<int> q;

        for (int s = 0; s < n; ++s) {
            if (color[s] != -1) continue;
            color[s] = 0;
            q.push_back(s);

            for (int head = 0; head < static_cast<int>(q.size()); ++head) {
                int u = q[head];
                for (auto [v, id] : adj[u]) {
                    if (u == v) {
                        return {true, {u}, {id}};
                    }
                    if (color[v] == -1) {
                        color[v] = color[u] ^ 1;
                        depth[v] = depth[u] + 1;
                        parent_v[v] = u;
                        parent_e[v] = id;
                        q.push_back(v);
                    } else if (color[v] == color[u]) {
                        // 奇閉路発見: u と v から LCA まで遡る
                        std::vector<int> vs_u, es_u, vs_v, es_v;
                        int cu = u, cv = v;
                        while (depth[cu] > depth[cv]) {
                            vs_u.push_back(cu);
                            es_u.push_back(parent_e[cu]);
                            cu = parent_v[cu];
                        }
                        while (depth[cv] > depth[cu]) {
                            vs_v.push_back(cv);
                            es_v.push_back(parent_e[cv]);
                            cv = parent_v[cv];
                        }
                        while (cu != cv) {
                            vs_u.push_back(cu);
                            es_u.push_back(parent_e[cu]);
                            cu = parent_v[cu];
                            vs_v.push_back(cv);
                            es_v.push_back(parent_e[cv]);
                            cv = parent_v[cv];
                        }
                        int lca = cu;
                        std::ranges::reverse(vs_u);
                        std::ranges::reverse(es_u);

                        std::vector<int> vs = {lca};
                        vs.insert(vs.end(), vs_u.begin(), vs_u.end());
                        vs.insert(vs.end(), vs_v.begin(), vs_v.end());

                        std::vector<int> es = std::move(es_u);
                        es.push_back(id);
                        es.insert(es.end(), es_v.begin(), es_v.end());
                        return {true, std::move(vs), std::move(es)};
                    }
                }
            }
        }
        return {false, {}, {}};
    }

    // 有向グラフの単純奇閉路 (Kosaraju SCC + 倍長グラフBFS + スタック縮小: O(N + M))
    CycleResult find_directed_odd_cycle() const {
        // 自己ループがあれば長さ1の奇閉路
        for (int u = 0; u < n; ++u) {
            for (auto [v, id] : adj[u]) {
                if (u == v) return {true, {u}, {id}};
            }
        }

        // 1. Kosaraju SCC (非再帰)
        std::vector<bool> vis(n, false);
        std::vector<int> iter(n, 0);
        std::vector<int> order, st;
        order.reserve(n);

        for (int s = 0; s < n; ++s) {
            if (vis[s]) continue;
            vis[s] = true;
            st.push_back(s);
            while (!st.empty()) {
                int u = st.back();
                if (iter[u] < static_cast<int>(adj[u].size())) {
                    int v = adj[u][iter[u]++].to;
                    if (!vis[v]) {
                        vis[v] = true;
                        st.push_back(v);
                    }
                } else {
                    st.pop_back();
                    order.push_back(u);
                }
            }
        }

        std::vector<std::vector<int>> radj(n);
        for (int u = 0; u < n; ++u) {
            for (auto [v, id] : adj[u])
                radj[v].push_back(u);
        }

        std::vector<int> comp(n, -1);
        std::vector<std::vector<int>> scc_nodes;
        std::vector<int> q;
        for (int i = n - 1; i >= 0; --i) {
            int s = order[i];
            if (comp[s] != -1) continue;
            int cid = static_cast<int>(scc_nodes.size());
            scc_nodes.emplace_back();
            comp[s] = cid;
            q.clear();
            q.push_back(s);
            for (int head = 0; head < static_cast<int>(q.size()); ++head) {
                int u = q[head];
                scc_nodes[cid].push_back(u);
                for (int v : radj[u]) {
                    if (comp[v] == -1) {
                        comp[v] = cid;
                        q.push_back(v);
                    }
                }
            }
        }

        // 2. 各SCC内で倍長グラフBFSを行い、奇閉路があれば単純閉路へ縮小する
        std::vector<int> dist(2 * n, -1);
        std::vector<int> prev_v(2 * n, -1);
        std::vector<int> prev_e(2 * n, -1);
        std::vector<int> pos(n, -1);

        for (int cid = 0; cid < static_cast<int>(scc_nodes.size()); ++cid) {
            if (scc_nodes[cid].size() <= 1) continue;
            int s = scc_nodes[cid][0];
            q.clear();
            dist[2 * s] = 0;
            q.push_back(2 * s);

            bool found = false;
            for (int head = 0; head < static_cast<int>(q.size()); ++head) {
                int cur = q[head];
                int u = cur / 2;
                int p = cur & 1;
                for (auto [v, id] : adj[u]) {
                    if (comp[v] != cid) continue;
                    int nxt = 2 * v + (p ^ 1);
                    if (dist[nxt] == -1) {
                        dist[nxt] = dist[cur] + 1;
                        prev_v[nxt] = cur;
                        prev_e[nxt] = id;
                        q.push_back(nxt);
                        if (nxt == 2 * s + 1) {
                            found = true;
                            break;
                        }
                    }
                }
                if (found) break;
            }

            if (found) {
                // (s, 0) -> (s, 1) の奇数長ウォークを復元
                std::vector<int> w_vs, w_es;
                int cur = 2 * s + 1;
                while (cur != 2 * s) {
                    w_vs.push_back(cur / 2);
                    w_es.push_back(prev_e[cur]);
                    cur = prev_v[cur];
                }
                w_vs.push_back(s);
                std::ranges::reverse(w_vs);
                std::ranges::reverse(w_es);

                // スタックで偶閉路を取り除き、単純奇閉路を抽出する
                std::vector<int> st_v, st_e;
                st_v.push_back(w_vs[0]);
                pos[w_vs[0]] = 0;
                int len_w = static_cast<int>(w_es.size());

                for (int i = 1; i <= len_w; ++i) {
                    int v = w_vs[i];
                    int e = w_es[i - 1];
                    st_e.push_back(e);
                    if (pos[v] != -1) {
                        int idx = pos[v];
                        int clen = static_cast<int>(st_v.size()) - idx;
                        if (clen % 2 == 1) {
                            std::vector<int> vs(st_v.begin() + idx, st_v.end());
                            std::vector<int> es(st_e.begin() + idx, st_e.end());
                            return {true, std::move(vs), std::move(es)};
                        }
                        // 偶数長の閉路は捨てて残りの奇数長ウォークを継続
                        while (static_cast<int>(st_v.size()) > idx + 1) {
                            pos[st_v.back()] = -1;
                            st_v.pop_back();
                        }
                        st_e.resize(idx);
                    } else {
                        pos[v] = static_cast<int>(st_v.size());
                        st_v.push_back(v);
                    }
                }
            }

            for (int u : scc_nodes[cid]) {
                dist[2 * u] = -1;
                dist[2 * u + 1] = -1;
            }
        }
        return {false, {}, {}};
    }

    // 無向グラフの単純偶閉路 (DFS木 + 木辺を共有する2本の後退辺からの3パス構築: O(N + M))
    CycleResult find_undirected_even_cycle() const {
        std::vector<bool> vis(n, false);
        std::vector<int> depth(n, 0);
        std::vector<int> parent_v(n, -1);
        std::vector<int> parent_e(n, -1);
        std::vector<int> iter(n, 0);
        std::vector<int> st;

        for (int s = 0; s < n; ++s) {
            if (vis[s]) continue;
            vis[s] = true;
            st.push_back(s);
            while (!st.empty()) {
                int u = st.back();
                if (iter[u] < static_cast<int>(adj[u].size())) {
                    auto [v, id] = adj[u][iter[u]++];
                    if (!vis[v]) {
                        vis[v] = true;
                        depth[v] = depth[u] + 1;
                        parent_v[v] = u;
                        parent_e[v] = id;
                        st.push_back(v);
                    }
                } else {
                    st.pop_back();
                }
            }
        }

        struct BackEdge {
            int u, v, id;
        };
        std::vector<BackEdge> cover(n, {-1, -1, -1});

        for (int u = 0; u < n; ++u) {
            for (auto [v, id] : adj[u]) {
                // 非木辺（自己ループ除く）を u が子孫、v が先祖の向きで1回ずつ処理
                if (id == parent_e[u] || depth[u] <= depth[v]) continue;

                // 1. 基本閉路が偶数長ならそのまま復元
                if ((depth[u] - depth[v] + 1) % 2 == 0) {
                    std::vector<int> vs, es;
                    int curr = u;
                    while (curr != v) {
                        vs.push_back(curr);
                        es.push_back(parent_e[curr]);
                        curr = parent_v[curr];
                    }
                    vs.push_back(v);
                    std::ranges::reverse(vs);
                    std::ranges::reverse(es);
                    es.push_back(id);
                    return {true, std::move(vs), std::move(es)};
                }

                // 2. 木辺を辿り、同じ木辺をカバーする別の後退辺があれば3本の素なパスから構築
                int curr = u;
                while (curr != v) {
                    if (cover[curr].id == -1) {
                        cover[curr] = {u, v, id};
                        curr = parent_v[curr];
                    } else {
                        auto [u1, v1, id1] = cover[curr];
                        auto [u2, v2, id2] = BackEdge{u, v, id};
                        if (depth[v1] > depth[v2]) {
                            std::swap(u1, u2);
                            std::swap(v1, v2);
                            std::swap(id1, id2);
                        }

                        // w = LCA(u1, u2)
                        int cu = u1, cv = u2;
                        while (depth[cu] > depth[cv])
                            cu = parent_v[cu];
                        while (depth[cv] > depth[cu])
                            cv = parent_v[cv];
                        while (cu != cv) {
                            cu = parent_v[cu];
                            cv = parent_v[cv];
                        }
                        int w = cu;

                        struct Path {
                            std::vector<int> vs, es;
                        };
                        auto climb = [&](int from, int to) -> Path {
                            Path p;
                            int x = from;
                            while (x != to) {
                                p.vs.push_back(x);
                                p.es.push_back(parent_e[x]);
                                x = parent_v[x];
                            }
                            p.vs.push_back(to);
                            return p;
                        };

                        // P1: w -> v2 (木上を登る)
                        Path p1 = climb(w, v2);

                        // P2: w -> u1 ->(id1)-> v1 -> v2
                        Path p2;
                        {
                            Path down1 = climb(u1, w);
                            std::ranges::reverse(down1.vs);
                            std::ranges::reverse(down1.es);
                            Path down2 = climb(v2, v1);
                            std::ranges::reverse(down2.vs);
                            std::ranges::reverse(down2.es);

                            p2.vs = std::move(down1.vs);
                            p2.es = std::move(down1.es);
                            p2.es.push_back(id1);
                            p2.vs.insert(p2.vs.end(), down2.vs.begin(), down2.vs.end());
                            p2.es.insert(p2.es.end(), down2.es.begin(), down2.es.end());
                        }

                        // P3: w -> u2 ->(id2)-> v2
                        Path p3;
                        {
                            Path down = climb(u2, w);
                            std::ranges::reverse(down.vs);
                            std::ranges::reverse(down.es);
                            p3.vs = std::move(down.vs);
                            p3.es = std::move(down.es);
                            p3.es.push_back(id2);
                            p3.vs.push_back(v2);
                        }

                        auto merge_even = [&](Path& a, Path& b) -> CycleResult {
                            std::vector<int> vs = std::move(a.vs);
                            std::vector<int> es = std::move(a.es);
                            std::ranges::reverse(b.vs);
                            std::ranges::reverse(b.es);
                            // b.vs は [v2, ..., w] なので先頭(v2)と末尾(w)を除いて繋ぐ
                            vs.insert(vs.end(), b.vs.begin() + 1, b.vs.end() - 1);
                            es.insert(es.end(), b.es.begin(), b.es.end());
                            return {true, std::move(vs), std::move(es)};
                        };

                        if (p1.es.size() % 2 == p2.es.size() % 2) return merge_even(p1, p2);
                        if (p1.es.size() % 2 == p3.es.size() % 2) return merge_even(p1, p3);
                        return merge_even(p2, p3);
                    }
                }
            }
        }
        return {false, {}, {}};
    }

   public:
    explicit CycleDetector(int n, bool directed = true) : n(n), directed(directed), adj(n) { assert(n >= 0); }

    int add_edge(int u, int v) {
        assert(0 <= u && u < n);
        assert(0 <= v && v < n);
        int id = edge_cnt++;
        adj[u].push_back({v, id});
        if (!directed && u != v) {
            adj[v].push_back({u, id});
        }
        return id;
    }

    // 任意の単純閉路を検出・復元する
    CycleResult find_cycle() const {
        std::vector<int> state(n, 0);
        std::vector<int> parent_v(n, -1);
        std::vector<int> parent_e(n, -1);
        std::vector<int> iter(n, 0);
        std::vector<int> st;

        for (int s = 0; s < n; ++s) {
            if (state[s] != 0) continue;
            state[s] = 1;
            st.push_back(s);

            while (!st.empty()) {
                int u = st.back();
                if (iter[u] < static_cast<int>(adj[u].size())) {
                    auto [v, id] = adj[u][iter[u]++];
                    if (!directed && id == parent_e[u]) continue;

                    if (state[v] == 0) {
                        state[v] = 1;
                        parent_v[v] = u;
                        parent_e[v] = id;
                        st.push_back(v);
                    } else if (state[v] == 1) {
                        std::vector<int> vs, es;
                        int curr = u;
                        while (curr != v) {
                            vs.push_back(curr);
                            es.push_back(parent_e[curr]);
                            curr = parent_v[curr];
                        }
                        vs.push_back(v);
                        std::ranges::reverse(vs);
                        std::ranges::reverse(es);
                        es.push_back(id);
                        return {true, std::move(vs), std::move(es)};
                    }
                } else {
                    state[u] = 2;
                    st.pop_back();
                }
            }
        }
        return {false, {}, {}};
    }

    // 単純奇閉路を検出・復元する (有向・無向両対応)
    CycleResult find_odd_cycle() const {
        if (directed) {
            return find_directed_odd_cycle();
        } else {
            return find_undirected_odd_cycle();
        }
    }

    // 単純偶閉路を検出・復元する (無向グラフ専用)
    CycleResult find_even_cycle() const {
        assert(!directed);
        return find_undirected_even_cycle();
    }
};
