#include <iostream>
#include <vector>
#include <queue>
#include <print>
#include <bit>
#include <algorithm>
#include <cassert>

struct FunctionalGraph {
    int n;
    std::vector<int> to;
    std::vector<int> in_deg;

    // 分解後の情報
    // cycle_id: 所属する連結成分のサイクル番号
    // pos_in_cycle: サイクル内の位置 (0-indexed, サイクル外は -1)
    // dist_to_cycle: サイクルまでの距離 (サイクル上の頂点は 0)
    // root: その頂点が最初に到達するサイクル上の頂点
    std::vector<int> cycle_id;
    std::vector<int> pos_in_cycle;
    std::vector<int> dist_to_cycle;
    std::vector<int> root;
    std::vector<std::vector<int>> cycles;

    // ダブリングテーブル (木部分のジャンプ用)
    int log_steps = 0;
    std::vector<std::vector<int>> up;

    FunctionalGraph(int n) 
        : n(n), to(n), in_deg(n, 0), cycle_id(n, -1),
          pos_in_cycle(n, -1), dist_to_cycle(n, 0), root(n, -1) {}

    // 有向辺 u -> v を追加
    void add_edge(int u, int v) {
        to[u] = v;
        in_deg[v]++;
    }

    // グラフの構築と成分分解
    void build() {
        // 1. 入次数を使ったトポロジカルソートで木部分を取り除く
        std::vector<int> deg = in_deg;
        std::queue<int> q;
        for (int i = 0; i < n; ++i) {
            if (deg[i] == 0) q.push(i);
        }
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            int v = to[u];
            if (--deg[v] == 0) {
                q.push(v);
            }
        }

        // 2. サイクルを検出
        std::vector<bool> visited(n, false);
        for (int i = 0; i < n; ++i) {
            if (deg[i] > 0 && !visited[i]) {
                std::vector<int> cycle;
                int curr = i;
                while (!visited[curr]) {
                    visited[curr] = true;
                    cycle.push_back(curr);
                    curr = to[curr];
                }
                int cid = static_cast<int>(cycles.size());
                for (int j = 0; j < static_cast<int>(cycle.size()); ++j) {
                    int u = cycle[j];
                    cycle_id[u] = cid;
                    pos_in_cycle[u] = j;
                    dist_to_cycle[u] = 0;
                    root[u] = u;
                }
                cycles.push_back(std::move(cycle));
            }
        }

        // 3. サイクルから逆方向に BFS して木部分の距離・根を決定
        std::vector<std::vector<int>> rev_to(n);
        for (int i = 0; i < n; ++i) {
            rev_to[to[i]].push_back(i);
        }

        for (const auto& cycle : cycles) {
            for (int u : cycle) {
                q.push(u);
            }
        }

        int max_dist = 0;
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (int p : rev_to[u]) {
                if (cycle_id[p] == -1) { // サイクル外の頂点
                    cycle_id[p] = cycle_id[u];
                    root[p] = root[u];
                    dist_to_cycle[p] = dist_to_cycle[u] + 1;
                    max_dist = std::max(max_dist, dist_to_cycle[p]);
                    q.push(p);
                }
            }
        }

        // 4. 木部分のジャンプに必要な分だけダブリングを構築
        if (max_dist > 0) {
            log_steps = std::bit_width(static_cast<unsigned int>(max_dist));
            up.assign(log_steps, std::vector<int>(n, -1));
            for (int i = 0; i < n; ++i) {
                up[0][i] = to[i];
            }
            for (int b = 1; b < log_steps; ++b) {
                for (int i = 0; i < n; ++i) {
                    up[b][i] = up[b - 1][up[b - 1][i]];
                }
            }
        }
    }

    // 頂点 u がサイクル上に存在するか判定
    bool in_cycle(int u) const {
        return dist_to_cycle[u] == 0;
    }

    // 頂点 u が所属するサイクルの長さを取得
    long cycle_length(int u) const {
        return static_cast<long>(cycles[cycle_id[u]].size());
    }

    // 頂点 u から k ステップ進んだ頂点を取得 (k >= 0)
    int jump(int u, long k) const {
        assert(k >= 0);
        // サイクルに到達する場合: O(1)
        if (k >= dist_to_cycle[u]) {
            k -= dist_to_cycle[u];
            int r = root[u];
            int cid = cycle_id[r];
            long len = static_cast<long>(cycles[cid].size());
            long p = (pos_in_cycle[r] + (k % len)) % len;
            return cycles[cid][p];
        }
        // 木の途中に留まる場合: ダブリングで O(log N)
        int curr = u;
        for (int b = 0; b < log_steps; ++b) {
            if ((k >> b) & 1) {
                curr = up[b][curr];
            }
        }
        return curr;
    }

    // 頂点 u から頂点 v へ到達可能なら最小ステップ数、不可能な場合 -1 を返す
    long dist(int u, int v) const {
        if (cycle_id[u] != cycle_id[v]) return -1;

        if (in_cycle(v)) {
            // v がサイクル上にある場合
            long d1 = dist_to_cycle[u];
            int r = root[u];
            long len = cycle_length(r);
            long d2 = (pos_in_cycle[v] - pos_in_cycle[r] + len) % len;
            return d1 + d2;
        } else {
            // v が木上にある場合: u は同じ部分木上にあり、深さが v 以上である必要がある
            if (root[u] != root[v]) return -1;
            long d = dist_to_cycle[u] - dist_to_cycle[v];
            if (d < 0) return -1;
            if (jump(u, d) == v) return d;
            return -1;
        }
    }
};
