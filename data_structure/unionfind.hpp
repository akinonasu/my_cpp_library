#include <numeric>
#include <print>
#include <utility>
#include <vector>

// 素集合データ構造 (Union-Find)
struct UnionFind {
    int n;
    // parent_or_size[i] < 0 のとき、i は根であり、その値の絶対値は連結成分のサイズを表す
    // parent_or_size[i] >= 0 のとき、parent_or_size[i] は親の頂点番号を表す
    std::vector<int> parent_or_size;

    explicit UnionFind(int n = 0) : n(n), parent_or_size(n, -1) {}

    // x が属する連結成分の代表元を返す（経路圧縮あり）
    int root(int x) {
        if (parent_or_size[x] < 0) return x;
        return parent_or_size[x] = root(parent_or_size[x]);
    }

    // x と y を同じ連結成分に結合する
    // 既に同じ集合なら false、結合されたなら true を返す
    bool unite(int x, int y) {
        x = root(x);
        y = root(y);
        if (x == y) return false;

        // サイズが大きい方に小さい方をマージ（Union by Size）
        if (-parent_or_size[x] < -parent_or_size[y]) {
            std::swap(x, y);
        }
        parent_or_size[x] += parent_or_size[y];
        parent_or_size[y] = x;
        return true;
    }

    // x と y が同じ連結成分に属するか判定
    bool same(int x, int y) { return root(x) == root(y); }

    // x が属する連結成分のサイズを返す
    int size(int x) { return -parent_or_size[root(x)]; }

    // 連結成分ごとに頂点のリストを分類して返す
    std::vector<std::vector<int>> groups() {
        std::vector<int> leader_buf(n);
        std::vector<int> group_size(n, 0);
        for (int i = 0; i < n; ++i) {
            leader_buf[i] = root(i);
            group_size[leader_buf[i]]++;
        }

        std::vector<std::vector<int>> result(n);
        for (int i = 0; i < n; ++i) {
            result[i].reserve(group_size[i]);
        }
        for (int i = 0; i < n; ++i) {
            result[leader_buf[i]].push_back(i);
        }

        std::vector<std::vector<int>> res;
        res.reserve(n);
        for (int i = 0; i < n; ++i) {
            if (!result[i].empty()) {
                res.push_back(std::move(result[i]));
            }
        }
        return res;
    }
};
