#include <algorithm>
#include <iostream>
#include <numeric>
#include <print>
#include <ranges>
#include <vector>

// --- Binary Indexed Tree (Fenwick Tree) ---
template <typename T>
struct FenwickTree {
    int n;
    std::vector<T> data;

    FenwickTree(int n) : n(n), data(n + 1, 0) {}

    // i番目 (0-indexed) に v を加算
    void add(int i, T v) {
        for (++i; i <= n; i += i & -i) {
            data[i] += v;
        }
    }

    // [0, i) の総和を取得
    T sum(int i) const {
        T s = 0;
        for (; i > 0; i -= i & -i) {
            s += data[i];
        }
        return s;
    }
};

/**
 * @brief 転倒数 (Inversion Number) を求める（隣接要素のswapによる最小ソート回数）
 * @note 最大値が N*(N-1)/2 に達し 32bit を超えるため戻り値は long 型
 * @complexity O(N log N)
 */
template <typename T>
long inversion_number(const std::vector<T>& a) {
    int n = a.size();
    if (n <= 1) return 0;

    // 座標圧縮
    std::vector<T> vals = a;
    std::ranges::sort(vals);
    auto ret = std::ranges::unique(vals);
    vals.erase(ret.begin(), ret.end());

    FenwickTree<int> bit(vals.size());
    long inv = 0;

    for (int i = 0; i < n; ++i) {
        // a[i] 以上の値が既に何個出現したかを数える
        int rank = std::ranges::lower_bound(vals, a[i]) - vals.begin();
        inv += i - bit.sum(rank + 1);
        bit.add(rank, 1);
    }
    return inv;
}

/**
 * @brief 任意の2要素のswapによる最小置換回数を求める
 * @note 各連結成分（サイクル）のサイズを L としたとき、L - 1 回のswapで整列可能。
 *       全体での最小回数は N - (サイクルの個数) となる。
 *       ※重複要素がある場合、同値の要素は元の相対順序を保つ（安定ソート）として計算します。
 * @complexity O(N log N)
 */
template <typename T>
int min_swaps_to_sort(const std::vector<T>& a) {
    int n = a.size();
    if (n <= 1) return 0;

    // (値, 元のインデックス)
    std::vector<std::pair<T, int>> p(n);
    for (int i = 0; i < n; ++i) {
        p[i] = {a[i], i};
    }
    // 安定ソート（値が同じなら元のインデックス順）
    std::ranges::stable_sort(p, [](const auto& left, const auto& right) { return left.first < right.first; });

    std::vector<bool> visited(n, false);
    int cycles = 0;

    for (int i = 0; i < n; ++i) {
        if (visited[i]) continue;

        // 新しいサイクルを検出
        cycles++;
        int cur = i;
        while (!visited[cur]) {
            visited[cur] = true;
            cur = p[cur].second;  // ソート前の位置へ移動
        }
    }

    return n - cycles;
}
