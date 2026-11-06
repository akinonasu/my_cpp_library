#include <iostream>
#include <print>
#include <queue>
#include <vector>

// ============================================================================
// ライブラリ: 隣接制約付き最大部分和 (Path Independent Set with Maximum Weight)
// ============================================================================
/**
 * @brief 列から隣接しない要素を k 個選ぶときの最大和を求める
 * @param a 各要素の重み (1-indexedではなく通常の0-indexedのvector)
 * @return std::vector<long> res
 *         res[k] := 隣接制約を満たして k 個選んだときの最大和 (1 <= k <= ceil(N/2))
 *         計算量: O(N log N)
 */
inline std::vector<long> max_nonadjacent_sum(const std::vector<long>& a) {
    int n = a.size();
    int max_k = (n + 1) / 2;
    std::vector<long> res(max_k + 1, 0);
    if (n == 0) return res;

    std::vector<long> val = a;
    std::vector<int> left(n), right(n);
    std::vector<bool> deleted(n, false);

    for (int i = 0; i < n; ++i) {
        left[i] = i - 1;
        right[i] = i + 1;
    }

    // (値, インデックス) の最大ヒープ
    std::priority_queue<std::pair<long, int>> pq;
    for (int i = 0; i < n; ++i) {
        pq.emplace(val[i], i);
    }

    long current_sum = 0;
    for (int k = 1; k <= max_k; ++k) {
        // 削除済み、または値が更新されている古い情報を読み飛ばす
        while (!pq.empty()) {
            auto [v, u] = pq.top();
            if (deleted[u] || v != val[u]) {
                pq.pop();
                continue;
            }
            break;
        }

        auto [v, u] = pq.top();
        pq.pop();

        current_sum += v;
        res[k] = current_sum;

        int l = left[u];
        int r = right[u];

        if (l != -1 && r != n) {
            // 両側に隣接要素がある場合: 両隣を削除し、差分値で置き換えてヒープに再追加
            deleted[l] = true;
            deleted[r] = true;
            val[u] = val[l] + val[r] - val[u];

            left[u] = left[l];
            right[u] = right[r];
            if (left[u] != -1) right[left[u]] = u;
            if (right[u] != n) left[right[u]] = u;

            pq.emplace(val[u], u);
        } else if (l != -1 && r == n) {
            // 右端の場合: 左隣のみ削除され、反転要素は作られない
            deleted[u] = true;
            deleted[l] = true;
            if (left[l] != -1) right[left[l]] = n;
        } else if (l == -1 && r != n) {
            // 左端の場合: 右隣のみ削除され、反転要素は作られない
            deleted[u] = true;
            deleted[r] = true;
            if (right[r] != n) left[right[r]] = -1;
        } else {
            // 孤立した1要素の場合
            deleted[u] = true;
        }
    }

    return res;
}
