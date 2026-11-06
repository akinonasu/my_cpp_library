#include <algorithm>
#include <print>
#include <vector>

/**
 * @brief 尺取り法（条件を満たす最大の右端を探索）
 * @details 各 l に対し、条件を満たす最大の半開区間 [l, r) を探索します。
 *          「区間を広げると条件に違反しやすくなる」問題（和が K 以下、種類数が K 以下など）に適しています。
 *
 * @param n 要素数
 * @param can_add (int r) -> bool : r を追加しても条件を満たすか判定
 * @param add (int r) -> void : r を区間に追加
 * @param remove (int l) -> void : l を区間から削除
 * @param process (int l, int r) -> void : 各 l に対する最大の r が決定したときのコールバック
 */
template <typename CanAdd, typename Add, typename Remove, typename Process>
void two_pointers(int n, CanAdd can_add, Add add, Remove remove, Process process) {
    int r = 0;
    for (int l = 0; l < n; ++l) {
        while (r < n && can_add(r)) {
            add(r);
            ++r;
        }
        process(l, r);
        if (l == r) {
            ++r;
        } else {
            remove(l);
        }
    }
}

/**
 * @brief 条件を満たす部分配列の総数をカウント
 * @return long 条件を満たす部分配列の個数 (最大 n*(n+1)/2 となり 32bit を超えるため long)
 */
template <typename CanAdd, typename Add, typename Remove>
long count_subarrays(int n, CanAdd can_add, Add add, Remove remove) {
    long count = 0;
    two_pointers(n, can_add, add, remove, [&](int l, int r) { count += (r - l); });
    return count;
}

/**
 * @brief 条件を満たす最長部分配列の長さを取得
 * @return int 最長長（条件を満たす部分配列が 1 つも存在しない場合は 0）
 */
template <typename CanAdd, typename Add, typename Remove>
int max_subarray_length(int n, CanAdd can_add, Add add, Remove remove) {
    int max_len = 0;
    two_pointers(n, can_add, add, remove, [&](int l, int r) { max_len = std::max(max_len, r - l); });
    return max_len;
}

/**
 * @brief 尺取り法（条件を満たす間、左端を縮める）
 * @details r を 1 つずつ進め、条件を満たしている間 l を縮めます。
 *          「区間を広げると条件を満たしやすくなる」問題（和が K 以上、全要素を含む最短区間など）に適しています。
 *
 * @param n 要素数
 * @param add (int r) -> void : r を区間に追加
 * @param is_valid () -> bool : 現在の区間 [l, r + 1) が条件を満たしているか判定
 * @param remove (int l) -> void : l を区間から削除
 * @param process (int l, int r) -> void : 条件を満たす極小な区間 [l, r) が見つかるたびに呼ばれるコールバック
 */
template <typename Add, typename IsValid, typename Remove, typename Process>
void two_pointers_shrink(int n, Add add, IsValid is_valid, Remove remove, Process process) {
    int l = 0;
    for (int r = 0; r < n; ++r) {
        add(r);
        while (l <= r && is_valid()) {
            process(l, r + 1);
            remove(l);
            ++l;
        }
    }
}

/**
 * @brief 条件を満たす最短部分配列の長さを取得
 * @return int 最短長（条件を満たす区間が存在しない場合は -1）
 */
template <typename Add, typename IsValid, typename Remove>
int min_subarray_length(int n, Add add, IsValid is_valid, Remove remove) {
    int min_len = n + 1;
    two_pointers_shrink(n, add, is_valid, remove, [&](int l, int r) { min_len = std::min(min_len, r - l); });
    return min_len > n ? -1 : min_len;
}

/**
 * @brief 固定長スライディングウィンドウ
 * @param n 要素数
 * @param k ウィンドウの幅 (1 <= k <= n)
 * @param add (int idx) -> void : 要素を追加
 * @param remove (int idx) -> void : 要素を削除
 * @param process (int l, int r) -> void : 幅 k の区間 [l, r) に対する処理
 */
template <typename Add, typename Remove, typename Process>
void sliding_window_fixed(int n, int k, Add add, Remove remove, Process process) {
    if (k <= 0 || k > n) return;
    for (int i = 0; i < k; ++i) {
        add(i);
    }
    process(0, k);
    for (int i = k; i < n; ++i) {
        add(i);
        remove(i - k);
        process(i - k + 1, i + 1);
    }
}
