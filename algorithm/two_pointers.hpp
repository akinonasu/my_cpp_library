#include <algorithm>
#include <cassert>
#include <concepts>
#include <print>
#include <type_traits>
#include <vector>

// ヘルパー: F が (l, r) を受け取るか、(index) のみを受け取るかを判定して呼び出す
template <class F>
auto invoke_with_indices(F& f, int l, int r, int target_idx) {
    if constexpr (std::invocable<F, int, int>) {
        return f(l, r);
    } else {
        return f(target_idx);
    }
}

template <class F>
auto invoke_condition(F& f, int l, int r) {
    if constexpr (std::invocable<F, int, int>) {
        return f(l, r);
    } else {
        return f();
    }
}

/**
 * @brief パターン1: 条件を満たす「最大」の区間 [l, r) を探索する
 *
 * 性質: 区間が条件を満たすとき、その部分区間も必ず条件を満たす（単調性）
 * 例: 区間の総和 <= K, 異なり数 <= K
 *
 * @param n 配列の長さ
 * @param can_push_right 右端 r を追加可能か判定する関数: bool(int r) または bool(int l, int r)
 * @param push_right 右端 r を追加する関数: void(int r) または void(int l, int r)
 * @param pop_left 左端 l を削除する関数: void(int l) または void(int l, int r)
 * @param report 各 l に対する最大の r を受け取る関数: void(int l, int r) または bool(int l, int r)
 *               (boolを返し、falseを返した場合は早期終了)
 */
template <class CanPush, class PushRight, class PopLeft, class Report>
void shakutori_max(int n, CanPush can_push_right, PushRight push_right, PopLeft pop_left, Report report) {
    assert(n >= 0);
    int r = 0;
    for (int l = 0; l < n; ++l) {
        while (r < n && invoke_with_indices(can_push_right, l, r, r)) {
            invoke_with_indices(push_right, l, r, r);
            ++r;
        }

        assert(l <= r);
        if constexpr (std::same_as<std::invoke_result_t<Report, int, int>, bool>) {
            if (!report(l, r)) break;
        } else {
            report(l, r);
        }

        if (l < r) {
            invoke_with_indices(pop_left, l, r, l);
        } else {
            // 区間が空 [l, l) の場合、要素 l は入っていないため pop_left は呼ばず
            // 次の l+1 に向けて r を揃える
            r = l + 1;
        }
    }
}

/**
 * @brief パターン2: 条件を満たす「最小」の区間 [l, r) を探索する
 *
 * 性質: 区間が条件を満たすとき、それを包含する超区間も必ず条件を満たす（単調性）
 * 例: 区間の総和 >= K, 異なり数 >= K
 *
 * @param n 配列の長さ
 * @param is_valid 現在の区間 [l, r) が条件を満たしているか判定する関数: bool() または bool(int l, int r)
 * @param push_right 右端 r を追加する関数: void(int r) または void(int l, int r)
 * @param pop_left 左端 l を削除する関数: void(int l) または void(int l, int r)
 * @param report 各 l に対する最小の r を受け取る関数: void(int l, int r) または bool(int l, int r)
 *               (満たす区間が存在しない場合は r = -1 が渡される)
 */
template <class IsValid, class PushRight, class PopLeft, class Report>
void shakutori_min(int n, IsValid is_valid, PushRight push_right, PopLeft pop_left, Report report) {
    assert(n >= 0);
    int r = 0;
    for (int l = 0; l < n; ++l) {
        while (r < n && !invoke_condition(is_valid, l, r)) {
            invoke_with_indices(push_right, l, r, r);
            ++r;
        }

        assert(l <= r);
        int valid_r = invoke_condition(is_valid, l, r) ? r : -1;

        if constexpr (std::same_as<std::invoke_result_t<Report, int, int>, bool>) {
            if (!report(l, valid_r)) break;
        } else {
            report(l, valid_r);
        }

        if (l < r) {
            invoke_with_indices(pop_left, l, r, l);
        } else {
            r = l + 1;
        }
    }
}

// -------------------------------------------------------------
// 頻出クエリ用の便利関数
// -------------------------------------------------------------

// 条件を満たす連続部分列の総数をカウント（総和 <= K など）
template <class CanPush, class PushRight, class PopLeft>
long count_subarrays_max(int n, CanPush can_push_right, PushRight push_right, PopLeft pop_left) {
    long count = 0;
    shakutori_max(n, can_push_right, push_right, pop_left, [&](int l, int r) { count += (r - l); });
    return count;
}

// 条件を満たす最長部分列の長さを取得（存在しない場合は 0）
template <class CanPush, class PushRight, class PopLeft>
int max_subarray_length(int n, CanPush can_push_right, PushRight push_right, PopLeft pop_left) {
    int max_len = 0;
    shakutori_max(n, can_push_right, push_right, pop_left, [&](int l, int r) { max_len = std::max(max_len, r - l); });
    return max_len;
}

// 条件を満たす連続部分列の総数をカウント（総和 >= K など）
template <class IsValid, class PushRight, class PopLeft>
long count_subarrays_min(int n, IsValid is_valid, PushRight push_right, PopLeft pop_left) {
    long count = 0;
    shakutori_min(n, is_valid, push_right, pop_left, [&](int l, int r) {
        if (r != -1) {
            count += (n - r + 1);
        }
    });
    return count;
}

// 条件を満たす最短部分列の長さを取得（存在しない場合は -1）
template <class IsValid, class PushRight, class PopLeft>
int min_subarray_length(int n, IsValid is_valid, PushRight push_right, PopLeft pop_left) {
    int min_len = -1;
    shakutori_min(n, is_valid, push_right, pop_left, [&](int l, int r) {
        if (r != -1) {
            int len = r - l;
            if (min_len == -1 || len < min_len) {
                min_len = len;
            }
        }
    });
    return min_len;
}

// -------------------------------------------------------------
// 手動で細かく制御したい場合のステートクラス
// -------------------------------------------------------------
template <class PushRight, class PopLeft>
class SlidingWindow {
    int n_;
    int l_ = 0;
    int r_ = 0;
    PushRight push_right_;
    PopLeft pop_left_;

   public:
    SlidingWindow(int n, PushRight push_right, PopLeft pop_left) : n_(n), push_right_(std::move(push_right)), pop_left_(std::move(pop_left)) { assert(n >= 0); }

    int left() const { return l_; }
    int right() const { return r_; }
    int length() const { return r_ - l_; }
    bool empty() const { return l_ == r_; }

    bool can_advance_right() const { return r_ < n_; }
    bool can_advance_left() const { return l_ < r_; }

    void advance_right() {
        assert(r_ < n_);
        invoke_with_indices(push_right_, l_, r_, r_);
        ++r_;
    }

    void advance_left() {
        assert(l_ < r_);
        invoke_with_indices(pop_left_, l_, r_, l_);
        ++l_;
    }
};

/*
int main() {
    // ---------------------------------------------------------
    // 例1: 出す/入れるで操作が異なる場合 (異なり数 <= K の最長区間)
    // ---------------------------------------------------------
    // 配列 a の中で、含まれる要素の種類数が K 個以下である最長部分列の長さを求める
    std::vector<int> a = {1, 2, 1, 3, 4, 3, 5};
    int n = static_cast<int>(a.size());
    int k = 2;

    std::vector<int> freq(10, 0);
    int distinct_count = 0;

    // 入れる操作: 頻度が 0 -> 1 になったら種類数 +1
    auto push_r = [&](int r) {
        if (freq[a[r]]++ == 0) {
            distinct_count++;
        }
    };

    // 出す操作: 頻度が 1 -> 0 になったら種類数 -1
    auto pop_l = [&](int l) {
        if (--freq[a[l]] == 0) {
            distinct_count--;
        }
    };

    // 入れても条件を満たすか？
    auto can_push = [&](int r) {
        return (freq[a[r]] > 0 || distinct_count + 1 <= k);
    };

    int max_len = shakutori::max_subarray_length(n, can_push, push_r, pop_l);
    std::println("異なり数 <= {} の最長連続部分列長: {}", k, max_len);

    // ---------------------------------------------------------
    // 例2: 総和が K 以上の最短部分列長 (存在しない場合は -1)
    // ---------------------------------------------------------
    std::vector<int> b = {2, 3, 1, 2, 4, 3};
    int n2 = static_cast<int>(b.size());
    long target_sum = 7;
    long current_sum = 0;

    auto push_sum = [&](int r) { current_sum += b[r]; };
    auto pop_sum = [&](int l) { current_sum -= b[l]; };
    auto is_valid = [&]() { return current_sum >= target_sum; };

    int min_len = shakutori::min_subarray_length(n2, is_valid, push_sum, pop_sum);
    std::println("総和 >= {} の最短連続部分列長: {}", target_sum, min_len);

    // ---------------------------------------------------------
    // 例3: 総和が K 以下の連続部分列の個数 (long型でオーバーフロー防止)
    // ---------------------------------------------------------
    std::vector<int> c = {10, 5, 2, 6};
    int n3 = static_cast<int>(c.size());
    long limit = 10;
    long win_sum = 0;

    auto can_push_sum = [&](int r) { return win_sum + c[r] <= limit; };
    auto push_c = [&](int r) { win_sum += c[r]; };
    auto pop_c = [&](int l) { win_sum -= c[l]; };

    long count = shakutori::count_subarrays_max(n3, can_push_sum, push_c, pop_c);
    std::println("総和 <= {} の連続部分列の総数: {}", limit, count);

    return 0;
}
*/
