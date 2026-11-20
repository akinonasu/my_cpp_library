#pragma once

#include <algorithm>
#include <cassert>
#include <limits>
#include <utility>
#include <vector>

// ==========================================
// 接頭辞最適値 DP (Best Prefix)
// IsMax = true  : 最大値を求める (BestMaxPrefix)
// IsMax = false : 最小値を求める (BestMinPrefix)
// 時間計算量: 構築 / reset / merge_* / skip_zero O(N), get / get_best O(1)
// 空間計算量: O(N)
// 事前条件:
//   - 0 <= n
//   - merge_* に渡す配列 v のサイズは n と等しいこと
// ==========================================
template <bool IsMax = true>
struct BestPrefix {
    static constexpr long inf = std::numeric_limits<long>::max() / 4;
    static constexpr long unused = IsMax ? -inf : inf;

    BestPrefix() : BestPrefix(0) {}

    explicit BestPrefix(long n) { reset(n); }

    template <class T>
    explicit BestPrefix(const std::vector<T>& init, bool force_one = false) {
        reset(std::ssize(init));
        if (force_one) {
            merge_with_one(init);
        } else {
            merge_with_empty(init);
        }
    }

    void reset(long n) {
        assert(0 <= n);
        n_ = n;
        pref_.assign(n_ + 1, unused);
        pref_[0] = 0;
    }

    // 配列 v から 0 個以上の要素を使う（空区間を許容）
    template <class T>
    void merge_with_empty(const std::vector<T>& v) {
        assert(std::ssize(v) == n_);
        for (long i = 0; i < n_; ++i) {
            if (pref_[i] == unused) continue;
            pref_[i + 1] = choose_better(pref_[i + 1], pref_[i] + static_cast<long>(v[i]));
        }
    }

    // 配列 v から 1 個以上の要素を使う（空区間を許容しない）
    template <class T>
    void merge_with_one(const std::vector<T>& v) {
        assert(std::ssize(v) == n_);
        std::vector<long> next(n_ + 1, unused);
        for (long i = 0; i < n_; ++i) {
            long val = static_cast<long>(v[i]);
            if (pref_[i] != unused) {
                next[i + 1] = choose_better(next[i + 1], pref_[i] + val);
            }
            if (next[i] != unused) {
                next[i + 1] = choose_better(next[i + 1], next[i] + val);
            }
        }
        pref_ = std::move(next);
    }

    // 値がすべて 0 の配列で空遷移（スキップ）を行う
    void skip_zero() {
        for (long i = 0; i < n_; ++i) {
            if (pref_[i] == unused) continue;
            pref_[i + 1] = choose_better(pref_[i + 1], pref_[i]);
        }
    }

    // 接頭辞 [0, i) までの最適値を取得
    long get(long i) const {
        assert(0 <= i && i <= n_);
        return pref_[i];
    }

    // 全体 [0, n) までの最適値を取得
    long get_best() const { return pref_[n_]; }

   private:
    long n_ = 0;
    std::vector<long> pref_;

    static long choose_better(long a, long b) {
        if constexpr (IsMax) {
            return std::max(a, b);
        } else {
            return std::min(a, b);
        }
    }
};

using BestMaxPrefix = BestPrefix<true>;
using BestMinPrefix = BestPrefix<false>;