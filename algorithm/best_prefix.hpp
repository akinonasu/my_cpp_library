#include <algorithm>
#include <cassert>
#include <iostream>
#include <print>
#include <span>
#include <vector>

// IsMax = true  : 最大値を求める (BestMaxPrefix)
// IsMax = false : 最小値を求める (BestMinPrefix)
template <bool IsMax = true>
class BestPrefix {
   public:
    static constexpr long INF = 1L << 60;  // 64bit整数の安全な無限大
    static constexpr long UNUSED = IsMax ? -INF : INF;

   private:
    int n = 0;
    std::vector<long> pref;  // pref[i] = i まで進んだときの最適値

    static constexpr long choose_better(long a, long b) noexcept {
        if constexpr (IsMax)
            return std::max(a, b);
        else
            return std::min(a, b);
    }

   public:
    BestPrefix() = default;

    explicit BestPrefix(int n_) { reset(n_); }

    explicit BestPrefix(std::span<const long> init, bool force_one = false) {
        reset(static_cast<int>(init.size()));
        if (force_one) {
            merge_with_one(init);
        } else {
            merge_with_empty(init);
        }
    }

    void reset(int n_) {
        n = n_;
        pref.assign(n + 1, UNUSED);
        pref[0] = 0L;
    }

    // 配列 v から 0 個以上の要素を使う（空区間を許容）
    void merge_with_empty(std::span<const long> v) {
        assert(static_cast<int>(v.size()) == n);
        for (int i = 0; i < n; ++i) {
            if (pref[i] == UNUSED) continue;
            pref[i + 1] = choose_better(pref[i + 1], pref[i] + v[i]);
        }
    }

    // 配列 v から 少なくとも1個以上の要素を使う（空区間を許容しない）
    void merge_with_one(std::span<const long> v) {
        assert(static_cast<int>(v.size()) == n);
        std::vector<long> next(n + 1, UNUSED);

        for (int i = 0; i < n; ++i) {
            // ここで v の採用を開始する
            if (pref[i] != UNUSED) {
                next[i + 1] = choose_better(next[i + 1], pref[i] + v[i]);
            }
            // 既に開始している v の採用を継続する
            if (next[i] != UNUSED) {
                next[i + 1] = choose_better(next[i + 1], next[i] + v[i]);
            }
        }
        pref = std::move(next);
    }

    // 便宜用：値がすべて 0 のゼロ配列で空遷移（スキップ）を行う
    void skip_zero() {
        for (int i = 0; i < n; ++i) {
            if (pref[i] == UNUSED) continue;
            pref[i + 1] = choose_better(pref[i + 1], pref[i]);
        }
    }

    // 途中のインデックス i までの最適値を取得
    long get(int i) const {
        assert(0 <= i && i <= n);
        return pref[i];
    }

    // 全体（インデックス n）まで割り当て終えたときの最適値
    long get_best() const { return pref[n]; }
};

// 型エイリアス
using BestMaxPrefix = BestPrefix<true>;
using BestMinPrefix = BestPrefix<false>;
