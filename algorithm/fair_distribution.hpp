#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

// 均等分配 (Round-Robin 分配)
// 構築・各種クエリ: 時間 O(1), 空間 O(1)
// to_vector: 時間 O(n), 空間 O(n)
// 事前条件: m >= 0, n > 0 (すべて 0-indexed, 半開区間 [l, r_idx))
struct FairDistribution {
    long m = 0;  // 配る物の総数
    long n = 1;  // 配る対象の人数
    long q = 0;  // 全員に配られる基本の個数 (m / n)
    long r = 0;  // 余り (m % n): 先頭 r 人が q + 1 個もらう

    FairDistribution() = default;

    FairDistribution(long total_items, long num_recipients) : m(total_items), n(num_recipients) {
        assert(m >= 0);
        assert(n > 0);
        q = m / n;
        r = m % n;
    }

    // i 番目の人がもらう個数
    long get(long i) const {
        assert(0 <= i && i < n);
        return q + (i < r ? 1 : 0);
    }

    long operator[](long i) const { return get(i); }

    // 先頭 k 人がもらう合計個数
    long prefix_sum(long k) const {
        assert(0 <= k && k <= n);
        return k * q + std::min(k, r);
    }

    // 区間 [l, r_idx) の人がもらう合計個数
    long range_sum(long l, long r_idx) const {
        assert(0 <= l && l <= r_idx && r_idx <= n);
        return prefix_sum(r_idx) - prefix_sum(l);
    }

    // 全員の個数を std::vector として取得
    template <class T = long>
    std::vector<T> to_vector() const {
        std::vector<T> res(n, static_cast<T>(q));
        std::fill_n(res.begin(), r, static_cast<T>(q + 1));
        return res;
    }
};
