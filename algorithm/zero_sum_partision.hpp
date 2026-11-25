#pragma once

#include <algorithm>
#include <bit>
#include <cassert>
#include <functional>
#include <utility>
#include <vector>

// 部分集合ごとの演算結果の全列挙
// 計算量: 時間 O(2^N), 空間 O(2^N)
// 事前条件: 0 <= N <= 25, op は結合的・可換で単位元 e を持つこと
template <class T, class Op = std::plus<>>
std::vector<T> subset_sums(const std::vector<T>& a, T e = T{}, Op op = {}) {
    long n = std::ssize(a);
    assert(0 <= n && n <= 25);
    std::vector<T> sum(1L << n, e);
    for (long mask = 1; mask < (1L << n); ++mask) {
        long i = std::countr_zero(static_cast<unsigned long>(mask));
        sum[mask] = op(sum[mask ^ (1L << i)], a[i]);
    }
    return sum;
}

namespace internal {

template <class T>
std::vector<long> build_zero_sum_dp(const std::vector<T>& sum, const T& e, long n) {
    long all = (1L << n) - 1;
    std::vector<long> dp(1L << n, 0);
    for (long mask = 0; mask < all; ++mask) {
        long cur = dp[mask];
        for (long rem = all ^ mask; rem > 0; rem &= rem - 1) {
            long bit = rem & -rem;
            long nmask = mask | bit;
            long val = cur + (sum[nmask] == e ? 1 : 0);
            if (val > dp[nmask]) {
                dp[nmask] = val;
            }
        }
    }
    return dp;
}

}  // namespace internal

// 可換群における零和部分集合への最大分割数（個数のみ）
// 計算量: 時間 O(N 2^N), 空間 O(2^N)
// 事前条件: 0 <= N <= 25, (T, op, e) は可換群をなすこと
template <class T, class Op = std::plus<>>
long zero_sum_partition_count(const std::vector<T>& a, T e = T{}, Op op = {}) {
    long n = std::ssize(a);
    auto sum = subset_sums(a, e, op);
    auto dp = internal::build_zero_sum_dp(sum, e, n);
    return dp.back();
}

// 可換群における零和部分集合への最大分割（0-indexed の添字リストとして復元）
// 計算量: 時間 O(N 2^N), 空間 O(2^N)
// 事前条件: zero_sum_partition_count と同じ
template <class T, class Op = std::plus<>>
std::vector<std::vector<long>> zero_sum_partition(const std::vector<T>& a, T e = T{}, Op op = {}) {
    long n = std::ssize(a);
    auto sum = subset_sums(a, e, op);
    auto dp = internal::build_zero_sum_dp(sum, e, n);

    std::vector<std::vector<long>> groups;
    long mask = (1L << n) - 1;
    while (mask > 0) {
        long is_zero = (sum[mask] == e ? 1 : 0);
        if (is_zero) {
            groups.emplace_back();
        }
        long pmask = -1;
        long picked = -1;
        for (long rem = mask; rem > 0; rem &= rem - 1) {
            long bit = rem & -rem;
            long cand = mask ^ bit;
            if (dp[mask] == dp[cand] + is_zero) {
                pmask = cand;
                picked = std::countr_zero(static_cast<unsigned long>(bit));
                break;
            }
        }
        assert(pmask != -1);
        if (!groups.empty()) {
            groups.back().push_back(picked);
        }
        mask = pmask;
    }
    for (auto& g : groups) {
        std::reverse(g.begin(), g.end());
    }
    std::reverse(groups.begin(), groups.end());
    return groups;
}

struct Transfer {
    long from;
    long to;
    long amount;
};

// 総和が等しい供給リスト {id, amount} と需要リスト {id, amount} を貪欲にマッチングする
// 計算量: 時間 O(|supply| + |demand|), 空間 O(|supply| + |demand|)
// 事前条件: 各 amount > 0、かつ supply の総和 == demand の総和
inline std::vector<Transfer> greedy_transport(const std::vector<std::pair<long, long>>& supply,
                                              const std::vector<std::pair<long, long>>& demand) {
    long ns = std::ssize(supply);
    long nd = std::ssize(demand);
    std::vector<Transfer> res;
    if (ns == 0 && nd == 0) {
        return res;
    }
    assert(ns > 0 && nd > 0);
    res.reserve(ns + nd - 1);

    long i = 0;
    long j = 0;
    long rem_s = supply[0].second;
    long rem_d = demand[0].second;
    assert(rem_s > 0 && rem_d > 0);

    while (i < ns && j < nd) {
        long z = std::min(rem_s, rem_d);
        res.push_back({supply[i].first, demand[j].first, z});
        rem_s -= z;
        rem_d -= z;
        if (rem_s == 0) {
            ++i;
            if (i < ns) {
                rem_s = supply[i].second;
                assert(rem_s > 0);
            }
        }
        if (rem_d == 0) {
            ++j;
            if (j < nd) {
                rem_d = demand[j].second;
                assert(rem_d > 0);
            }
        }
    }
    assert(i == ns && j == nd);
    return res;
}

// 総和が 0 の収支列 d から、全員の収支を 0 にする最小回数の移動列を構築する
// 計算量: 時間 O(N 2^N), 空間 O(2^N)
// 事前条件: 0 <= N <= 25, sum(d) == 0
inline std::vector<Transfer> min_balance_transfers(const std::vector<long>& d) {
    long n = std::ssize(d);
    assert(0 <= n && n <= 25);
    __int128_t total = 0;
    for (long x : d) {
        total += x;
    }
    assert(total == 0);

    auto groups = zero_sum_partition(d);
    std::vector<Transfer> res;
    for (const auto& g : groups) {
        std::vector<std::pair<long, long>> supply, demand;
        for (long idx : g) {
            if (d[idx] > 0) {
                supply.push_back({idx, d[idx]});
            } else if (d[idx] < 0) {
                demand.push_back({idx, -d[idx]});
            }
        }
        auto part = greedy_transport(supply, demand);
        res.insert(res.end(), part.begin(), part.end());
    }
    return res;
}