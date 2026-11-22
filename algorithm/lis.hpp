#pragma once

#include <algorithm>
#include <functional>
#include <iterator>
#include <utility>
#include <vector>

// 最長増加部分列 (LIS) の長さを求める
// 時間計算量: O(N log N), 空間計算量: O(N)
// 事前条件: Comp は strict weak ordering を満たすこと
template <class T, class Comp = std::less<T>>
long lis_length(const std::vector<T>& a, bool strict = true, Comp comp = Comp{}) {
    std::vector<T> tails;
    for (const auto& x : a) {
        auto it = strict ? std::lower_bound(tails.begin(), tails.end(), x, comp) : std::upper_bound(tails.begin(), tails.end(), x, comp);
        if (it == tails.end()) {
            tails.push_back(x);
        } else {
            *it = x;
        }
    }
    return std::ssize(tails);
}

// 各位置 i を末尾とする最長増加部分列 (LIS) の長さを求める
// 時間計算量: O(N log N), 空間計算量: O(N)
// 事前条件: Comp は strict weak ordering を満たすこと
template <class T, class Comp = std::less<T>>
std::vector<long> lis_lengths(const std::vector<T>& a, bool strict = true, Comp comp = Comp{}) {
    long n = std::ssize(a);
    std::vector<long> res(n);
    std::vector<T> tails;
    for (long i = 0; i < n; ++i) {
        auto it = strict ? std::lower_bound(tails.begin(), tails.end(), a[i], comp) : std::upper_bound(tails.begin(), tails.end(), a[i], comp);
        long pos = it - tails.begin();
        if (it == tails.end()) {
            tails.push_back(a[i]);
        } else {
            *it = a[i];
        }
        res[i] = pos + 1;
    }
    return res;
}

// 各位置 i を先頭とする最長増加部分列 (LIS) の長さを求める
// 時間計算量: O(N log N), 空間計算量: O(N)
// 事前条件: Comp は strict weak ordering を満たすこと
template <class T, class Comp = std::less<T>>
std::vector<long> lis_suffix_lengths(const std::vector<T>& a, bool strict = true, Comp comp = Comp{}) {
    std::vector<T> rev(a.rbegin(), a.rend());
    auto rev_comp = [&](const T& x, const T& y) { return comp(y, x); };
    auto res = lis_lengths(rev, strict, rev_comp);
    std::reverse(res.begin(), res.end());
    return res;
}

// 最長増加部分列 (LIS) を構成する要素の添字列 (0-indexed) を復元する
// 時間計算量: O(N log N), 空間計算量: O(N)
// 事前条件: Comp は strict weak ordering を満たすこと
template <class T, class Comp = std::less<T>>
std::vector<long> lis_restore(const std::vector<T>& a, bool strict = true, Comp comp = Comp{}) {
    long n = std::ssize(a);
    if (n == 0) return {};
    std::vector<T> tails_val;
    std::vector<long> tails_idx;
    std::vector<long> prev(n, -1);
    for (long i = 0; i < n; ++i) {
        auto it = strict ? std::lower_bound(tails_val.begin(), tails_val.end(), a[i], comp) : std::upper_bound(tails_val.begin(), tails_val.end(), a[i], comp);
        long pos = it - tails_val.begin();
        if (it == tails_val.end()) {
            tails_val.push_back(a[i]);
            tails_idx.push_back(i);
        } else {
            *it = a[i];
            tails_idx[pos] = i;
        }
        if (pos > 0) prev[i] = tails_idx[pos - 1];
    }
    std::vector<long> res;
    res.reserve(tails_idx.size());
    for (long cur = tails_idx.back(); cur != -1; cur = prev[cur]) {
        res.push_back(cur);
    }
    std::reverse(res.begin(), res.end());
    return res;
}

// 各要素がいずれかの最長増加部分列 (LIS) に含まれるかを判定する
// 時間計算量: O(N log N), 空間計算量: O(N)
// 事前条件: Comp は strict weak ordering を満たすこと
template <class T, class Comp = std::less<T>>
std::vector<bool> lis_in_any(const std::vector<T>& a, bool strict = true, Comp comp = Comp{}) {
    long n = std::ssize(a);
    if (n == 0) return {};
    auto fwd = lis_lengths(a, strict, comp);
    auto bwd = lis_suffix_lengths(a, strict, comp);
    long max_len = *std::max_element(fwd.begin(), fwd.end());
    std::vector<bool> res(n, false);
    for (long i = 0; i < n; ++i) {
        if (fwd[i] + bwd[i] - 1 == max_len) res[i] = true;
    }
    return res;
}

// 各要素が全ての最長増加部分列 (LIS) に含まれるかを判定する
// 時間計算量: O(N log N), 空間計算量: O(N)
// 事前条件: Comp は strict weak ordering を満たすこと
template <class T, class Comp = std::less<T>>
std::vector<bool> lis_in_all(const std::vector<T>& a, bool strict = true, Comp comp = Comp{}) {
    long n = std::ssize(a);
    if (n == 0) return {};
    auto fwd = lis_lengths(a, strict, comp);
    auto bwd = lis_suffix_lengths(a, strict, comp);
    long max_len = *std::max_element(fwd.begin(), fwd.end());
    std::vector<long> cnt(max_len + 1, 0);
    for (long i = 0; i < n; ++i) {
        if (fwd[i] + bwd[i] - 1 == max_len) ++cnt[fwd[i]];
    }
    std::vector<bool> res(n, false);
    for (long i = 0; i < n; ++i) {
        if (fwd[i] + bwd[i] - 1 == max_len && cnt[fwd[i]] == 1) res[i] = true;
    }
    return res;
}

// 最長増加部分列 (LIS) の個数を求める（a が空の場合は 0）
// 時間計算量: O(N log N), 空間計算量: O(N)
// 事前条件: Comp は strict weak ordering を満たすこと、Mint は加減算が可能であること
template <class Mint = long, class T, class Comp = std::less<T>>
Mint lis_count(const std::vector<T>& a, bool strict = true, Comp comp = Comp{}) {
    if (a.empty()) return Mint(0);
    std::vector<T> tails;
    std::vector<std::vector<std::pair<T, Mint>>> buckets;
    for (const auto& x : a) {
        auto it = strict ? std::lower_bound(tails.begin(), tails.end(), x, comp) : std::upper_bound(tails.begin(), tails.end(), x, comp);
        long pos = it - tails.begin();
        if (it == tails.end()) {
            tails.push_back(x);
            buckets.emplace_back();
        } else {
            *it = x;
        }
        Mint ways(0);
        if (pos == 0) {
            ways = Mint(1);
        } else {
            const auto& prev = buckets[pos - 1];
            auto bound = std::partition_point(prev.begin(), prev.end(), [&](const auto& p) { return strict ? !comp(p.first, x) : comp(x, p.first); });
            long idx = bound - prev.begin();
            Mint total = prev.back().second;
            ways = (idx == 0) ? total : total - prev[idx - 1].second;
        }
        Mint sum = buckets[pos].empty() ? ways : buckets[pos].back().second + ways;
        buckets[pos].emplace_back(x, sum);
    }
    return buckets.back().back().second;
};
