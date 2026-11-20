#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

namespace internal {

template <class T>
struct FenwickTree {
    FenwickTree(long n) : n_(n), data_(n + 1, T{0}) {}

    void add(long idx, T val) {
        for (long i = idx + 1; i <= n_; i += i & -i) {
            data_[i] += val;
        }
    }

    // [0, r)
    T sum(long r) const {
        T res = 0;
        for (long i = r; i > 0; i -= i & -i) {
            res += data_[i];
        }
        return res;
    }

    // [l, r)
    T sum(long l, long r) const {
        if (l >= r) return T{0};
        return sum(r) - sum(l);
    }

   private:
    long n_;
    std::vector<T> data_;
};

template <class T, class Op>
struct SegTree {
    SegTree(long n, T id, Op op) : n_(n), id_(id), op_(op), data_(2 * n, id) {}

    void merge_at(long p, T val) {
        p += n_;
        data_[p] = op_(data_[p], val);
        for (p >>= 1; p > 0; p >>= 1) {
            data_[p] = op_(data_[2 * p], data_[2 * p + 1]);
        }
    }

    // [l, r)
    T prod(long l, long r) const {
        T res_l = id_;
        T res_r = id_;
        for (l += n_, r += n_; l < r; l >>= 1, r >>= 1) {
            if (l & 1) res_l = op_(res_l, data_[l++]);
            if (r & 1) res_r = op_(data_[--r], res_r);
        }
        return op_(res_l, res_r);
    }

   private:
    long n_;
    T id_;
    Op op_;
    std::vector<T> data_;
};

struct CompressedPrefix {
    long n;
    long m;
    std::vector<__int128_t> s;
    std::vector<__int128_t> vals;
    std::vector<long> comp;
};

template <class T>
CompressedPrefix build_compressed(const std::vector<T>& a) {
    long n = std::ssize(a);
    std::vector<__int128_t> s(n + 1, 0);
    for (long i = 0; i < n; ++i) {
        s[i + 1] = s[i] + static_cast<__int128_t>(a[i]);
    }

    std::vector<__int128_t> vals = s;
    std::sort(vals.begin(), vals.end());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
    long m = std::ssize(vals);

    std::vector<long> comp(n + 1);
    for (long i = 0; i <= n; ++i) {
        comp[i] = std::lower_bound(vals.begin(), vals.end(), s[i]) - vals.begin();
    }
    return {n, m, std::move(s), std::move(vals), std::move(comp)};
}

// s_r - sum_r <= vals[idx] <= s_r - sum_l を満たす座標圧縮インデックスの範囲 [lo, hi)
inline std::pair<long, long> query_range(const std::vector<__int128_t>& vals, __int128_t s_r, long sum_l, long sum_r) {
    __int128_t lo_val = s_r - static_cast<__int128_t>(sum_r);
    __int128_t hi_val = s_r - static_cast<__int128_t>(sum_l);
    long lo = std::lower_bound(vals.begin(), vals.end(), lo_val) - vals.begin();
    long hi = std::upper_bound(vals.begin(), vals.end(), hi_val) - vals.begin();
    return {lo, hi};
}

}  // namespace internal

// 総和が [sum_l, sum_r] かつ 長さが [min_len, max_len] の連続部分列の個数
// 計算量: 時間 O(N log N), 空間 O(N)
// 事前条件: sum_l <= sum_r, 0 <= min_len <= max_len
template <class T>
long count_subarray_sum_in_range(const std::vector<T>& a, long sum_l, long sum_r, long min_len = 1, long max_len = std::numeric_limits<long>::max()) {
    assert(sum_l <= sum_r);
    assert(0 <= min_len && min_len <= max_len);

    auto cp = internal::build_compressed(a);
    if (min_len > cp.n) return 0;
    max_len = std::min(max_len, cp.n);

    if (sum_l == sum_r) {
        std::vector<long> cnt(cp.m, 0);
        long ans = 0;
        for (long r = 0; r <= cp.n; ++r) {
            if (long add_l = r - min_len; add_l >= 0) {
                cnt[cp.comp[add_l]]++;
            }
            if (long rem_l = r - max_len - 1; rem_l >= 0) {
                cnt[cp.comp[rem_l]]--;
            }
            __int128_t target = cp.s[r] - static_cast<__int128_t>(sum_l);
            auto it = std::lower_bound(cp.vals.begin(), cp.vals.end(), target);
            if (it != cp.vals.end() && *it == target) {
                ans += cnt[it - cp.vals.begin()];
            }
        }
        return ans;
    }

    internal::FenwickTree<long> bit(cp.m);
    long ans = 0;
    for (long r = 0; r <= cp.n; ++r) {
        if (long add_l = r - min_len; add_l >= 0) {
            bit.add(cp.comp[add_l], 1);
        }
        if (long rem_l = r - max_len - 1; rem_l >= 0) {
            bit.add(cp.comp[rem_l], -1);
        }
        auto [lo, hi] = internal::query_range(cp.vals, cp.s[r], sum_l, sum_r);
        ans += bit.sum(lo, hi);
    }
    return ans;
}

// 総和がちょうど target かつ 長さが [min_len, max_len] の連続部分列の個数
// 計算量: 時間 O(N log N), 空間 O(N)
// 事前条件: 0 <= min_len <= max_len
template <class T>
long count_subarray_sum(const std::vector<T>& a, long target = 0, long min_len = 1, long max_len = std::numeric_limits<long>::max()) {
    return count_subarray_sum_in_range(a, target, target, min_len, max_len);
}

// 総和が [sum_l, sum_r] かつ 長さが [min_len, max_len] の連続部分列の長さの総和
// 計算量: 時間 O(N log N), 空間 O(N)
// 事前条件: sum_l <= sum_r, 0 <= min_len <= max_len
template <class T>
long sum_len_subarray_sum_in_range(const std::vector<T>& a, long sum_l, long sum_r, long min_len = 1, long max_len = std::numeric_limits<long>::max()) {
    assert(sum_l <= sum_r);
    assert(0 <= min_len && min_len <= max_len);

    auto cp = internal::build_compressed(a);
    if (min_len > cp.n) return 0;
    max_len = std::min(max_len, cp.n);

    internal::FenwickTree<long> bit_cnt(cp.m);
    internal::FenwickTree<long> bit_l(cp.m);
    long ans = 0;
    for (long r = 0; r <= cp.n; ++r) {
        if (long add_l = r - min_len; add_l >= 0) {
            bit_cnt.add(cp.comp[add_l], 1);
            bit_l.add(cp.comp[add_l], add_l);
        }
        if (long rem_l = r - max_len - 1; rem_l >= 0) {
            bit_cnt.add(cp.comp[rem_l], -1);
            bit_l.add(cp.comp[rem_l], -rem_l);
        }
        auto [lo, hi] = internal::query_range(cp.vals, cp.s[r], sum_l, sum_r);
        long cnt = bit_cnt.sum(lo, hi);
        long sum_idx = bit_l.sum(lo, hi);
        ans += cnt * r - sum_idx;
    }
    return ans;
}

// 総和が [sum_l, sum_r] かつ 長さが min_len 以上の連続部分列の最大長 (存在しない場合は -1)
// 計算量: 時間 O(N log N), 空間 O(N)
// 事前条件: sum_l <= sum_r, min_len >= 0
template <class T>
long max_len_subarray_sum_in_range(const std::vector<T>& a, long sum_l, long sum_r, long min_len = 1) {
    assert(sum_l <= sum_r);
    assert(min_len >= 0);

    auto cp = internal::build_compressed(a);
    if (min_len > cp.n) return -1;

    constexpr long inf = std::numeric_limits<long>::max() / 4;
    internal::SegTree seg(cp.m, inf, [](long x, long y) { return std::min(x, y); });
    long ans = -1;
    for (long r = 0; r <= cp.n; ++r) {
        if (long add_l = r - min_len; add_l >= 0) {
            seg.merge_at(cp.comp[add_l], add_l);
        }
        auto [lo, hi] = internal::query_range(cp.vals, cp.s[r], sum_l, sum_r);
        if (lo < hi) {
            long min_l = seg.prod(lo, hi);
            if (min_l != inf) {
                ans = std::max(ans, r - min_l);
            }
        }
    }
    return ans;
}

// 総和が [sum_l, sum_r] かつ 長さが min_len 以上の連続部分列の最小長 (存在しない場合は -1)
// 計算量: 時間 O(N log N), 空間 O(N)
// 事前条件: sum_l <= sum_r, min_len >= 0
template <class T>
long min_len_subarray_sum_in_range(const std::vector<T>& a, long sum_l, long sum_r, long min_len = 1) {
    assert(sum_l <= sum_r);
    assert(min_len >= 0);

    auto cp = internal::build_compressed(a);
    if (min_len > cp.n) return -1;

    constexpr long inf = std::numeric_limits<long>::max() / 4;
    internal::SegTree seg(cp.m, -1, [](long x, long y) { return std::max(x, y); });
    long ans = inf;
    for (long r = 0; r <= cp.n; ++r) {
        if (long add_l = r - min_len; add_l >= 0) {
            seg.merge_at(cp.comp[add_l], add_l);
        }
        auto [lo, hi] = internal::query_range(cp.vals, cp.s[r], sum_l, sum_r);
        if (lo < hi) {
            long max_l = seg.prod(lo, hi);
            if (max_l != -1) {
                ans = std::min(ans, r - max_l);
            }
        }
    }
    return (ans == inf) ? -1 : ans;
}

// 総和が [sum_l, sum_r] かつ 長さが [min_len, max_len] の連続部分列 [left, right) を 1 つ復元 (存在しない場合は {-1, -1})
// 計算量: 時間 O(N log N), 空間 O(N)
// 事前条件: sum_l <= sum_r, 0 <= min_len <= max_len
template <class T>
std::pair<long, long> find_subarray_sum_in_range(const std::vector<T>& a, long sum_l, long sum_r, long min_len = 1, long max_len = std::numeric_limits<long>::max()) {
    assert(sum_l <= sum_r);
    assert(0 <= min_len && min_len <= max_len);

    auto cp = internal::build_compressed(a);
    if (min_len > cp.n) return {-1, -1};
    max_len = std::min(max_len, cp.n);

    internal::FenwickTree<long> bit(cp.m);
    for (long r = 0; r <= cp.n; ++r) {
        if (long add_l = r - min_len; add_l >= 0) {
            bit.add(cp.comp[add_l], 1);
        }
        if (long rem_l = r - max_len - 1; rem_l >= 0) {
            bit.add(cp.comp[rem_l], -1);
        }
        auto [lo, hi] = internal::query_range(cp.vals, cp.s[r], sum_l, sum_r);
        if (bit.sum(lo, hi) > 0) {
            long l_start = std::max(0L, r - max_len);
            long l_end = r - min_len;
            for (long left = l_start; left <= l_end; ++left) {
                __int128_t diff = cp.s[r] - cp.s[left];
                if (sum_l <= diff && diff <= sum_r) {
                    return {left, r};
                }
            }
        }
    }
    return {-1, -1};
}