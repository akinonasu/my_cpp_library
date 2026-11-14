#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <limits>
#include <ranges>
#include <utility>
#include <vector>

namespace zsr_detail {

template <typename T>
struct FenwickTree {
    int n;
    std::vector<T> data;

    explicit FenwickTree(int n) : n(n), data(n + 1, T{0}) {}

    void add(int idx, T val) {
        for (int i = idx + 1; i <= n; i += i & -i) {
            data[i] += val;
        }
    }

    // [0, r)
    T sum(int r) const {
        T res = 0;
        for (int i = r; i > 0; i -= i & -i) {
            res += data[i];
        }
        return res;
    }

    // [l, r)
    T sum(int l, int r) const {
        if (l >= r) return T{0};
        return sum(r) - sum(l);
    }
};

template <typename T, typename Op>
struct SegTree {
    int n;
    T id;
    Op op;
    std::vector<T> data;

    SegTree(int n, T id, Op op) : n(n), id(id), op(op), data(2 * n, id) {}

    void merge_at(int p, T val) {
        p += n;
        data[p] = op(data[p], val);
        for (p >>= 1; p > 0; p >>= 1) {
            data[p] = op(data[2 * p], data[2 * p + 1]);
        }
    }

    // [l, r)
    T prod(int l, int r) const {
        T res_l = id;
        T res_r = id;
        for (l += n, r += n; l < r; l >>= 1, r >>= 1) {
            if (l & 1) res_l = op(res_l, data[l++]);
            if (r & 1) res_r = op(data[--r], res_r);
        }
        return op(res_l, res_r);
    }
};

struct CompressedPrefix {
    int n;
    int m;
    std::vector<long> s;
    std::vector<long> vals;
    std::vector<int> comp;
};

template <std::ranges::input_range Rng>
    requires std::convertible_to<std::ranges::range_value_t<Rng>, long>
inline CompressedPrefix build_compressed(const Rng& a) {
    std::vector<long> s;
    if constexpr (std::ranges::sized_range<Rng>) {
        s.reserve(std::ranges::size(a) + 1);
    }
    s.push_back(0);
    for (auto&& x : a) {
        s.push_back(s.back() + static_cast<long>(x));
    }
    assert(s.size() - 1 <= static_cast<std::size_t>(std::numeric_limits<int>::max()));
    int n = static_cast<int>(s.size()) - 1;

    std::vector<long> vals = s;
    std::ranges::sort(vals);
    auto ret = std::ranges::unique(vals);
    vals.erase(ret.begin(), ret.end());
    int m = static_cast<int>(vals.size());

    std::vector<int> comp(n + 1);
    for (int i = 0; i <= n; ++i) {
        comp[i] = static_cast<int>(std::ranges::lower_bound(vals, s[i]) - vals.begin());
    }
    return {n, m, std::move(s), std::move(vals), std::move(comp)};
}

// a - b をオーバーフロー安全に計算する
// 戻り値: {dir, val}
// dir == -1: a - b < LONG_MIN (アンダーフロー)
// dir ==  1: a - b > LONG_MAX (オーバーフロー)
// dir ==  0: 正常値 val = a - b
inline std::pair<int, long> safe_sub(long a, long b) {
    constexpr long LMIN = std::numeric_limits<long>::min();
    constexpr long LMAX = std::numeric_limits<long>::max();
    if (b > 0 && a < LMIN + b) return {-1, 0};
    if (b < 0 && a > LMAX + b) return {1, 0};
    return {0, a - b};
}

// s_r - sum_r <= vals[idx] <= s_r - sum_l を満たす座標圧縮インデックスの範囲 [lo, hi)
inline std::pair<int, int> query_range(const std::vector<long>& vals, long s_r, long sum_l, long sum_r) {
    auto [lo_ov, lo_val] = safe_sub(s_r, sum_r);
    auto [hi_ov, hi_val] = safe_sub(s_r, sum_l);

    if (lo_ov > 0 || hi_ov < 0) {
        return {0, 0};
    }
    if ((lo_ov == 0 && lo_val > vals.back()) || (hi_ov == 0 && hi_val < vals.front())) {
        return {0, 0};
    }

    int m = static_cast<int>(vals.size());
    int lo = (lo_ov < 0 || lo_val <= vals.front()) ? 0 : static_cast<int>(std::ranges::lower_bound(vals, lo_val) - vals.begin());
    int hi = (hi_ov > 0 || hi_val >= vals.back()) ? m : static_cast<int>(std::ranges::upper_bound(vals, hi_val) - vals.begin());
    return {lo, hi};
}

inline bool in_sum_range(long s_l, long s_r, long sum_l, long sum_r) {
    auto [ov, diff] = safe_sub(s_r, s_l);
    if (ov != 0) return false;
    return sum_l <= diff && diff <= sum_r;
}

}  // namespace zsr_detail

// 総和が [sum_l, sum_r] かつ 長さが [min_len, max_len] の連続部分列の個数
template <std::ranges::input_range Rng>
    requires std::convertible_to<std::ranges::range_value_t<Rng>, long>
inline long count_subarray_sum_in_range(const Rng& a, long sum_l, long sum_r, int min_len = 1, int max_len = std::numeric_limits<int>::max()) {
    assert(sum_l <= sum_r);
    assert(0 <= min_len && min_len <= max_len);

    auto cp = zsr_detail::build_compressed(a);
    if (min_len > cp.n) return 0;
    max_len = std::min(max_len, cp.n);

    if (sum_l == sum_r) {
        std::vector<int> cnt(cp.m, 0);
        long ans = 0;
        for (int r = 0; r <= cp.n; ++r) {
            if (int add_l = r - min_len; add_l >= 0) {
                cnt[cp.comp[add_l]]++;
            }
            if (int rem_l = r - max_len - 1; rem_l >= 0) {
                cnt[cp.comp[rem_l]]--;
            }
            auto [ov, t] = zsr_detail::safe_sub(cp.s[r], sum_l);
            if (ov == 0 && t >= cp.vals.front() && t <= cp.vals.back()) {
                auto it = std::ranges::lower_bound(cp.vals, t);
                if (it != cp.vals.end() && *it == t) {
                    ans += cnt[static_cast<int>(it - cp.vals.begin())];
                }
            }
        }
        return ans;
    }

    zsr_detail::FenwickTree<int> bit(cp.m);
    long ans = 0;
    for (int r = 0; r <= cp.n; ++r) {
        if (int add_l = r - min_len; add_l >= 0) {
            bit.add(cp.comp[add_l], 1);
        }
        if (int rem_l = r - max_len - 1; rem_l >= 0) {
            bit.add(cp.comp[rem_l], -1);
        }
        auto [lo, hi] = zsr_detail::query_range(cp.vals, cp.s[r], sum_l, sum_r);
        ans += bit.sum(lo, hi);
    }
    return ans;
}

// 総和がちょうど target かつ 長さが [min_len, max_len] の連続部分列の個数
template <std::ranges::input_range Rng>
    requires std::convertible_to<std::ranges::range_value_t<Rng>, long>
inline long count_subarray_sum(const Rng& a, long target = 0, int min_len = 1, int max_len = std::numeric_limits<int>::max()) {
    return count_subarray_sum_in_range(a, target, target, min_len, max_len);
}

// 総和が [sum_l, sum_r] かつ 長さが [min_len, max_len] の連続部分列の長さの総和
template <std::ranges::input_range Rng>
    requires std::convertible_to<std::ranges::range_value_t<Rng>, long>
inline long sum_len_subarray_sum_in_range(const Rng& a, long sum_l, long sum_r, int min_len = 1, int max_len = std::numeric_limits<int>::max()) {
    assert(sum_l <= sum_r);
    assert(0 <= min_len && min_len <= max_len);

    auto cp = zsr_detail::build_compressed(a);
    if (min_len > cp.n) return 0;
    max_len = std::min(max_len, cp.n);

    zsr_detail::FenwickTree<int> bit_cnt(cp.m);
    zsr_detail::FenwickTree<long> bit_l(cp.m);
    long ans = 0;
    for (int r = 0; r <= cp.n; ++r) {
        if (int add_l = r - min_len; add_l >= 0) {
            bit_cnt.add(cp.comp[add_l], 1);
            bit_l.add(cp.comp[add_l], add_l);
        }
        if (int rem_l = r - max_len - 1; rem_l >= 0) {
            bit_cnt.add(cp.comp[rem_l], -1);
            bit_l.add(cp.comp[rem_l], -rem_l);
        }
        auto [lo, hi] = zsr_detail::query_range(cp.vals, cp.s[r], sum_l, sum_r);
        long cnt = bit_cnt.sum(lo, hi);
        long sum_idx = bit_l.sum(lo, hi);
        ans += cnt * r - sum_idx;
    }
    return ans;
}

// 総和が [sum_l, sum_r] かつ 長さが min_len 以上の連続部分列の最大長 (存在しない場合は -1)
template <std::ranges::input_range Rng>
    requires std::convertible_to<std::ranges::range_value_t<Rng>, long>
inline int max_len_subarray_sum_in_range(const Rng& a, long sum_l, long sum_r, int min_len = 1) {
    assert(sum_l <= sum_r);
    assert(min_len >= 0);

    auto cp = zsr_detail::build_compressed(a);
    if (min_len > cp.n) return -1;

    constexpr int INF = std::numeric_limits<int>::max();
    zsr_detail::SegTree seg(cp.m, INF, [](int x, int y) { return std::min(x, y); });
    int ans = -1;
    for (int r = 0; r <= cp.n; ++r) {
        if (int add_l = r - min_len; add_l >= 0) {
            seg.merge_at(cp.comp[add_l], add_l);
        }
        auto [lo, hi] = zsr_detail::query_range(cp.vals, cp.s[r], sum_l, sum_r);
        if (lo < hi) {
            int min_l = seg.prod(lo, hi);
            if (min_l != INF) {
                ans = std::max(ans, r - min_l);
            }
        }
    }
    return ans;
}

// 総和が [sum_l, sum_r] かつ 長さが min_len 以上の連続部分列の最小長 (存在しない場合は -1)
template <std::ranges::input_range Rng>
    requires std::convertible_to<std::ranges::range_value_t<Rng>, long>
inline int min_len_subarray_sum_in_range(const Rng& a, long sum_l, long sum_r, int min_len = 1) {
    assert(sum_l <= sum_r);
    assert(min_len >= 0);

    auto cp = zsr_detail::build_compressed(a);
    if (min_len > cp.n) return -1;

    constexpr int INF = std::numeric_limits<int>::max();
    zsr_detail::SegTree seg(cp.m, -1, [](int x, int y) { return std::max(x, y); });
    int ans = INF;
    for (int r = 0; r <= cp.n; ++r) {
        if (int add_l = r - min_len; add_l >= 0) {
            seg.merge_at(cp.comp[add_l], add_l);
        }
        auto [lo, hi] = zsr_detail::query_range(cp.vals, cp.s[r], sum_l, sum_r);
        if (lo < hi) {
            int max_l = seg.prod(lo, hi);
            if (max_l != -1) {
                ans = std::min(ans, r - max_l);
            }
        }
    }
    return (ans == INF) ? -1 : ans;
}

// 総和が [sum_l, sum_r] かつ 長さが [min_len, max_len] の連続部分列 [left, right) を 1 つ復元 (存在しない場合は {-1, -1})
template <std::ranges::input_range Rng>
    requires std::convertible_to<std::ranges::range_value_t<Rng>, long>
inline std::pair<int, int> find_subarray_sum_in_range(const Rng& a, long sum_l, long sum_r, int min_len = 1, int max_len = std::numeric_limits<int>::max()) {
    assert(sum_l <= sum_r);
    assert(0 <= min_len && min_len <= max_len);

    auto cp = zsr_detail::build_compressed(a);
    if (min_len > cp.n) return {-1, -1};
    max_len = std::min(max_len, cp.n);

    zsr_detail::FenwickTree<int> bit(cp.m);
    for (int r = 0; r <= cp.n; ++r) {
        if (int add_l = r - min_len; add_l >= 0) {
            bit.add(cp.comp[add_l], 1);
        }
        if (int rem_l = r - max_len - 1; rem_l >= 0) {
            bit.add(cp.comp[rem_l], -1);
        }
        auto [lo, hi] = zsr_detail::query_range(cp.vals, cp.s[r], sum_l, sum_r);
        if (bit.sum(lo, hi) > 0) {
            int l_start = std::max(0, r - max_len);
            int l_end = r - min_len;
            for (int left = l_start; left <= l_end; ++left) {
                if (zsr_detail::in_sum_range(cp.s[left], cp.s[r], sum_l, sum_r)) {
                    return {left, r};
                }
            }
        }
    }
    return {-1, -1};
}