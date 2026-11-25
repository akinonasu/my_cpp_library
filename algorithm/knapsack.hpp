#pragma once

#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

struct KnapsackItem {
    long w, v, m = 1;
};

struct KnapsackPowItem {
    long t, v;
};

struct KnapsackScaledItem {
    long t, s, v, m = 1;
};

namespace internal {

template <class T>
void concave_max_plus(long step, const std::vector<T>& gains, std::vector<T>& dp, T neg_inf) {
    assert(step > 0);
    assert(!gains.empty());
    long z = std::ssize(dp) - 1;
    T uninit = neg_inf * 4;
    for (long s = 0; s < step && s <= z; s++) {
        long n = (z - s) / step + 1;
        std::vector<long> idx(n + 1);
        std::vector<T> d0(n), d1(n, uninit);
        for (long i = 0; i < n; i++)
            d0[i] = dp[i * step + s];
        d1[0] = d0[0] + gains[0];
        idx[n] = n - 1;
        long bi = 1;
        while (bi < n)
            bi *= 2;
        for (long q = bi / 2; q > 0; q /= 2) {
            for (long m = q; m < n; m += q * 2) {
                long l = m - q;
                long r = std::min(m + q, n);
                idx[m] = idx[l];
                for (long t = idx[l]; t <= idx[r]; t++) {
                    if (t <= m && m - t < std::ssize(gains)) {
                        T val = d0[t] + gains[m - t];
                        if (d1[m] < val) {
                            d1[m] = val;
                            idx[m] = t;
                        }
                    }
                }
            }
        }
        for (long i = 0; i < n; i++)
            dp[i * step + s] = std::max(neg_inf, d1[i]);
    }
}

}  // namespace internal

// [knapsack_1] 個数制限付きナップサック（愚直展開）
// 時間計算量: O(Z \sum min(m_i, Z / w_i)), 空間計算量: O(Z)
// 事前条件: z >= 0, w_i >= 0, m_i >= 0
inline long knapsack_bounded_naive(long z, const std::vector<KnapsackItem>& items) {
    assert(z >= 0);
    long base = 0;
    std::vector<long> dp(z + 1, 0);
    for (const auto& a : items) {
        assert(a.w >= 0 && a.m >= 0);
        if (a.m == 0 || a.v <= 0 || a.w > z) continue;
        if (a.w == 0) {
            base += a.v * a.m;
            continue;
        }
        long cnt = std::min(a.m, z / a.w);
        for (long mi = 0; mi < cnt; mi++) {
            for (long i = z; i >= a.w; i--) {
                dp[i] = std::max(dp[i], dp[i - a.w] + a.v);
            }
        }
    }
    return base + dp[z];
}

// [knapsack_2] 個数無制限ナップサック（m は無視）
// 時間計算量: O(N Z), 空間計算量: O(Z)
// 事前条件: z >= 0, w_i >= 0, (w_i == 0 かつ v_i > 0 は禁止)
inline long knapsack_unbounded(long z, const std::vector<KnapsackItem>& items) {
    assert(z >= 0);
    std::vector<long> dp(z + 1, 0);
    for (const auto& a : items) {
        assert(a.w >= 0);
        assert(!(a.w == 0 && a.v > 0));
        if (a.v <= 0 || a.w > z) continue;
        for (long i = a.w; i <= z; i++) {
            dp[i] = std::max(dp[i], dp[i - a.w] + a.v);
        }
    }
    return dp[z];
}

// [knapsack_3] 個数制限付きナップサック（ブロック分割）
// 時間計算量: O(N Z), 空間計算量: O(Z)
// 事前条件: z >= 0, w_i >= 0, m_i >= 0
inline long knapsack_bounded(long z, const std::vector<KnapsackItem>& items) {
    assert(z >= 0);
    long base = 0;
    std::vector<long> dp(z + 1, 0);
    for (const auto& a : items) {
        assert(a.w >= 0 && a.m >= 0);
        if (a.m == 0 || a.v <= 0 || a.w > z) continue;
        if (a.w == 0) {
            base += a.v * a.m;
            continue;
        }
        long m = std::min(a.m, z / a.w);
        std::vector<long> buf = dp;
        for (long s = 0; s * a.w <= z; s += m) {
            long l = s * a.w;
            long r = std::min(z + 1, (s + m) * a.w);
            for (long i = l; i < r - a.w; i++) {
                dp[i + a.w] = std::max(dp[i + a.w], dp[i] + a.v);
            }
            for (long i = r - a.w - 1; i >= l; i--) {
                buf[i] = std::max(buf[i], buf[i + a.w] - a.v);
            }
        }
        for (long i = a.w * m; i <= z; i++) {
            dp[i] = std::max(dp[i], buf[i - a.w * m] + a.v * m);
        }
    }
    return base + dp[z];
}

// [knapsack_4] 個数制限付きナップサック（価値総和 DP）
// 時間計算量: O(Y \sum min(m_i, Z / w_i)) (Y = \sum v_i min(m_i, Z / w_i)), 空間計算量: O(Y)
// 事前条件: z >= 0, w_i >= 0, v_i >= 0, m_i >= 0
inline long knapsack_bounded_by_value(long z, const std::vector<KnapsackItem>& items) {
    assert(z >= 0);
    long base = 0;
    long total_v = 0;
    std::vector<KnapsackItem> cand;
    for (const auto& a : items) {
        assert(a.w >= 0 && a.v >= 0 && a.m >= 0);
        if (a.m == 0 || a.v == 0 || a.w > z) continue;
        if (a.w == 0) {
            base += a.v * a.m;
            continue;
        }
        long cnt = std::min(a.m, z / a.w);
        total_v += a.v * cnt;
        cand.push_back({a.w, a.v, cnt});
    }
    std::vector<long> dp(total_v + 1, z + 1);
    dp[0] = 0;
    for (const auto& a : cand) {
        for (long mi = 0; mi < a.m; mi++) {
            for (long i = total_v; i >= a.v; i--) {
                if (dp[i - a.v] <= z - a.w) {
                    dp[i] = std::min(dp[i], dp[i - a.v] + a.w);
                }
            }
        }
    }
    long ans = 0;
    for (long y = 0; y <= total_v; y++) {
        if (dp[y] <= z) ans = y;
    }
    return base + ans;
}

// [knapsack_5] 0-1 ナップサック（半分全列挙・m は無視）
// 時間計算量: O(2^{N/2}), 空間計算量: O(2^{N/2})
// 事前条件: z >= 0, w_i >= 0
inline long knapsack_01_meet_in_the_middle(long z, const std::vector<KnapsackItem>& items) {
    assert(z >= 0);
    long base = 0;
    std::vector<std::pair<long, long>> wv1 = {{0, 0}}, wv2 = {{0, 0}};
    for (const auto& a : items) {
        assert(a.w >= 0);
        if (a.v <= 0 || a.w > z) continue;
        if (a.w == 0) {
            base += a.v;
            continue;
        }
        std::vector<std::pair<long, long>> tmp;
        long n = std::ssize(wv1);
        long p = 0, q = 0, xv = -1;
        while (q < n) {
            long pw = 0, pv = 0;
            if (p < n && wv1[p].first <= wv1[q].first + a.w) {
                pw = wv1[p].first;
                pv = wv1[p].second;
                p++;
            } else {
                pw = wv1[q].first + a.w;
                pv = wv1[q].second + a.v;
                q++;
            }
            if (pw <= z && xv < pv) {
                xv = pv;
                tmp.push_back({pw, pv});
            }
        }
        std::swap(wv1, tmp);
        std::swap(wv1, wv2);
    }
    long ans = 0;
    for (auto [w, v] : wv1) {
        while (wv2.back().first + w > z)
            wv2.pop_back();
        ans = std::max(ans, v + wv2.back().second);
    }
    return base + ans;
}

// [knapsack_6] 重みが k^t の 0-1 ナップサック
// 時間計算量: O(N \log N + \log_k Z), 空間計算量: O(N + \log_k Z)
// 事前条件: z >= 0, k >= 2, t_i >= 0
inline long knapsack_01_pow_weight(long z, long k, const std::vector<KnapsackPowItem>& items) {
    assert(z >= 0 && k >= 2);
    long s_cnt = 0;
    for (long x = z; x > 0; x /= k)
        s_cnt++;
    std::vector<std::vector<long>> buckets(s_cnt + 1);
    for (const auto& a : items) {
        assert(a.t >= 0);
        if (a.v > 0 && a.t < s_cnt) buckets[a.t].push_back(a.v);
    }
    long ans = 0;
    for (long s = 0; s < s_cnt && z > 0; s++) {
        auto& cur = buckets[s];
        std::sort(cur.begin(), cur.end());
        while (z % k > 0 && !cur.empty()) {
            ans += cur.back();
            cur.pop_back();
            z--;
        }
        while (!cur.empty()) {
            long sum = 0;
            for (long i = 0; i < k && !cur.empty(); i++) {
                sum += cur.back();
                cur.pop_back();
            }
            buckets[s + 1].push_back(sum);
        }
        z /= k;
    }
    return ans;
}

// [knapsack_7] 重みが s * 2^t の個数制限付きナップサック
// 時間計算量: O((\sum s_i) \log Z) 程度, 空間計算量: O(\sum s_i + N)
// 事前条件: z >= 0, t_i >= 0, s_i >= 0, m_i >= 0
inline long knapsack_bounded_scaled_pow2(long z, const std::vector<KnapsackScaledItem>& items) {
    assert(z >= 0);
    constexpr long k = 2;
    long base = 0;
    long max_t = -1;
    for (const auto& a : items) {
        assert(a.t >= 0 && a.s >= 0 && a.m >= 0);
        if (a.m == 0 || a.v <= 0) continue;
        if (a.s == 0) {
            base += a.v * a.m;
            continue;
        }
        max_t = std::max(max_t, a.t);
    }
    if (max_t < 0) return base;

    std::vector<std::vector<KnapsackItem>> ranked(max_t + 1);
    for (const auto& a : items) {
        if (a.m > 0 && a.v > 0 && a.s > 0) {
            ranked[a.t].push_back({a.s, a.v, a.m});
        }
    }

    std::vector<long> dp = {0};
    std::vector<KnapsackItem> carry;
    for (long t = 0; z > 0 && (t <= max_t || !carry.empty()); t++) {
        if (t <= max_t) {
            for (const auto& a : ranked[t])
                carry.push_back(a);
        }
        std::vector<KnapsackItem> next_carry;
        for (auto [w, v, m] : carry) {
            if (w > z) continue;
            m = std::min(m, z / w);
            long xm = m;
            if (m >= k - 1) xm = k - 1 + (m - (k - 1)) % k;
            for (long mi = 0; mi < xm; mi++) {
                long cur_z = std::ssize(dp) - 1;
                for (long i = 0; i < w; i++)
                    dp.push_back(dp.back());
                for (long i = cur_z; i >= 0; i--) {
                    dp[i + w] = std::max(dp[i + w], dp[i] + v);
                }
            }
            long rem_m = (m - xm) / k;
            if (rem_m > 0) next_carry.push_back({w, v * k, rem_m});
        }
        carry = std::move(next_carry);
        std::vector<long> nx;
        for (long i = z % k; i < std::ssize(dp); i += k)
            nx.push_back(dp[i]);
        nx.push_back(dp.back());
        std::swap(dp, nx);
        z /= k;
    }
    return base + dp[std::min(z, std::ssize(dp) - 1)];
}

// [knapsack_8] 巨大 Z・小重み向け個数無制限ナップサック（m は無視）
// 時間計算量: O(N + W_{top}^2), 空間計算量: O(W_{top})
// 事前条件: z >= 0, w_i >= 0, (w_i == 0 かつ v_i > 0 は禁止), かつ z >= (W_{top} - 1) * max(w_i)
inline long knapsack_unbounded_large_z(long z, const std::vector<KnapsackItem>& items) {
    assert(z >= 0);
    std::vector<KnapsackItem> cand;
    for (const auto& a : items) {
        assert(a.w >= 0);
        assert(!(a.w == 0 && a.v > 0));
        if (a.v <= 0 || a.w > z) continue;
        cand.push_back(a);
    }
    if (cand.empty()) return 0;

    for (auto& a : cand) {
        __int128_t lhs = __int128_t(cand.back().v) * a.w;
        __int128_t rhs = __int128_t(a.v) * cand.back().w;
        if (lhs < rhs || (lhs == rhs && a.w < cand.back().w)) {
            std::swap(cand.back(), a);
        }
    }
    auto top = cand.back();
    cand.pop_back();

    std::vector<long> gain(top.w, 0);
    for (const auto& a : cand) {
        long r = a.w % top.w;
        long val = a.v - top.v * (a.w / top.w);
        gain[r] = std::max(gain[r], val);
    }

    constexpr long inf = std::numeric_limits<long>::max() / 4;
    std::vector<long> dp(top.w, -inf);
    dp[0] = 0;
    for (long r = 1; r < top.w; r++) {
        if (gain[r] <= 0) continue;
        long g = std::gcd(r, top.w);
        long xg = top.w / g;
        for (long s = 0; s < g; s++) {
            long p = s;
            for (long i = 0; i < xg * 2; i++) {
                long q = p + r;
                if (q < top.w) {
                    if (dp[p] > -inf) dp[q] = std::max(dp[q], dp[p] + gain[r]);
                } else {
                    q -= top.w;
                    if (dp[p] > -inf) dp[q] = std::max(dp[q], dp[p] + gain[r] - top.v);
                }
                p = q;
            }
        }
    }

    long q = z / top.w;
    long rem = z % top.w;
    long ans = 0;
    for (long r = 0; r < top.w; r++) {
        if (dp[r] == -inf) continue;
        if (r <= rem) {
            ans = std::max(ans, dp[r] + q * top.v);
        } else if (q >= 1) {
            ans = std::max(ans, dp[r] + (q - 1) * top.v);
        }
    }
    return ans;
}

// [knapsack_9] 個数制限付きナップサック（重み別グルーピング＋単調性 Max-Plus）
// 時間計算量: O(N \log N + D Z \log Z) (D は相異なる w_i の個数), 空間計算量: O(Z + N)
// 事前条件: z >= 0, w_i >= 0, m_i >= 0
inline long knapsack_bounded_monotone(long z, const std::vector<KnapsackItem>& items) {
    assert(z >= 0);
    constexpr long inf = std::numeric_limits<long>::max() / 8;
    long base = 0;
    std::vector<std::vector<KnapsackItem>> ranked(z + 1);
    for (const auto& a : items) {
        assert(a.w >= 0 && a.m >= 0);
        if (a.m == 0 || a.v <= 0 || a.w > z) continue;
        if (a.w == 0) {
            base += a.v * a.m;
            continue;
        }
        ranked[a.w].push_back(a);
    }
    std::vector<long> dp(z + 1, 0);
    for (long w = 1; w <= z; w++) {
        if (ranked[w].empty()) continue;
        std::sort(ranked[w].begin(), ranked[w].end(), [](const KnapsackItem& l, const KnapsackItem& r) { return l.v > r.v; });
        long sz = z / w;
        std::vector<long> gains = {0};
        for (const auto& a : ranked[w]) {
            for (long i = 0; i < a.m && std::ssize(gains) <= sz; i++) {
                gains.push_back(gains.back() + a.v);
            }
        }
        internal::concave_max_plus(w, gains, dp, -inf);
    }
    return base + dp[z];
}

// [knapsack_10] 小重み・巨大 Z, m_i 向け個数制限付きナップサック
// 時間計算量: O(N \log N + W_{max}^3 \log W_{max}), 空間計算量: O(W_{max}^2 + N)
// 事前条件: z >= 0, w_i >= 0, m_i >= 0
inline long knapsack_bounded_small_weight(long z, const std::vector<KnapsackItem>& items) {
    assert(z >= 0);
    long base = 0;
    long max_w = 0;
    std::vector<KnapsackItem> cand;
    for (const auto& a : items) {
        assert(a.w >= 0 && a.m >= 0);
        if (a.m == 0 || a.v <= 0) continue;
        if (a.w == 0) {
            base += a.v * a.m;
            continue;
        }
        if (a.w > z) continue;
        max_w = std::max(max_w, a.w);
        cand.push_back(a);
    }
    if (cand.empty()) return base;

    std::sort(cand.begin(), cand.end(),
              [](const KnapsackItem& l, const KnapsackItem& r) { return __int128_t(l.v) * r.w < __int128_t(l.w) * r.v; });

    std::vector<std::vector<KnapsackItem>> A(max_w + 1), B(max_w + 1);
    __int128_t off = 0;
    while (!cand.empty()) {
        auto& a = cand.back();
        long n = std::min(z / a.w, a.m);
        off += __int128_t(n) * a.v;
        z -= n * a.w;
        a.m -= n;
        if (n > 0) A[a.w].push_back({a.w, a.v, n});
        if (a.m > 0) break;
        cand.pop_back();
    }
    if (cand.empty()) return base + static_cast<long>(off);
    for (const auto& a : cand)
        B[a.w].push_back(a);

    long z2 = max_w * (max_w + 1);
    constexpr __int128_t inf = __int128_t(1) << 100;
    std::vector<__int128_t> dpal(z2 + 1, -inf), dpbe(z2 + 1, 0);
    dpal[0] = 0;

    for (long w = 1; w <= max_w; w++) {
        if (A[w].empty()) continue;
        std::sort(A[w].begin(), A[w].end(), [](const KnapsackItem& l, const KnapsackItem& r) { return l.v < r.v; });
        long sz = z2 / w;
        std::vector<__int128_t> gains = {0};
        for (const auto& a : A[w]) {
            for (long i = 0; i < a.m && std::ssize(gains) <= sz; i++) {
                gains.push_back(gains.back() - a.v);
            }
        }
        internal::concave_max_plus(w, gains, dpal, -inf);
    }
    for (long w = 1; w <= max_w; w++) {
        if (B[w].empty()) continue;
        std::sort(B[w].begin(), B[w].end(), [](const KnapsackItem& l, const KnapsackItem& r) { return l.v > r.v; });
        long sz = z2 / w;
        std::vector<__int128_t> gains = {0};
        for (const auto& a : B[w]) {
            for (long i = 0; i < a.m && std::ssize(gains) <= sz; i++) {
                gains.push_back(gains.back() + a.v);
            }
        }
        internal::concave_max_plus(w, gains, dpbe, -inf);
    }

    for (long i = z2 - 1; i >= 0; i--)
        dpal[i] = std::max(dpal[i], dpal[i + 1]);
    __int128_t ans = off;
    for (long de = 0; de <= max_w * (max_w - 1); de++) {
        if (dpal[de] > -inf) {
            ans = std::max(ans, off + dpal[de] + dpbe[de + z]);
        }
    }
    return base + static_cast<long>(ans);
}

// [knapsack_11] 小価値・巨大 Z, W_i, m_i 向け個数制限付きナップサック
// 時間計算量: O(N \log N + V_{max}^3 \log V_{max}), 空間計算量: O(V_{max}^2 + N)
// 事前条件: z >= 0, w_i >= 0, m_i >= 0
inline long knapsack_bounded_small_value(long z, const std::vector<KnapsackItem>& items) {
    assert(z >= 0);
    long base = 0;
    long max_v = 0;
    std::vector<KnapsackItem> cand;
    for (const auto& a : items) {
        assert(a.w >= 0 && a.m >= 0);
        if (a.m == 0 || a.v <= 0) continue;
        if (a.w == 0) {
            base += a.v * a.m;
            continue;
        }
        if (a.w > z) continue;
        max_v = std::max(max_v, a.v);
        cand.push_back(a);
    }
    if (cand.empty()) return base;

    std::sort(cand.begin(), cand.end(),
              [](const KnapsackItem& l, const KnapsackItem& r) { return __int128_t(l.v) * r.w < __int128_t(l.w) * r.v; });

    std::vector<std::vector<KnapsackItem>> A(max_v + 1), B(max_v + 1);
    __int128_t off = 0;
    while (!cand.empty()) {
        auto& a = cand.back();
        long n = std::min(z / a.w, a.m);
        off += __int128_t(n) * a.v;
        z -= n * a.w;
        a.m -= n;
        if (n > 0) A[a.v].push_back({a.w, a.v, n});
        if (a.m > 0) break;
        cand.pop_back();
    }
    if (cand.empty()) return base + static_cast<long>(off);
    for (const auto& a : cand)
        B[a.v].push_back(a);

    long z2 = max_v * (max_v + 1);
    constexpr __int128_t inf = __int128_t(1) << 100;
    std::vector<__int128_t> dpal(z2 + 1, 0), dpbe(z2 + 1, -inf);
    dpbe[0] = 0;

    for (long v = 1; v <= max_v; v++) {
        if (A[v].empty()) continue;
        std::sort(A[v].begin(), A[v].end(), [](const KnapsackItem& l, const KnapsackItem& r) { return l.w > r.w; });
        long sz = z2 / v;
        std::vector<__int128_t> gains = {0};
        for (const auto& a : A[v]) {
            for (long i = 0; i < a.m && std::ssize(gains) <= sz; i++) {
                gains.push_back(gains.back() + a.w);
            }
        }
        internal::concave_max_plus(v, gains, dpal, -inf);
    }
    for (long v = 1; v <= max_v; v++) {
        if (B[v].empty()) continue;
        std::sort(B[v].begin(), B[v].end(), [](const KnapsackItem& l, const KnapsackItem& r) { return l.w < r.w; });
        long sz = z2 / v;
        std::vector<__int128_t> gains = {0};
        for (const auto& a : B[v]) {
            for (long i = 0; i < a.m && std::ssize(gains) <= sz; i++) {
                gains.push_back(gains.back() - a.w);
            }
        }
        internal::concave_max_plus(v, gains, dpbe, -inf);
    }

    for (long i = z2 - 1; i >= 0; i--)
        dpbe[i] = std::max(dpbe[i], dpbe[i + 1]);
    __int128_t ans = off;
    long q = z2;
    for (long de = max_v * (max_v - 1); de >= 0; de--) {
        while (dpal[de] + z < -dpbe[q])
            q--;
        ans = std::max(ans, off - de + q);
    }
    return base + static_cast<long>(ans);
}
