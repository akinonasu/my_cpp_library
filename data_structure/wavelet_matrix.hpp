#pragma once

#include <algorithm>
#include <bit>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <queue>
#include <type_traits>
#include <utility>
#include <vector>

namespace wm_detail {

template <typename T>
struct DefaultSumTypeHelper {
    using type = void;
};

template <std::integral T>
struct DefaultSumTypeHelper<T> {
    using type = long;
};

template <std::floating_point T>
struct DefaultSumTypeHelper<T> {
    using type = T;
};

template <typename T>
using DefaultSumType = typename DefaultSumTypeHelper<T>::type;

// 簡潔ビットベクトル (Succinct BitVector)
class BitVector {
   public:
    BitVector() = default;

    explicit BitVector(int n) : n_(n), bits_((n >> 6) + 1, 0ULL), sum_((n >> 6) + 2, 0) { assert(n >= 0); }

    void set(int i) {
        assert(0 <= i && i < n_);
        bits_[i >> 6] |= (1ULL << (i & 63));
    }

    void build() {
        const int blocks = static_cast<int>(bits_.size());
        sum_[0] = 0;
        for (int i = 0; i < blocks; ++i) {
            sum_[i + 1] = sum_[i] + std::popcount(bits_[i]);
        }
    }

    bool access(int i) const {
        assert(0 <= i && i < n_);
        return ((bits_[i >> 6] >> (i & 63)) & 1ULL) != 0ULL;
    }

    // [0, i) に含まれる 1 の個数
    int rank1(int i) const {
        assert(0 <= i && i <= n_);
        const int blk = i >> 6;
        const int rem = i & 63;
        const std::uint64_t mask = (1ULL << rem) - 1ULL;
        return sum_[blk] + std::popcount(bits_[blk] & mask);
    }

    // [0, i) に含まれる 0 の個数
    int rank0(int i) const {
        assert(0 <= i && i <= n_);
        return i - rank1(i);
    }

   private:
    int n_ = 0;
    std::vector<std::uint64_t> bits_;
    std::vector<int> sum_;
};

template <typename SumT>
struct SumStorage {
    std::vector<std::vector<SumT>> pref_sum;
};

template <>
struct SumStorage<void> {};

}  // namespace wm_detail

template <typename T, typename SumT = wm_detail::DefaultSumType<T>>
    requires std::totally_ordered<T> && std::default_initializable<T>
class WaveletMatrix {
   public:
    static constexpr bool has_sum = !std::is_void_v<SumT>;

    struct IntersectResult {
        T val;
        int cnt1;
        int cnt2;
    };

    WaveletMatrix() = default;

    // 値のみを指定して構築（SumT が void でない場合、値をそのまま重みとして使用）
    explicit WaveletMatrix(const std::vector<T>& a) {
        if constexpr (has_sum) {
            static_assert(std::convertible_to<T, SumT>, "T must be convertible to SumT, or pass weights explicitly / set SumT = void.");
            std::vector<SumT> w(a.size());
            for (int i = 0; i < static_cast<int>(a.size()); ++i) {
                w[i] = static_cast<SumT>(a[i]);
            }
            build_internal(a, w);
        } else {
            build_internal_no_sum(a);
        }
    }

    // 値と重みを個別に指定して構築
    WaveletMatrix(const std::vector<T>& a, const std::vector<SumT>& w)
        requires has_sum
    {
        assert(a.size() == w.size());
        build_internal(a, w);
    }

    int size() const noexcept { return n_; }
    bool empty() const noexcept { return n_ == 0; }
    int distinct_size() const noexcept { return sigma_; }

    // i 番目の要素を取得: O(log σ)
    T access(int i) const {
        assert(0 <= i && i < n_);
        int id = 0;
        for (int d = bit_len_ - 1; d >= 0; --d) {
            if (bv_[d].access(i)) {
                id |= (1 << d);
                i = mid_[d] + bv_[d].rank1(i);
            } else {
                i = bv_[d].rank0(i);
            }
        }
        return vals_[id];
    }

    T operator[](int i) const { return access(i); }

    // i 番目の重みを取得: O(1)
    SumT weight_at(int i) const
        requires has_sum
    {
        assert(0 <= i && i < n_);
        return sums_.pref_sum[bit_len_][i + 1] - sums_.pref_sum[bit_len_][i];
    }

    // [l, r) に含まれる x の個数: O(log σ)
    int rank(const T& x, int l, int r) const {
        assert(0 <= l && l <= r && r <= n_);
        if (l == r || sigma_ == 0) return 0;
        const int id = find_exact_id(x);
        if (id == -1) return 0;
        for (int d = bit_len_ - 1; d >= 0; --d) {
            const int l0 = bv_[d].rank0(l);
            const int r0 = bv_[d].rank0(r);
            if ((id >> d) & 1) {
                l = mid_[d] + (l - l0);
                r = mid_[d] + (r - r0);
            } else {
                l = l0;
                r = r0;
            }
        }
        return r - l;
    }

    // [0, r) に含まれる x の個数: O(log σ)
    int rank(const T& x, int r) const { return rank(x, 0, r); }

    // [l, r) で k 番目 (0-indexed) に小さい値: O(log σ)
    T kth_smallest(int l, int r, int k) const {
        assert(0 <= l && l <= r && r <= n_);
        assert(0 <= k && k < r - l);
        int id = 0;
        for (int d = bit_len_ - 1; d >= 0; --d) {
            const int l0 = bv_[d].rank0(l);
            const int r0 = bv_[d].rank0(r);
            const int cnt0 = r0 - l0;
            if (k < cnt0) {
                l = l0;
                r = r0;
            } else {
                id |= (1 << d);
                k -= cnt0;
                l = mid_[d] + (l - l0);
                r = mid_[d] + (r - r0);
            }
        }
        return vals_[id];
    }

    // [l, r) で k 番目 (0-indexed) に大きい値: O(log σ)
    T kth_largest(int l, int r, int k) const {
        assert(0 <= l && l <= r && r <= n_);
        assert(0 <= k && k < r - l);
        return kth_smallest(l, r, (r - l) - 1 - k);
    }

    // [l, r) で lower <= x < upper を満たす要素の個数: O(log σ)
    int range_freq(int l, int r, const T& lower, const T& upper) const {
        assert(0 <= l && l <= r && r <= n_);
        assert(!(upper < lower));
        if (l == r || lower == upper || sigma_ == 0) return 0;
        const int lid = lower_id(lower);
        const int uid = lower_id(upper);
        return count_less_id(l, r, uid) - count_less_id(l, r, lid);
    }

    int count_less(int l, int r, const T& upper) const {
        assert(0 <= l && l <= r && r <= n_);
        return count_less_id(l, r, lower_id(upper));
    }

    int count_less_eq(int l, int r, const T& upper) const {
        assert(0 <= l && l <= r && r <= n_);
        return count_less_id(l, r, upper_id(upper));
    }

    int count_greater(int l, int r, const T& lower) const {
        assert(0 <= l && l <= r && r <= n_);
        return (r - l) - count_less_id(l, r, upper_id(lower));
    }

    int count_greater_eq(int l, int r, const T& lower) const {
        assert(0 <= l && l <= r && r <= n_);
        return (r - l) - count_less_id(l, r, lower_id(lower));
    }

    // [l, r) で x < upper を満たす最大値（存在しない場合は {false, T{}}）: O(log σ)
    std::pair<bool, T> prev_value(int l, int r, const T& upper) const {
        assert(0 <= l && l <= r && r <= n_);
        const int cnt = count_less_id(l, r, lower_id(upper));
        if (cnt == 0) return {false, T{}};
        return {true, kth_smallest(l, r, cnt - 1)};
    }

    // [l, r) で x <= upper を満たす最大値（存在しない場合は {false, T{}}）: O(log σ)
    std::pair<bool, T> prev_value_le(int l, int r, const T& upper) const {
        assert(0 <= l && l <= r && r <= n_);
        const int cnt = count_less_id(l, r, upper_id(upper));
        if (cnt == 0) return {false, T{}};
        return {true, kth_smallest(l, r, cnt - 1)};
    }

    // [l, r) で lower <= x を満たす最小値（存在しない場合は {false, T{}}）: O(log σ)
    std::pair<bool, T> next_value(int l, int r, const T& lower) const {
        assert(0 <= l && l <= r && r <= n_);
        const int cnt = count_less_id(l, r, lower_id(lower));
        if (cnt == r - l) return {false, T{}};
        return {true, kth_smallest(l, r, cnt)};
    }

    // [l, r) で lower < x を満たす最小値（存在しない場合は {false, T{}}）: O(log σ)
    std::pair<bool, T> next_value_gt(int l, int r, const T& lower) const {
        assert(0 <= l && l <= r && r <= n_);
        const int cnt = count_less_id(l, r, upper_id(lower));
        if (cnt == r - l) return {false, T{}};
        return {true, kth_smallest(l, r, cnt)};
    }

    // --- 重み・和クエリ ---

    // [l, r) の重みの総和: O(1)
    SumT range_sum(int l, int r) const
        requires has_sum
    {
        assert(0 <= l && l <= r && r <= n_);
        return sums_.pref_sum[bit_len_][r] - sums_.pref_sum[bit_len_][l];
    }

    // [l, r) で lower <= x < upper を満たす要素の重み和: O(log σ)
    SumT range_sum(int l, int r, const T& lower, const T& upper) const
        requires has_sum
    {
        assert(0 <= l && l <= r && r <= n_);
        assert(!(upper < lower));
        if (l == r || lower == upper || sigma_ == 0) return SumT{};
        const int lid = lower_id(lower);
        const int uid = lower_id(upper);
        return sum_less_id(l, r, uid) - sum_less_id(l, r, lid);
    }

    // [l, r) で x < upper を満たす要素の重み和: O(log σ)
    SumT sum_less(int l, int r, const T& upper) const
        requires has_sum
    {
        assert(0 <= l && l <= r && r <= n_);
        return sum_less_id(l, r, lower_id(upper));
    }

    // [l, r) で値が小さい方から k 個の重み和 (0 <= k <= r - l): O(log σ)
    // ※同値が複数ある場合、元の配列で左にある要素が優先されます
    SumT kth_smallest_sum(int l, int r, int k) const
        requires has_sum
    {
        assert(0 <= l && l <= r && r <= n_);
        assert(0 <= k && k <= r - l);
        if (k == 0) return SumT{};
        if (k == r - l) return range_sum(l, r);
        SumT res{};
        for (int d = bit_len_ - 1; d >= 0; --d) {
            const int l0 = bv_[d].rank0(l);
            const int r0 = bv_[d].rank0(r);
            const int cnt0 = r0 - l0;
            if (k <= cnt0) {
                l = l0;
                r = r0;
            } else {
                res += sums_.pref_sum[d][r0] - sums_.pref_sum[d][l0];
                k -= cnt0;
                l = mid_[d] + (l - l0);
                r = mid_[d] + (r - r0);
            }
        }
        res += sums_.pref_sum[0][l + k] - sums_.pref_sum[0][l];
        return res;
    }

    // [l, r) で値が大きい方から k 個の重み和 (0 <= k <= r - l): O(log σ)
    // ※同値が複数ある場合、元の配列で左にある要素が優先されます
    SumT kth_largest_sum(int l, int r, int k) const
        requires has_sum
    {
        assert(0 <= l && l <= r && r <= n_);
        assert(0 <= k && k <= r - l);
        if (k == 0) return SumT{};
        if (k == r - l) return range_sum(l, r);
        SumT res{};
        for (int d = bit_len_ - 1; d >= 0; --d) {
            const int l0 = bv_[d].rank0(l);
            const int r0 = bv_[d].rank0(r);
            const int l1 = mid_[d] + (l - l0);
            const int r1 = mid_[d] + (r - r0);
            const int cnt1 = r1 - l1;
            if (k <= cnt1) {
                l = l1;
                r = r1;
            } else {
                res += sums_.pref_sum[d][r1] - sums_.pref_sum[d][l1];
                k -= cnt1;
                l = l0;
                r = r0;
            }
        }
        res += sums_.pref_sum[0][l + k] - sums_.pref_sum[0][l];
        return res;
    }

    // --- 応用列挙クエリ ---

    // [l, r) 内で lower <= x < upper を満たす {値, 個数} を値の昇順に列挙: O(K log σ)
    std::vector<std::pair<T, int>> range_list(int l, int r, const T& lower, const T& upper) const {
        assert(0 <= l && l <= r && r <= n_);
        assert(!(upper < lower));
        std::vector<std::pair<T, int>> res;
        if (l == r || lower == upper || sigma_ == 0) return res;
        const int lid = lower_id(lower);
        const int uid = lower_id(upper);
        if (lid >= uid) return res;

        auto dfs = [&](auto& self, int d, int cl, int cr, int cur_id) -> void {
            if (cl == cr) return;
            const long seg_l = cur_id;
            const long seg_r = cur_id + (1L << (d + 1));
            if (seg_r <= lid || uid <= seg_l) return;
            if (d == -1) {
                res.emplace_back(vals_[cur_id], cr - cl);
                return;
            }
            const int l0 = bv_[d].rank0(cl);
            const int r0 = bv_[d].rank0(cr);
            self(self, d - 1, l0, r0, cur_id);
            const int l1 = mid_[d] + (cl - l0);
            const int r1 = mid_[d] + (cr - r0);
            self(self, d - 1, l1, r1, cur_id | (1 << d));
        };
        dfs(dfs, bit_len_ - 1, l, r, 0);
        return res;
    }

    // [l, r) 内で min_count 回以上出現する {値, 個数} を値の昇順に列挙: O(((r - l) / min_count) log σ)
    std::vector<std::pair<T, int>> heavy_hitters(int l, int r, int min_count) const {
        assert(0 <= l && l <= r && r <= n_);
        assert(min_count >= 1);
        std::vector<std::pair<T, int>> res;
        if (r - l < min_count || sigma_ == 0) return res;

        auto dfs = [&](auto& self, int d, int cl, int cr, int cur_id) -> void {
            if (cr - cl < min_count) return;
            if (d == -1) {
                res.emplace_back(vals_[cur_id], cr - cl);
                return;
            }
            const int l0 = bv_[d].rank0(cl);
            const int r0 = bv_[d].rank0(cr);
            self(self, d - 1, l0, r0, cur_id);
            const int l1 = mid_[d] + (cl - l0);
            const int r1 = mid_[d] + (cr - r0);
            self(self, d - 1, l1, r1, cur_id | (1 << d));
        };
        dfs(dfs, bit_len_ - 1, l, r, 0);
        return res;
    }

    // [l, r) 内の出現回数上位 k 個の {値, 個数} を出現回数降順（同数は値昇順）で列挙
    std::vector<std::pair<T, int>> top_k(int l, int r, int k) const {
        assert(0 <= l && l <= r && r <= n_);
        assert(k >= 0);
        std::vector<std::pair<T, int>> res;
        if (l == r || k == 0 || sigma_ == 0) return res;

        struct Node {
            int cnt;
            int depth;
            int cl;
            int cr;
            int id;
            bool operator<(const Node& o) const {
                if (cnt != o.cnt) return cnt < o.cnt;
                if (depth != o.depth) return depth > o.depth;  // -1 (葉) を優先
                return id > o.id;                              // 小さい ID を優先
            }
        };

        std::priority_queue<Node> pq;
        pq.push({r - l, bit_len_ - 1, l, r, 0});
        while (!pq.empty() && static_cast<int>(res.size()) < k) {
            auto [cnt, d, cl, cr, cur_id] = pq.top();
            pq.pop();
            if (d == -1) {
                res.emplace_back(vals_[cur_id], cnt);
                continue;
            }
            const int l0 = bv_[d].rank0(cl);
            const int r0 = bv_[d].rank0(cr);
            if (r0 > l0) {
                pq.push({r0 - l0, d - 1, l0, r0, cur_id});
            }
            const int l1 = mid_[d] + (cl - l0);
            const int r1 = mid_[d] + (cr - r0);
            if (r1 > l1) {
                pq.push({r1 - l1, d - 1, l1, r1, cur_id | (1 << d)});
            }
        }
        return res;
    }

    // [l, r) の最頻値 {値, 個数} を返す（空区間なら {false, {}}）
    std::pair<bool, std::pair<T, int>> mode(int l, int r) const {
        assert(0 <= l && l <= r && r <= n_);
        if (l == r) return {false, {}};
        auto tk = top_k(l, r, 1);
        return {true, tk[0]};
    }

    // [l1, r1) と [l2, r2) の両方に含まれる値とそれぞれの個数を値の昇順に列挙
    std::vector<IntersectResult> intersect(int l1, int r1, int l2, int r2) const {
        assert(0 <= l1 && l1 <= r1 && r1 <= n_);
        assert(0 <= l2 && l2 <= r2 && r2 <= n_);
        std::vector<IntersectResult> res;
        if (l1 == r1 || l2 == r2 || sigma_ == 0) return res;

        auto dfs = [&](auto& self, int d, int cl1, int cr1, int cl2, int cr2, int cur_id) -> void {
            if (cl1 == cr1 || cl2 == cr2) return;
            if (d == -1) {
                res.push_back({vals_[cur_id], cr1 - cl1, cr2 - cl2});
                return;
            }
            const int l1_0 = bv_[d].rank0(cl1);
            const int r1_0 = bv_[d].rank0(cr1);
            const int l2_0 = bv_[d].rank0(cl2);
            const int r2_0 = bv_[d].rank0(cr2);
            self(self, d - 1, l1_0, r1_0, l2_0, r2_0, cur_id);

            const int l1_1 = mid_[d] + (cl1 - l1_0);
            const int r1_1 = mid_[d] + (cr1 - r1_0);
            const int l2_1 = mid_[d] + (cl2 - l2_0);
            const int r2_1 = mid_[d] + (cr2 - r2_0);
            self(self, d - 1, l1_1, r1_1, l2_1, r2_1, cur_id | (1 << d));
        };
        dfs(dfs, bit_len_ - 1, l1, r1, l2, r2, 0);
        return res;
    }

   private:
    int n_ = 0;
    int sigma_ = 0;
    int bit_len_ = 0;
    std::vector<T> vals_;
    std::vector<wm_detail::BitVector> bv_;
    std::vector<int> mid_;
    wm_detail::SumStorage<SumT> sums_;

    void init_compression(const std::vector<T>& a, std::vector<int>& ids) {
        n_ = static_cast<int>(a.size());
        if (n_ == 0) {
            sigma_ = 0;
            bit_len_ = 0;
            return;
        }
        vals_ = a;
        std::ranges::sort(vals_);
        auto ret = std::ranges::unique(vals_);
        vals_.erase(ret.begin(), ret.end());
        sigma_ = static_cast<int>(vals_.size());
        bit_len_ = std::max(1, static_cast<int>(std::bit_width(static_cast<unsigned int>(sigma_ - 1))));

        ids.resize(n_);
        for (int i = 0; i < n_; ++i) {
            ids[i] = static_cast<int>(std::ranges::lower_bound(vals_, a[i]) - vals_.begin());
        }
        bv_.assign(bit_len_, wm_detail::BitVector(n_));
        mid_.assign(bit_len_, 0);
    }

    void build_internal_no_sum(const std::vector<T>& a) {
        std::vector<int> cur_id;
        init_compression(a, cur_id);
        if (n_ == 0) return;

        std::vector<int> nxt_id(n_);
        for (int d = bit_len_ - 1; d >= 0; --d) {
            int cnt0 = 0;
            for (int i = 0; i < n_; ++i) {
                if ((cur_id[i] >> d) & 1) {
                    bv_[d].set(i);
                } else {
                    ++cnt0;
                }
            }
            bv_[d].build();
            mid_[d] = cnt0;

            int p0 = 0, p1 = cnt0;
            for (int i = 0; i < n_; ++i) {
                if ((cur_id[i] >> d) & 1) {
                    nxt_id[p1++] = cur_id[i];
                } else {
                    nxt_id[p0++] = cur_id[i];
                }
            }
            std::swap(cur_id, nxt_id);
        }
    }

    void build_internal(const std::vector<T>& a, const std::vector<SumT>& w)
        requires has_sum
    {
        std::vector<int> cur_id;
        init_compression(a, cur_id);
        sums_.pref_sum.assign(bit_len_ + 1, std::vector<SumT>(n_ + 1, SumT{}));
        for (int i = 0; i < n_; ++i) {
            sums_.pref_sum[bit_len_][i + 1] = sums_.pref_sum[bit_len_][i] + w[i];
        }
        if (n_ == 0) return;

        std::vector<int> nxt_id(n_);
        std::vector<SumT> cur_w = w;
        std::vector<SumT> nxt_w(n_);

        for (int d = bit_len_ - 1; d >= 0; --d) {
            int cnt0 = 0;
            for (int i = 0; i < n_; ++i) {
                if ((cur_id[i] >> d) & 1) {
                    bv_[d].set(i);
                } else {
                    ++cnt0;
                }
            }
            bv_[d].build();
            mid_[d] = cnt0;

            int p0 = 0, p1 = cnt0;
            for (int i = 0; i < n_; ++i) {
                if ((cur_id[i] >> d) & 1) {
                    nxt_id[p1] = cur_id[i];
                    nxt_w[p1] = cur_w[i];
                    ++p1;
                } else {
                    nxt_id[p0] = cur_id[i];
                    nxt_w[p0] = cur_w[i];
                    ++p0;
                }
            }
            for (int i = 0; i < n_; ++i) {
                sums_.pref_sum[d][i + 1] = sums_.pref_sum[d][i] + nxt_w[i];
            }
            std::swap(cur_id, nxt_id);
            std::swap(cur_w, nxt_w);
        }
    }

    int lower_id(const T& x) const { return static_cast<int>(std::ranges::lower_bound(vals_, x) - vals_.begin()); }

    int upper_id(const T& x) const { return static_cast<int>(std::ranges::upper_bound(vals_, x) - vals_.begin()); }

    int find_exact_id(const T& x) const {
        const auto it = std::ranges::lower_bound(vals_, x);
        if (it == vals_.end() || *it != x) return -1;
        return static_cast<int>(it - vals_.begin());
    }

    int count_less_id(int l, int r, int uid) const {
        if (uid <= 0 || l == r) return 0;
        if (uid >= sigma_) return r - l;
        int cnt = 0;
        for (int d = bit_len_ - 1; d >= 0; --d) {
            const int l0 = bv_[d].rank0(l);
            const int r0 = bv_[d].rank0(r);
            if ((uid >> d) & 1) {
                cnt += r0 - l0;
                l = mid_[d] + (l - l0);
                r = mid_[d] + (r - r0);
            } else {
                l = l0;
                r = r0;
            }
        }
        return cnt;
    }

    SumT sum_less_id(int l, int r, int uid) const
        requires has_sum
    {
        if (uid <= 0 || l == r) return SumT{};
        if (uid >= sigma_) return range_sum(l, r);
        SumT res{};
        for (int d = bit_len_ - 1; d >= 0; --d) {
            const int l0 = bv_[d].rank0(l);
            const int r0 = bv_[d].rank0(r);
            if ((uid >> d) & 1) {
                res += sums_.pref_sum[d][r0] - sums_.pref_sum[d][l0];
                l = mid_[d] + (l - l0);
                r = mid_[d] + (r - r0);
            } else {
                l = l0;
                r = r0;
            }
        }
        return res;
    }
};

// CTAD (クラステンプレート引数推論) ガイド
template <typename T>
WaveletMatrix(std::vector<T>) -> WaveletMatrix<T>;

template <typename T, typename SumT>
WaveletMatrix(std::vector<T>, std::vector<SumT>) -> WaveletMatrix<T, SumT>;
