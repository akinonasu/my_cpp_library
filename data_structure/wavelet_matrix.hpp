#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include "components/quad_bit_vector.hpp"

template <class T = std::int64_t,
          class Sum = std::conditional_t<std::is_integral_v<T>, std::conditional_t<(sizeof(T) < 8), std::int64_t, __int128_t>, T>>
class WaveletMatrix {
    using BitVector = QuadBitVector;
    int n = 0, levels = 0, ordered = 0;
    std::vector<T> values;
    std::vector<std::uint32_t> original;
    std::vector<int> positions, offsets;
    bool dense = false;
    std::vector<BitVector> bits;
    std::vector<std::array<int, 4>> middle;
    static constexpr bool can_narrow = std::is_integral_v<T> && sizeof(T) <= 8 && std::is_same_v<Sum, __int128_t>;
    using CompactSum = std::conditional_t<can_narrow, std::int64_t, Sum>;
    mutable std::vector<std::vector<Sum>> sum_prefix;
    mutable std::vector<std::vector<CompactSum>> compact_prefix;
    mutable bool sum_ready = false, narrow_sums = false;
    static constexpr bool supports_sum = requires(const T& x, Sum s) {
        Sum(x);
        s + s;
        s - s;
        s * 1;
    };
    std::array<std::uint32_t, 8> small_codes(int l, int r) const {
        std::array<std::uint32_t, 8> a;
        a.fill(std::numeric_limits<std::uint32_t>::max());
        for (int i = 0; i < r - l; ++i)
            a[i] = original[l + i];
        auto swap = [&](int i, int j) {
            auto lo = std::min(a[i], a[j]), hi = std::max(a[i], a[j]);
            a[i] = lo;
            a[j] = hi;
        };
        swap(0, 1), swap(2, 3);
        if (r - l <= 4) {
            swap(0, 2), swap(1, 3), swap(1, 2);
        } else {
            swap(4, 5), swap(6, 7);
            swap(0, 2), swap(1, 3), swap(4, 6), swap(5, 7);
            swap(1, 2), swap(5, 6), swap(0, 4), swap(3, 7);
            swap(1, 5), swap(2, 6);
            swap(1, 4), swap(3, 6);
            swap(2, 4), swap(3, 5);
            swap(3, 4);
        }
        return a;
    }
    void check(int l, int r) const { assert(0 <= l && l <= r && r <= n); }
    std::uint32_t code(const T& x) const {
        if (values.empty()) return 0;
        if constexpr (std::is_integral_v<T> && sizeof(T) <= 8)
            if (dense) {
                if (x <= values.front()) return 0;
                if (values.back() < x) return static_cast<int>(values.size());
                return std::uint64_t(x) - std::uint64_t(values.front());
            }
        int p = 0, length = static_cast<int>(values.size());
        while (length > 1) {
            int half = length / 2;
            p += (values[p + half] < x) * half;
            length -= half;
        }
        return p + (values[p] < x);
    }
    std::pair<int, int> interval(int l, int r, std::uint32_t x) const {
        for (int d = 0; d < levels; ++d) {
            int digit = (x >> (2 * (levels - 1 - d))) & 3;
            l = middle[d][digit] + bits[d].rank(l, digit);
            r = middle[d][digit] + bits[d].rank(r, digit);
        }
        return {l, r};
    }
    int less_than(int l, int r, std::uint32_t x) const {
        if (!x || l == r) return 0;
        if (x >= values.size()) return r - l;
        if (ordered) return std::clamp(ordered > 0 ? offsets[x] - l : r - n + offsets[x], 0, r - l);
        int answer = 0;
        if (r - l <= 8) {
            for (int i = l; i < r; ++i)
                answer += original[i] < x;
            return answer;
        }
        for (int d = 0; d < levels; ++d) {
            int digit = (x >> (2 * (levels - 1 - d))) & 3;
            int a = bits[d].rank(l, digit), b = bits[d].rank(r, digit);
            if (digit == 3)
                answer += r - l - b + a;
            else
                answer += bits[d].below(r, digit) - bits[d].below(l, digit);
            l = middle[d][digit] + a;
            r = middle[d][digit] + b;
        }
        return answer;
    }
    std::vector<std::uint32_t> compress(const std::vector<T>& a) {
        std::vector<std::uint32_t> current(n);
        bool ascending = std::is_sorted(a.begin(), a.end());
        if (ascending || std::is_sorted(a.rbegin(), a.rend())) {
            ordered = ascending ? 1 : -1;
            values.clear();
            values.reserve(n);
            for (int i = 0; i < n; ++i) {
                int j = ascending ? i : n - 1 - i;
                if (values.empty() || values.back() != a[j]) values.push_back(a[j]);
                current[j] = static_cast<int>(values.size()) - 1;
            }
            values.shrink_to_fit();
            return current;
        }
        if constexpr (std::is_integral_v<T> && !std::is_same_v<T, bool> && sizeof(T) <= 8) {
            if (n >= 1024) {
                using U = std::make_unsigned_t<T>;
                struct Item {
                    U key;
                    std::uint32_t index;
                };
                std::vector<Item> order(n);
                U lo = std::numeric_limits<U>::max(), hi = 0;
                constexpr U sign = std::is_signed_v<T> ? U(1) << (sizeof(T) * 8 - 1) : U(0);
                for (int i = 0; i < n; ++i) {
                    U key = U(a[i]) ^ sign;
                    order[i] = {key, std::uint32_t(i)};
                    lo = std::min(lo, key);
                    hi = std::max(hi, key);
                }
                if (hi - lo < 65536 && std::uint64_t(hi - lo) < unsigned(n)) {
                    std::vector<int> codes(std::size_t(hi - lo) + 1, -1);
                    for (const auto& v : order)
                        codes[v.key - lo] = 0;
                    values.clear();
                    values.reserve(codes.size());
                    for (std::size_t i = 0; i < codes.size(); ++i)
                        if (codes[i] == 0) {
                            codes[i] = static_cast<int>(values.size());
                            values.push_back(T(U(lo + U(i)) ^ sign));
                        }
                    for (const auto& v : order)
                        current[v.index] = codes[v.key - lo];
                    values.shrink_to_fit();
                    return current;
                }
                std::vector<Item> next(n);
                for (int shift = 0, width = std::bit_width(U(lo ^ hi)); shift < width; shift += 11) {
                    std::array<int, 2048> offsets{};
                    for (const auto& v : order)
                        ++offsets[(v.key >> shift) & 2047];
                    int total = 0;
                    for (int& c : offsets) {
                        int t = c;
                        c = total;
                        total += t;
                    }
                    for (const auto& v : order)
                        next[offsets[(v.key >> shift) & 2047]++] = v;
                    order.swap(next);
                }
                values.clear();
                values.reserve(n);
                U previous = 0;
                for (const auto& v : order) {
                    if (values.empty() || v.key != previous) values.push_back(a[v.index]);
                    previous = v.key;
                    current[v.index] = static_cast<int>(values.size()) - 1;
                }
                values.shrink_to_fit();
                return current;
            }
        }
        values = a;
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()), values.end());
        values.shrink_to_fit();
        for (int i = 0; i < n; ++i)
            current[i] = code(a[i]);
        return current;
    }
    void partition(std::vector<std::uint32_t>& current, std::vector<std::uint32_t>& next, int d) const {
        auto offset = middle[d];
        int shift = 2 * (levels - 1 - d);
        for (int i = 0; i < n; ++i)
            next[offset[(current[i] >> shift) & 3]++] = current[i];
        current.swap(next);
    }
    template <class S>
    void prefix_sum(std::vector<S>& prefix, const std::vector<std::uint32_t>& current, int l, int r) const
        requires supports_sum
    {
        prefix.resize(std::size_t(r - l) + 1);
        for (int i = 0; i < r - l; ++i)
            prefix[i + 1] = prefix[i] + S(values[current[l + i]]);
    }
    template <class S>
    void partition_sum(std::vector<std::vector<S>>& prefix, std::vector<std::uint32_t>& current, std::vector<std::uint32_t>& next,
                       int d) const
        requires supports_sum
    {
        int zero = 0, one = middle[d][2], shift = 2 * (levels - 1 - d);
        for (std::uint32_t x : current) {
            bool high = (x >> (shift + 1)) & 1;
            int p = high ? one : zero;
            one += high;
            zero += !high;
            next[p] = x;
        }
        int base = 3 * d + 1;
        prefix_sum(prefix[base], next, 0, middle[d][2]);
        current.swap(next);
        partition(current, next, d);
        prefix_sum(prefix[base + 1], current, 0, middle[d][1]);
        prefix_sum(prefix[base + 2], current, middle[d][2], middle[d][3]);
    }
    struct Route {
        int digit, a, b, k, high_l, high_r, zero_l, zero_r;
    };
    [[gnu::always_inline]] Route route(int l, int r, int k, int d) const {
        const auto &left = bits[d].blocks[l >> 6], &right = bits[d].blocks[r >> 6];
        std::uint64_t lm = (std::uint64_t(1) << (l & 63)) - 1, rm = (std::uint64_t(1) << (r & 63)) - 1;
        std::uint64_t lh = left.high & lm, rh = right.high & rm;
        int ha = left.prefix[2] + left.prefix[3] + std::popcount(lh);
        int hb = right.prefix[2] + right.prefix[3] + std::popcount(rh);
        int high_l = l - ha, high_r = r - hb;
        bool high = k >= high_r - high_l;
        k -= high ? high_r - high_l : 0;
        int a = left.prefix[2 * high] + std::popcount(~left.low & (high ? lh : lm ^ lh));
        int b = right.prefix[2 * high] + std::popcount(~right.low & (high ? rh : rm ^ rh));
        bool low = k >= b - a;
        k -= low ? b - a : 0;
        return {2 * high + low, low ? (high ? ha : high_l) - a : a, low ? (high ? hb : high_r) - b : b, k, high_l, high_r, a, b};
    }
    void init_sum() const
        requires supports_sum
    {
        if constexpr (can_narrow) {
            Sum positive = 0, negative = 0;
            for (int i = 0; i < static_cast<int>(values.size()); ++i) {
                Sum v = Sum(values[i]) * (offsets[i + 1] - offsets[i]);
                if (v > 0)
                    positive += v;
                else
                    negative += v;
            }
            narrow_sums = positive <= std::numeric_limits<std::int64_t>::max() && negative >= std::numeric_limits<std::int64_t>::min();
        }
        if (narrow_sums) {
            compact_prefix.resize(ordered ? 1 : 3 * levels + 1);
            prefix_sum(compact_prefix[0], original, 0, n);
        } else {
            sum_prefix.resize(ordered ? 1 : 3 * levels + 1);
            prefix_sum(sum_prefix[0], original, 0, n);
        }
    }
    void construct(const std::vector<T>& a, bool with_sum) {
        assert(a.size() <= std::size_t(std::numeric_limits<int>::max()));
        sum_prefix.clear();
        compact_prefix.clear();
        sum_ready = narrow_sums = false;
        n = static_cast<int>(a.size());
        dense = false;
        ordered = 0;
        if (!n || std::all_of(a.begin() + 1, a.end(), [&](const T& x) { return x == a.front(); })) {
            values.clear();
            if (n) values.push_back(a.front());
            values.shrink_to_fit();
            original.clear();
            positions.clear();
            offsets = {0, n};
            levels = 0;
            bits.clear();
            middle.clear();
            sum_ready = with_sum;
            return;
        }
        original = compress(a);
        if constexpr (std::is_integral_v<T> && sizeof(T) <= 8)
            dense = std::uint64_t(values.back()) - std::uint64_t(values.front()) == values.size() - 1;
        offsets.assign(values.size() + 1, 0);
        levels = (std::bit_width(unsigned(static_cast<int>(values.size()) - 1)) + 1) / 2;
        for (std::uint32_t x : original)
            ++offsets[x + 1];
        for (int i = 0; i < static_cast<int>(values.size()); ++i)
            offsets[i + 1] += offsets[i];
        if (ordered) {
            positions.clear();
            bits.clear();
            middle.clear();
            if constexpr (supports_sum)
                if (with_sum) init_sum();
            sum_ready = with_sum;
            return;
        }
        positions.resize(n);
        auto cursor = offsets;
        for (int i = 0; i < n; ++i)
            positions[cursor[original[i]]++] = i;
        auto current = original;
        bits.clear();
        bits.reserve(levels);
        middle.resize(levels);
        if constexpr (supports_sum)
            if (with_sum) init_sum();
        std::vector<std::uint32_t> next(n);
        for (int d = 0; d < levels; ++d) {
            bits.emplace_back(n);
            bits.back().build(current, 2 * (levels - 1 - d));
            auto totals = bits.back().ranks(n);
            middle[d] = {0, totals[0], totals[0] + totals[1], totals[0] + totals[1] + totals[2]};
            if constexpr (supports_sum) {
                if (with_sum) {
                    if (narrow_sums)
                        partition_sum(compact_prefix, current, next, d);
                    else
                        partition_sum(sum_prefix, current, next, d);
                } else if (d + 1 < levels)
                    partition(current, next, d);
            } else if (d + 1 < levels)
                partition(current, next, d);
        }
        sum_ready = with_sum;
    }
    template <class S>
    S sum_smallest_impl(const std::vector<std::vector<S>>& prefix, int l, int r, int k) const
        requires supports_sum
    {
        S answer = 0;
        std::uint32_t x = 0;
        for (int d = 0; d < levels; ++d) {
            auto step = route(l, r, k, d);
            int base = 3 * d + 1;
            if (step.digit & 2) answer += prefix[base][step.high_r] - prefix[base][step.high_l];
            if (step.digit & 1) {
                const auto& row = prefix[base + 1 + (step.digit >> 1)];
                answer += row[step.zero_r] - row[step.zero_l];
            }
            k = step.k;
            if (!k) return answer;
            x = (x << 2) | step.digit;
            l = middle[d][step.digit] + step.a;
            r = middle[d][step.digit] + step.b;
        }
        return answer + S(values[x]) * k;
    }

   public:
    using sum_type = Sum;
    WaveletMatrix() = default;
    explicit WaveletMatrix(const std::vector<T>& a) { build(a); }
    WaveletMatrix(const std::vector<T>& a, bool with_sum)
        requires supports_sum
    {
        construct(a, with_sum);
    }
    void build(const std::vector<T>& a) { construct(a, false); }
    int size() const { return n; }
    bool empty() const { return n == 0; }
    T access(int i) const {
        assert(0 <= i && i < n);
        return values[levels ? original[i] : 0];
    }
    T operator[](int i) const { return access(i); }
    int count(int l, int r, const T& x) const {
        check(l, r);
        if (l == r) return 0;
        std::uint32_t c = code(x);
        if (c == values.size() || values[c] != x) return 0;
        if (!levels) return r - l;
        int start = offsets[c], end = offsets[c + 1];
        if (ordered) {
            int first = ordered > 0 ? start : n - end;
            return std::max(0, std::min(r, first + (end - start)) - std::max(l, first));
        }
        if (l == 0 && r == n) return end - start;
        if (end - start == 1) return int(l <= positions[start]) & int(positions[start] < r);
        if (end - start <= 16) {
            int answer = 0;
            for (int i = start; i < end; ++i)
                answer += unsigned(positions[i] - l) < unsigned(r - l);
            return answer;
        }
        if (r - l <= 8) {
            int answer = 0;
            for (int i = l; i < r; ++i)
                answer += original[i] == c;
            return answer;
        }
        auto [a, b] = interval(l, r, c);
        return b - a;
    }
    int rank(const T& x, int r) const { return count(0, r, x); }
    int select(const T& x, int k) const {
        std::uint32_t c = code(x);
        if (k < 0 || c == values.size() || values[c] != x || k >= offsets[c + 1] - offsets[c]) return -1;
        if (ordered) return (ordered > 0 ? offsets[c] : n - offsets[c + 1]) + k;
        return levels ? positions[offsets[c] + k] : k;
    }
    T kth_smallest(int l, int r, int k) const {
        check(l, r);
        assert(0 <= k && k < r - l);
        if (!levels) return values[0];
        if (ordered) return values[original[ordered > 0 ? l + k : r - 1 - k]];
        if (r - l <= 8) return values[small_codes(l, r)[k]];
        std::uint32_t x = 0;
        for (int d = 0; d < levels; ++d) {
            auto step = route(l, r, k, d);
            x = (x << 2) | step.digit;
            k = step.k;
            l = middle[d][step.digit] + step.a;
            r = middle[d][step.digit] + step.b;
        }
        return values[x];
    }
    T kth_largest(int l, int r, int k) const {
        check(l, r);
        assert(0 <= k && k < r - l);
        return kth_smallest(l, r, r - l - 1 - k);
    }
    int range_freq(int l, int r, const T& upper) const {
        check(l, r);
        return less_than(l, r, code(upper));
    }
    int range_freq(int l, int r, const T& lower, const T& upper) const {
        check(l, r);
        if (!(lower < upper)) return 0;
        return less_than(l, r, code(upper)) - less_than(l, r, code(lower));
    }
    void prepare_sum() const
        requires supports_sum
    {
        if (sum_ready) return;
        if (!levels) {
            sum_ready = true;
            return;
        }
        init_sum();
        if (ordered) {
            sum_ready = true;
            return;
        }
        std::vector<std::uint32_t> current = original, next(n);
        if (narrow_sums)
            for (int d = 0; d < levels; ++d)
                partition_sum(compact_prefix, current, next, d);
        else
            for (int d = 0; d < levels; ++d)
                partition_sum(sum_prefix, current, next, d);
        sum_ready = true;
    }
    Sum range_sum(int l, int r) const
        requires supports_sum
    {
        check(l, r);
        if (l == r) return Sum(0);
        if (!levels) return Sum(values[0]) * (r - l);
        if (!sum_ready && r - l <= 8) {
            Sum answer = 0;
            for (int i = l; i < r; ++i)
                answer += Sum(values[original[i]]);
            return answer;
        }
        prepare_sum();
        if (narrow_sums) return compact_prefix[0][r] - compact_prefix[0][l];
        return sum_prefix[0][r] - sum_prefix[0][l];
    }
    Sum sum_smallest(int l, int r, int k) const
        requires supports_sum
    {
        check(l, r);
        assert(0 <= k && k <= r - l);
        if (!k) return Sum(0);
        if (k == r - l) return range_sum(l, r);
        if (!levels) return Sum(values[0]) * k;
        if (ordered) return ordered > 0 ? range_sum(l, l + k) : range_sum(r - k, r);
        if (k == 1) return Sum(kth_smallest(l, r, 0));
        if (r - l <= 8) {
            auto a = small_codes(l, r);
            Sum answer = 0;
            for (int i = 0; i < k; ++i)
                answer += Sum(values[a[i]]);
            return answer;
        }
        prepare_sum();
        if (narrow_sums) return sum_smallest_impl(compact_prefix, l, r, k);
        return sum_smallest_impl(sum_prefix, l, r, k);
    }
    Sum sum_largest(int l, int r, int k) const
        requires supports_sum
    {
        check(l, r);
        assert(0 <= k && k <= r - l);
        if (!k) return Sum(0);
        if (!levels) return Sum(values[0]) * k;
        if (ordered) return ordered > 0 ? range_sum(r - k, r) : range_sum(l, l + k);
        if (k == 1) return Sum(kth_largest(l, r, 0));
        if (k == r - l) return range_sum(l, r);
        if (r - l <= 8) {
            auto a = small_codes(l, r);
            Sum answer = 0;
            for (int i = r - l - k; i < r - l; ++i)
                answer += Sum(values[a[i]]);
            return answer;
        }
        return range_sum(l, r) - sum_smallest(l, r, r - l - k);
    }
    std::optional<T> prev_value(int l, int r, const T& upper) const {
        int k = range_freq(l, r, upper);
        if (!k) return std::nullopt;
        return kth_smallest(l, r, k - 1);
    }
    std::optional<T> next_value(int l, int r, const T& lower) const {
        int k = range_freq(l, r, lower);
        if (k == r - l) return std::nullopt;
        return kth_smallest(l, r, k);
    }
};

template <class T>
WaveletMatrix(const std::vector<T>&, bool) -> WaveletMatrix<T>;
