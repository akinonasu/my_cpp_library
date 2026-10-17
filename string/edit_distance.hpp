#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string_view>
#include <utility>
#include <vector>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

namespace edit_distance_detail {
using Word = std::uint64_t;

#if defined(__x86_64__)
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize("O3")
#endif
__attribute__((target("avx2"))) inline void avx2(std::string_view s, const std::array<const Word*, 256>& masks, Word* __restrict positive,
                                                 Word* __restrict negative, int words) {
    const __m256i ones = _mm256_set1_epi64x(-1), one = _mm256_set1_epi64x(1), positions = _mm256_setr_epi64x(0, 1, 2, 3);
    for (unsigned char c : s) {
        const Word* __restrict mask = masks[c];
        unsigned carry = 0;
        __m256i previous_p = _mm256_set1_epi64x(std::numeric_limits<long long>::min()), previous_m = _mm256_setzero_si256();
        for (int j = 0; j < words; j += 4) {
            __m256i eq = _mm256_load_si256(reinterpret_cast<const __m256i*>(mask + j));
            __m256i p = _mm256_load_si256(reinterpret_cast<const __m256i*>(positive + j));
            __m256i m = _mm256_load_si256(reinterpret_cast<const __m256i*>(negative + j));
            __m256i x = _mm256_and_si256(eq, p), sum = _mm256_add_epi64(x, p);
            unsigned g = unsigned(_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_or_si256(x, _mm256_andnot_si256(sum, p)))));
            unsigned propagation = unsigned(_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpeq_epi64(sum, ones))));
            unsigned incoming = (propagation + (g << 1) + carry) ^ propagation;
            carry = incoming >> 4;
            __m256i increments = _mm256_and_si256(_mm256_srlv_epi64(_mm256_set1_epi64x(incoming), positions), one);
            __m256i xv = _mm256_or_si256(eq, m);
            __m256i xh = _mm256_or_si256(_mm256_xor_si256(_mm256_add_epi64(sum, increments), p), eq);
            __m256i ph = _mm256_or_si256(m, _mm256_andnot_si256(_mm256_or_si256(xh, p), ones)), mh = _mm256_and_si256(p, xh);
            __m256i in_p = _mm256_blend_epi32(_mm256_permute4x64_epi64(ph, 0x90), _mm256_permute4x64_epi64(previous_p, 0xff), 3);
            __m256i in_m = _mm256_blend_epi32(_mm256_permute4x64_epi64(mh, 0x90), _mm256_permute4x64_epi64(previous_m, 0xff), 3);
            previous_p = ph;
            previous_m = mh;
            ph = _mm256_or_si256(_mm256_slli_epi64(ph, 1), _mm256_srli_epi64(in_p, 63));
            mh = _mm256_or_si256(_mm256_slli_epi64(mh, 1), _mm256_srli_epi64(in_m, 63));
            _mm256_store_si256(reinterpret_cast<__m256i*>(positive + j),
                               _mm256_or_si256(mh, _mm256_andnot_si256(_mm256_or_si256(xv, ph), ones)));
            _mm256_store_si256(reinterpret_cast<__m256i*>(negative + j), _mm256_and_si256(ph, xv));
        }
    }
}

__attribute__((target("avx512f"))) inline void avx512(std::string_view s, const std::array<const Word*, 256>& masks,
                                                      Word* __restrict positive, Word* __restrict negative, int words) {
    const __m512i ones = _mm512_set1_epi64(-1), one = _mm512_set1_epi64(1);
    for (unsigned char c : s) {
        const Word* __restrict mask = masks[c];
        unsigned carry = 0;
        __m512i previous_p = _mm512_set1_epi64(std::numeric_limits<long long>::min()), previous_m = _mm512_setzero_si512();
        for (int j = 0; j < words; j += 8) {
            __m512i eq = _mm512_load_si512(mask + j), p = _mm512_load_si512(positive + j), m = _mm512_load_si512(negative + j);
            __m512i x = _mm512_and_si512(eq, p), sum = _mm512_add_epi64(x, p);
            unsigned g = unsigned(_mm512_cmp_epu64_mask(sum, x, _MM_CMPINT_LT));
            unsigned propagation = unsigned(_mm512_cmpeq_epi64_mask(sum, ones));
            unsigned incoming = (propagation + (g << 1) + carry) ^ propagation;
            carry = incoming >> 8;
            sum = _mm512_mask_add_epi64(sum, __mmask8(incoming), sum, one);
            __m512i xv = _mm512_or_si512(eq, m), xh = _mm512_ternarylogic_epi64(sum, p, eq, 0xbe);
            __m512i ph = _mm512_ternarylogic_epi64(m, xh, p, 0xf1), mh = _mm512_and_si512(p, xh);
            __m512i in_p = _mm512_alignr_epi64(ph, previous_p, 7), in_m = _mm512_alignr_epi64(mh, previous_m, 7);
            previous_p = ph;
            previous_m = mh;
            ph = _mm512_or_si512(_mm512_slli_epi64(ph, 1), _mm512_srli_epi64(in_p, 63));
            mh = _mm512_or_si512(_mm512_slli_epi64(mh, 1), _mm512_srli_epi64(in_m, 63));
            _mm512_store_si512(positive + j, _mm512_ternarylogic_epi64(mh, xv, ph, 0xf1));
            _mm512_store_si512(negative + j, _mm512_and_si512(ph, xv));
        }
    }
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif
#endif

inline int count(Word* positive, Word* negative, int n, int m) {
    int words = (m - 1) / 64 + 1;
    if (m % 64) {
        Word valid = (Word(1) << (m % 64)) - 1;
        positive[words - 1] &= valid;
        negative[words - 1] &= valid;
    }
    int result = n;
    for (int j = 0; j < words; ++j)
        result += std::popcount(positive[j]) - std::popcount(negative[j]);
    return result;
}

template <int backend = 0>
inline int bit_parallel(std::string_view s, std::string_view t) {
    if (t.empty()) return int(s.size());
    if (t.size() <= 64) {
        std::array<Word, 256> masks{};
        for (int j = 0; j < int(t.size()); ++j)
            masks[static_cast<unsigned char>(t[j])] |= Word(1) << j;
        Word p = ~Word(0), m = 0;
        for (unsigned char c : s) {
            Word eq = masks[c], xv = eq | m, xh = (((eq & p) + p) ^ p) | eq;
            Word ph = (m | ~(xh | p)) << 1 | 1, mh = (p & xh) << 1;
            p = mh | ~(xv | ph);
            m = ph & xv;
        }
        Word valid = ~Word(0) >> (64 - t.size());
        return int(s.size()) + (std::popcount(p & valid) - std::popcount(m & valid));
    }
    const int words = int((t.size() + 63) / 64);
    std::array<int, 256> ids;
    ids.fill(0);
    int alphabet = 1;
    for (unsigned char c : t)
        if (!ids[c]) ids[c] = alphabet++;
#if defined(__x86_64__)
    int kernel = backend;
    if (backend == 0 && words >= 8) {
        if (__builtin_cpu_supports("avx512f"))
            kernel = 3;
        else if (__builtin_cpu_supports("avx2"))
            kernel = 2;
    }
    if (kernel >= 2) {
        int stride = (words + 7) & ~7;
        std::vector<Word> storage(std::size_t(alphabet + 2) * stride + 7);
        Word* data = storage.data() + ((-reinterpret_cast<std::uintptr_t>(storage.data()) & 63) / sizeof(Word));
        Word *positive = data + std::size_t(alphabet) * stride, *negative = positive + stride;
        std::fill_n(positive, stride, ~Word(0));
        std::array<const Word*, 256> masks;
        for (int c = 0; c < 256; ++c)
            masks[c] = data + std::size_t(ids[c]) * stride;
        for (int j = 0, m = int(t.size()); j < m; ++j)
            data[std::size_t(ids[static_cast<unsigned char>(t[j])]) * stride + j / 64] |= Word(1) << (j % 64);
        if (kernel == 3)
            avx512(s, masks, positive, negative, stride);
        else
            avx2(s, masks, positive, negative, stride);
        return count(positive, negative, int(s.size()), int(t.size()));
    }
#endif
    std::vector<Word> masks(std::size_t(alphabet) * words), positive(words, ~Word(0)), negative(words);
    for (int j = 0, m = int(t.size()); j < m; ++j)
        masks[std::size_t(ids[static_cast<unsigned char>(t[j])]) * words + j / 64] |= Word(1) << (j % 64);
    for (unsigned char c : s) {
        const Word* mask = masks.data() + std::size_t(ids[c]) * words;
        Word carry_p = 1, carry_m = 0, carry_add = 0;
        for (int j = 0; j < words; ++j) {
            Word p = positive[j], m = negative[j], eq = mask[j], xv = eq | m;
            __uint128_t sum = __uint128_t(eq & p) + p + carry_add;
            carry_add = Word(sum >> 64);
            Word xh = (Word(sum) ^ p) | eq;
            Word ph = m | ~(xh | p), mh = p & xh;
            Word next_p = ph >> 63, next_m = mh >> 63;
            ph = (ph << 1) | carry_p;
            mh = (mh << 1) | carry_m;
            positive[j] = mh | ~(xv | ph);
            negative[j] = ph & xv;
            carry_p = next_p;
            carry_m = next_m;
        }
    }
    return count(positive.data(), negative.data(), int(s.size()), int(t.size()));
}

inline int extend(std::string_view s, std::string_view t, int x, int y, int& budget) {
    int limit = std::min(int(s.size()) - x, int(t.size()) - y);
    while (limit >= 8) {
        if (--budget < 0) return -1;
        Word a, b;
        std::memcpy(&a, s.data() + x, 8);
        std::memcpy(&b, t.data() + y, 8);
        Word diff = a ^ b;
        if (diff) {
            if constexpr (std::endian::native == std::endian::little)
                return x + std::countr_zero(diff) / 8;
            else
                return x + std::countl_zero(diff) / 8;
        }
        x += 8;
        y += 8;
        limit -= 8;
    }
    while (limit--) {
        if (--budget < 0) return -1;
        if (s[x] != t[y]) break;
        ++x;
        ++y;
    }
    return x;
}

inline int wavefront(std::string_view s, std::string_view t, int cap, int budget) {
    int n = int(s.size()), m = int(t.size()), base = cap + 1;
    std::vector<int> previous(2 * cap + 3, -1), current(2 * cap + 3, -1);
    previous[base] = extend(s, t, 0, 0, budget);
    if (previous[base] < 0) return -1;
    if (previous[base] == n && n == m) return 0;
    for (int d = 1; d <= cap; ++d) {
        for (int k = std::max(-d, -m), last = std::min(d, n); k <= last; ++k) {
            if (--budget < 0) return -1;
            int x = -1, a = previous[base + k - 1], b = previous[base + k], c = previous[base + k + 1];
            if (a >= 0 && a < n) x = a + 1;
            if (b >= 0 && b < n && b - k < m) x = std::max(x, b + 1);
            if (c >= 0 && c - (k + 1) < m) x = std::max(x, c);
            if (x >= 0) {
                x = extend(s, t, x, x - k, budget);
                if (x < 0) return -1;
                if (x == n && x - k == m) return d;
            }
            current[base + k] = x;
        }
        previous.swap(current);
    }
    return -1;
}
}  // namespace edit_distance_detail

inline int edit_distance(std::string_view s, std::string_view t) {
    assert(s.size() <= std::size_t(std::numeric_limits<int>::max()) && t.size() <= std::size_t(std::numeric_limits<int>::max()));
    std::size_t prefix = 0, limit = std::min(s.size(), t.size());
    while (prefix < limit && s[prefix] == t[prefix])
        ++prefix;
    s.remove_prefix(prefix);
    t.remove_prefix(prefix);
    std::size_t suffix = 0;
    limit -= prefix;
    while (suffix < limit && s[s.size() - 1 - suffix] == t[t.size() - 1 - suffix])
        ++suffix;
    s.remove_suffix(suffix);
    t.remove_suffix(suffix);
    if (s.size() < t.size()) std::swap(s, t);
    if (t.empty()) return int(s.size());
    if (s.size() != t.size()) {
        std::size_t matched = 0;
        for (char c : s)
            if (c == t[matched] && ++matched == t.size()) return int(s.size() - t.size());
    }
    std::array<bool, 256> present{};
    for (unsigned char c : t)
        present[c] = true;
    bool common = false;
    for (unsigned char c : s)
        if (present[c]) {
            common = true;
            break;
        }
    if (!common) return int(s.size());
    std::size_t work = s.size() * ((t.size() + 63) / 64);
    if (work >= 65536 && s.size() - t.size() <= 1024) {
        int result = edit_distance_detail::wavefront(s, t, 1024, int(std::min<std::size_t>(1 << 22, work / 8)));
        if (result >= 0) return result;
    }
    return edit_distance_detail::bit_parallel(s, t);
}
