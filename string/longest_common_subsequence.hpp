#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

namespace lcs_detail {
using Word = std::uint64_t;

#if defined(__x86_64__)
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize("O3")
#endif
template <bool reverse>
__attribute__((target("avx2"))) inline void avx2(std::string_view s, const std::array<const Word*, 256>& masks, Word* __restrict row,
                                                 int words) {
    const __m256i ones = _mm256_set1_epi64x(-1), one = _mm256_set1_epi64x(1), positions = _mm256_setr_epi64x(0, 1, 2, 3);
    for (int i = 0, n = int(s.size()); i < n; ++i) {
        const Word* __restrict mask = masks[static_cast<unsigned char>(s[reverse ? n - 1 - i : i])];
        if (!mask) continue;
        unsigned carry = 0;
        for (int j = 0; j < words; j += 4) {
            __m256i v = _mm256_load_si256(reinterpret_cast<const __m256i*>(row + j));
            __m256i u = _mm256_and_si256(v, _mm256_load_si256(reinterpret_cast<const __m256i*>(mask + j))), sum = _mm256_add_epi64(v, u);
            unsigned g = unsigned(_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_or_si256(u, _mm256_andnot_si256(sum, v)))));
            unsigned propagation = unsigned(_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpeq_epi64(sum, ones))));
            unsigned incoming = (propagation + (g << 1) + carry) ^ propagation;
            carry = incoming >> 4;
            __m256i increments = _mm256_and_si256(_mm256_srlv_epi64(_mm256_set1_epi64x(incoming), positions), one);
            sum = _mm256_add_epi64(sum, increments);
            _mm256_store_si256(reinterpret_cast<__m256i*>(row + j), _mm256_or_si256(sum, _mm256_xor_si256(v, u)));
        }
    }
}

template <bool reverse>
__attribute__((target("avx512f"))) inline void avx512(std::string_view s, const std::array<const Word*, 256>& masks, Word* __restrict row,
                                                      int words) {
    const __m512i ones = _mm512_set1_epi64(-1), one = _mm512_set1_epi64(1);
    for (int i = 0, n = int(s.size()); i < n; ++i) {
        const Word* __restrict mask = masks[static_cast<unsigned char>(s[reverse ? n - 1 - i : i])];
        if (!mask) continue;
        unsigned carry = 0;
        for (int j = 0; j < words; j += 8) {
            __m512i v = _mm512_load_si512(row + j), u = _mm512_and_si512(v, _mm512_load_si512(mask + j)), sum = _mm512_add_epi64(v, u);
            unsigned g = unsigned(_mm512_cmp_epu64_mask(sum, u, _MM_CMPINT_LT)), propagation = unsigned(_mm512_cmpeq_epi64_mask(sum, ones));
            unsigned incoming = (propagation + (g << 1) + carry) ^ propagation;
            carry = incoming >> 8;
            sum = _mm512_mask_add_epi64(sum, __mmask8(incoming), sum, one);
            _mm512_store_si512(row + j, _mm512_ternarylogic_epi64(sum, v, u, 0xf6));
        }
    }
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif
#endif

inline bool subsequence(std::string_view s, std::string_view t) {
    if (s.empty()) return true;
    std::size_t i = 0;
    for (char c : t)
        if (s[i] == c && ++i == s.size()) return true;
    return false;
}

inline std::pair<std::string_view, std::string_view> trim(std::string_view& s, std::string_view& t) {
    std::size_t prefix = 0, limit = std::min(s.size(), t.size());
    while (prefix < limit && s[prefix] == t[prefix])
        ++prefix;
    std::string_view head = s.substr(0, prefix);
    s.remove_prefix(prefix);
    t.remove_prefix(prefix);
    std::size_t suffix = 0;
    limit -= prefix;
    while (suffix < limit && s[s.size() - 1 - suffix] == t[t.size() - 1 - suffix])
        ++suffix;
    std::string_view tail = s.substr(s.size() - suffix);
    s.remove_suffix(suffix);
    t.remove_suffix(suffix);
    return {head, tail};
}

inline void filter(std::string_view& s, std::string_view& t, std::string& a, std::string& b) {
    std::array<unsigned char, 256> used{};
    for (unsigned char c : s)
        used[c] |= 1;
    for (unsigned char c : t)
        used[c] |= 2;
    bool change_s = false, change_t = false;
    for (int x : used) {
        change_s |= x == 1;
        change_t |= x == 2;
    }
    if (change_s) {
        a.reserve(s.size());
        for (unsigned char c : s)
            if (used[c] == 3) a.push_back(char(c));
        s = a;
    }
    if (change_t) {
        b.reserve(t.size());
        for (unsigned char c : t)
            if (used[c] == 3) b.push_back(char(c));
        t = b;
    }
}

template <bool reverse = false, int backend = 0>
inline std::vector<Word> profile(std::string_view s, std::string_view t) {
    if (t.empty()) return {};
    const int m = int(t.size()), words = int((t.size() + 63) / 64);
    std::array<int, 256> ids;
    ids.fill(-1);
    int alphabet = 0;
    for (unsigned char c : t)
        if (ids[c] == -1) ids[c] = alphabet++;
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
        std::vector<Word> storage(std::size_t(alphabet + 1) * stride + 7);
        Word* data = storage.data() + ((-reinterpret_cast<std::uintptr_t>(storage.data()) & 63) / sizeof(Word));
        Word* row = data + std::size_t(alphabet) * stride;
        std::fill_n(row, stride, ~Word(0));
        std::array<const Word*, 256> masks{};
        for (int c = 0; c < 256; ++c)
            if (ids[c] != -1) masks[c] = data + std::size_t(ids[c]) * stride;
        for (int j = 0; j < m; ++j) {
            unsigned char c = t[reverse ? m - 1 - j : j];
            data[std::size_t(ids[c]) * stride + j / 64] |= Word(1) << (j % 64);
        }
        if (kernel == 3)
            avx512<reverse>(s, masks, row, stride);
        else
            avx2<reverse>(s, masks, row, stride);
        return std::vector<Word>(row, row + words);
    }
#endif
    std::vector<Word> masks(std::size_t(alphabet) * words), row(words, ~Word(0));
    for (int j = 0; j < m; ++j) {
        unsigned char c = t[reverse ? m - 1 - j : j];
        masks[std::size_t(ids[c]) * words + j / 64] |= Word(1) << (j % 64);
    }
    if (words == 1) {
        Word v = ~Word(0);
        for (int i = 0, n = int(s.size()); i < n; ++i) {
            int id = ids[static_cast<unsigned char>(s[reverse ? n - 1 - i : i])];
            if (id == -1) continue;
            Word u = v & masks[id];
            v = (v + u) | (v ^ u);
        }
        row[0] = v;
        return row;
    }
    for (int i = 0, n = int(s.size()); i < n; ++i) {
        int id = ids[static_cast<unsigned char>(s[reverse ? n - 1 - i : i])];
        if (id == -1) continue;
        const Word* mask = masks.data() + std::size_t(id) * words;
        Word carry = 0;
        for (int j = 0; j < words; ++j) {
            Word v = row[j], u = v & mask[j];
            __uint128_t sum = __uint128_t(v) + u + carry;
            row[j] = Word(sum) | (v ^ u);
            carry = Word(sum >> 64);
        }
    }
    return row;
}

inline int count(const std::vector<Word>& row) {
    int result = 0;
    for (Word v : row)
        result += std::popcount(~v);
    return result;
}

template <int backend = 0>
inline void restore(std::string_view s, std::string_view t, std::string& out) {
    auto [head, tail] = trim(s, t);
    out.append(head);
    if (s.size() < t.size()) std::swap(s, t);
    if (t.empty()) {
        out.append(tail);
        return;
    }
    if (subsequence(t, s))
        out.append(t);
    else if (s.size() <= 4096 / t.size()) {
        int n = int(s.size()), m = int(t.size()), stride = m + 1;
        std::vector<std::uint16_t> dp((n + 1) * stride);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                dp[(i + 1) * stride + j + 1] =
                    s[i] == t[j] ? dp[i * stride + j] + 1 : std::max(dp[i * stride + j + 1], dp[(i + 1) * stride + j]);
        std::size_t start = out.size();
        out.resize(start + dp[n * stride + m]);
        std::size_t pos = out.size();
        while (n && m) {
            if (s[n - 1] == t[m - 1])
                out[--pos] = s[--n], --m;
            else if (dp[(n - 1) * stride + m] >= dp[n * stride + m - 1])
                --n;
            else
                --m;
        }
    } else if (t.size() > 1) {
        int middle = int(s.size() / 2), split = 0, best;
        {
            auto left = profile<false, backend>(s.substr(0, middle), t);
            auto right = profile<true, backend>(s.substr(middle), t);
            int a = 0, b = count(right);
            best = b;
            for (int j = 0, m = int(t.size()); j < m; ++j) {
                a += int((~left[j / 64] >> (j % 64)) & 1);
                int k = m - 1 - j;
                b -= int((~right[k / 64] >> (k % 64)) & 1);
                if (a + b > best) best = a + b, split = j + 1;
            }
        }
        if (best) {
            restore<backend>(s.substr(0, middle), t.substr(0, split), out);
            restore<backend>(s.substr(middle), t.substr(split), out);
        }
    }
    out.append(tail);
}
}  // namespace lcs_detail

inline int lcs_length(std::string_view s, std::string_view t) {
    assert(s.size() <= std::size_t(std::numeric_limits<int>::max()) && t.size() <= std::size_t(std::numeric_limits<int>::max()));
    auto [head, tail] = lcs_detail::trim(s, t);
    int result = int(head.size() + tail.size());
    if (s.size() < t.size()) std::swap(s, t);
    if (t.empty() || lcs_detail::subsequence(t, s)) return result + int(t.size());
    std::string a, b;
    lcs_detail::filter(s, t, a, b);
    if (s.size() < t.size()) std::swap(s, t);
    return result + lcs_detail::count(lcs_detail::profile(s, t));
}

inline std::string lcs(std::string_view s, std::string_view t) {
    assert(s.size() <= std::size_t(std::numeric_limits<int>::max()) && t.size() <= std::size_t(std::numeric_limits<int>::max()));
    std::string a, b, result;
    result.reserve(std::min(s.size(), t.size()));
    auto [head, tail] = lcs_detail::trim(s, t);
    result.append(head);
    if (!s.empty() && !t.empty()) {
        lcs_detail::filter(s, t, a, b);
        lcs_detail::restore(s, t, result);
    }
    result.append(tail);
    return result;
}
