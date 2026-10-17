#ifndef NTT_MOD_HPP
#define NTT_MOD_HPP

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>
#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

namespace ntt_mod_detail {
template <std::uint32_t mod>
std::uint32_t power(std::uint32_t a, std::uint32_t k) {
	std::uint32_t result = 1;
	for (; k; k >>= 1, a = std::uint64_t(a) * a % mod)
		if (k & 1) result = std::uint64_t(result) * a % mod;
	return result;
}

template <std::uint32_t mod>
void butterfly(std::uint32_t *a, std::size_t n, std::size_t half,
			   const std::uint32_t *roots, const std::uint32_t *quotients) {
	for (std::size_t start = 0; start < n; start += 2 * half)
		for (std::size_t j = 0; j < half; ++j) {
			std::uint32_t x = a[start + j], v = a[start + j + half];
			std::uint32_t q = std::uint64_t(v) * quotients[j] >> 32;
			std::uint32_t y = v * roots[j] - q * mod;
			if (y >= mod) y -= mod;
			a[start + j] = x + y < mod ? x + y : x + y - mod;
			a[start + j + half] = x >= y ? x - y : x + mod - y;
		}
}

#if defined(__x86_64__) || defined(__i386__)
template <std::uint32_t mod>
__attribute__((target("avx2")))
void butterfly_avx2(std::uint32_t *a, std::size_t n, std::size_t half,
					const std::uint32_t *roots, const std::uint32_t *quotients) {
	const auto p = _mm256_set1_epi32(mod);
	for (std::size_t start = 0; start < n; start += 2 * half)
		for (std::size_t j = 0; j < half; j += 8) {
			auto x = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(a + start + j));
			auto y = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(a + start + j + half));
			auto r = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(roots + j));
			auto q = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(quotients + j));
			auto even = _mm256_mul_epu32(y, q);
			auto odd = _mm256_mul_epu32(_mm256_srli_epi64(y, 32), _mm256_srli_epi64(q, 32));
			q = _mm256_blend_epi32(_mm256_srli_epi64(even, 32), odd, 0xaa);
			y = _mm256_sub_epi32(_mm256_mullo_epi32(y, r), _mm256_mullo_epi32(q, p));
			y = _mm256_min_epu32(y, _mm256_sub_epi32(y, p));
			auto sum = _mm256_add_epi32(x, y), difference = _mm256_sub_epi32(x, y);
			sum = _mm256_min_epu32(sum, _mm256_sub_epi32(sum, p));
			difference = _mm256_min_epu32(difference, _mm256_add_epi32(difference, p));
			_mm256_storeu_si256(reinterpret_cast<__m256i *>(a + start + j), sum);
			_mm256_storeu_si256(reinterpret_cast<__m256i *>(a + start + j + half), difference);
		}
}
#endif

template <std::uint32_t mod, std::uint32_t root>
void transform(std::vector<std::uint32_t> &a, bool inverse, [[maybe_unused]] bool vectorized) {
	static_assert(mod > 2 && mod < (std::uint32_t(1) << 31));
	std::size_t n = a.size();
	assert(std::has_single_bit(n) && (mod - 1) % n == 0);
	for (std::size_t i = 1, j = 0; i < n; ++i) {
		std::size_t bit = n >> 1;
		for (; j & bit; bit >>= 1) j ^= bit;
		j ^= bit;
		if (i < j) std::swap(a[i], a[j]);
	}
	std::vector<std::uint32_t> roots(n > 16 ? n / 2 : 0), quotients(roots.size());
	for (std::size_t width = 2; width <= n; width *= 2) {
		std::uint32_t step = ntt_mod_detail::power<mod>(root, (mod - 1) / width);
		if (inverse) step = ntt_mod_detail::power<mod>(step, mod - 2);
		if (n <= 16) {
			for (std::size_t start = 0; start < n; start += width) {
				std::uint32_t w = 1;
				for (std::size_t j = 0; j < width / 2; ++j) {
					std::uint32_t x = a[start + j], y = std::uint64_t(a[start + j + width / 2]) * w % mod;
					a[start + j] = x + y < mod ? x + y : x + y - mod;
					a[start + j + width / 2] = x >= y ? x - y : x + mod - y;
					w = std::uint64_t(w) * step % mod;
				}
			}
			continue;
		}
		std::uint32_t w = 1;
		for (std::size_t j = 0; j < width / 2; ++j) {
			roots[j] = w;
			quotients[j] = (std::uint64_t(w) << 32) / mod;
			w = std::uint64_t(w) * step % mod;
		}
#if defined(__x86_64__) || defined(__i386__)
		if (vectorized && width >= 16) ntt_mod_detail::butterfly_avx2<mod>(a.data(), n, width / 2, roots.data(), quotients.data());
		else
#endif
		ntt_mod_detail::butterfly<mod>(a.data(), n, width / 2, roots.data(), quotients.data());
	}
	if (inverse) {
		std::uint32_t inv = ntt_mod_detail::power<mod>(n, mod - 2);
		for (auto &x : a) x = std::uint64_t(x) * inv % mod;
	}
}
}

template <std::uint32_t mod, std::uint32_t root>
void ntt_mod(std::vector<std::uint32_t> &a, bool inverse = false) {
	bool vectorized = false;
#if defined(__x86_64__) || defined(__i386__)
	vectorized = __builtin_cpu_supports("avx2");
#endif
	ntt_mod_detail::transform<mod, root>(a, inverse, vectorized);
}

#endif
