#ifndef LIMB_ARITHMETIC_HPP
#define LIMB_ARITHMETIC_HPP

#include <cstddef>
#include <cstdint>
#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize("O3")
#endif

namespace limb_arithmetic {
using U = std::uint64_t;

#if defined(__x86_64__) || defined(__i386__)
template <U base, bool subtract>
__attribute__((target("avx2"))) inline U blocks(U *a, const U *b, std::size_t &offset, std::size_t n, U carry) {
	const auto limit = _mm256_set1_epi64x(base - 1), modulus = _mm256_set1_epi64x(base);
	const auto one = _mm256_set1_epi64x(1), zero = _mm256_setzero_si256();
	const auto indices = _mm256_setr_epi64x(0, 1, 2, 3);
	for (; offset + 16 <= n; offset += 16) {
		__m256i values[4];
		unsigned generate = 0, propagate = 0;
		for (int j = 0; j < 4; ++j) {
			auto x = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(a + offset + 4 * j));
			auto y = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(b + offset + 4 * j));
			__m256i g, p;
			if constexpr (subtract) {
				values[j] = _mm256_sub_epi64(x, y);
				g = _mm256_cmpgt_epi64(y, x);
				p = _mm256_cmpeq_epi64(x, y);
			} else {
				values[j] = _mm256_add_epi64(x, y);
				g = _mm256_cmpgt_epi64(values[j], limit);
				p = _mm256_cmpeq_epi64(values[j], limit);
			}
			generate |= unsigned(_mm256_movemask_pd(_mm256_castsi256_pd(g))) << (4 * j);
			propagate |= unsigned(_mm256_movemask_pd(_mm256_castsi256_pd(p))) << (4 * j);
		}
		unsigned incoming = (propagate + (generate << 1) + unsigned(carry)) ^ propagate;
		for (int j = 0; j < 4; ++j) {
			auto in = _mm256_and_si256(_mm256_srlv_epi64(_mm256_set1_epi64x(incoming >> (4 * j)), indices), one);
			auto out = _mm256_and_si256(_mm256_srlv_epi64(_mm256_set1_epi64x(incoming >> (4 * j + 1)), indices), one);
			auto correction = _mm256_and_si256(_mm256_sub_epi64(zero, out), modulus);
			auto result = subtract ? _mm256_add_epi64(_mm256_sub_epi64(values[j], in), correction)
								   : _mm256_sub_epi64(_mm256_add_epi64(values[j], in), correction);
			_mm256_storeu_si256(reinterpret_cast<__m256i *>(a + offset + 4 * j), result);
		}
		carry = incoming >> 16;
	}
	return carry;
}
#endif

template <U base>
inline U add(U *a, const U *b, std::size_t n, U carry = 0) {
	static_assert(base > 1 && base <= (U(1) << 62));
	std::size_t i = 0;
#if defined(__x86_64__) || defined(__i386__)
	if (n >= 32 && __builtin_cpu_supports("avx2")) carry = blocks<base, false>(a, b, i, n, carry);
#endif
	for (; i < n; ++i) {
		U value = a[i] + b[i] + carry;
		carry = value >= base;
		a[i] = value - carry * base;
	}
	return carry;
}

template <U base>
inline U subtract(U *a, const U *b, std::size_t n, U borrow = 0) {
	static_assert(base > 1 && base <= (U(1) << 62));
	std::size_t i = 0;
#if defined(__x86_64__) || defined(__i386__)
	if (n >= 32 && __builtin_cpu_supports("avx2")) borrow = blocks<base, true>(a, b, i, n, borrow);
#endif
	for (; i < n; ++i) {
		U value = b[i] + borrow;
		borrow = a[i] < value;
		a[i] = a[i] + borrow * base - value;
	}
	return borrow;
}
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

#endif
