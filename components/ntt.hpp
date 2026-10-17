#ifndef NTT_998244353_HPP
#define NTT_998244353_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <memory>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>
#if defined(__ARM_NEON)
#include <arm_neon.h>
#endif

#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#include "../math/mint.hpp"

#if defined(__GNUC__) && !defined(__clang__) && (defined(__x86_64__) || defined(__i386__))
#pragma GCC push_options
#pragma GCC optimize("O3")
#endif

namespace fps_detail {
struct NTT {
	using T = Mint;
	typedef unsigned int U __attribute__((may_alias));
	using V = std::vector<unsigned int>;
	static constexpr U mod = 998244353, twice = 2 * mod, one = 301989884;
	V roots, inverse_roots, twists, inverse_twists;
	V root_quotients, inverse_root_quotients;
	mutable V cube_roots, inverse_cube_roots, cube_quotients, inverse_cube_quotients;
	void prepare_cubes(int n) const {
		int begin = cube_roots.size();
		if (begin >= n) return;
		cube_roots.resize(n);
		inverse_cube_roots.resize(n);
		cube_quotients.resize(n);
		inverse_cube_quotients.resize(n);
		for (int k = begin; k < n; ++k) {
			U x = (unsigned long long)roots[k] * roots[2 * k] % mod;
			U y = (unsigned long long)inverse_roots[k] * inverse_roots[2 * k] % mod;
			cube_roots[k] = x;
			inverse_cube_roots[k] = y;
			cube_quotients[k] = quotient(x);
			inverse_cube_quotients[k] = quotient(y);
		}
	}
	static U reduce(unsigned long long x) {
		return (x + (unsigned long long)(U(x) * 998244351u) * mod) >> 32;
	}
	static U mul(U a, U b) { return reduce((unsigned long long)a * b); }
	static U add(U a, U b) {
		U c = a + b;
		return c < twice ? c : c - twice;
	}
	static U sub(U a, U b) {
		U c = a - b;
		return c < twice ? c : c + twice;
	}
	static U in(U a) { return mul(a, 932051910); }
	static U out(U a) {
		a = reduce(a);
		return a < mod ? a : a - mod;
	}
	static U quotient(U r) { return ((unsigned long long)r << 32) / mod; }
	static U fixed_mul(U a, U r, U q) {
		return a * r - U(((unsigned long long)a * q) >> 32) * mod;
	}
	NTT(int size, bool composition = true)
		: roots(size, one),
		  inverse_roots(size, one),
		  twists(composition ? size : 0, one),
		  inverse_twists(composition ? size : 0, one),
		  root_quotients(size),
		  inverse_root_quotients(size) {
		for (int h = 1; h < size; h <<= 1) {
			U r = in(T(3).pow((mod - 1) / (4 * h)).val());
			U ir = in(T(3).pow(mod - 1 - (mod - 1) / (4 * h)).val());
			U t = in(T(3).pow((mod - 1) / (2 * h)).val());
			U it = in(T(3).pow(mod - 1 - (mod - 1) / (2 * h)).val());
			U a = in(mod - (mod - 1) / h), b = a;
			for (int i = 0; i < h; ++i) {
				roots[h + i] = mul(roots[i], r);
				inverse_roots[h + i] = mul(inverse_roots[i], ir);
				if (composition) twists[h + i] = a, inverse_twists[h + i] = b;
				a = mul(a, t);
				b = mul(b, it);
			}
		}
		for (int i = 0; i < size; ++i) {
			root_quotients[i] = (roots[i] < mod ? roots[i] : roots[i] - mod) * 998244351u;
			inverse_root_quotients[i] = (inverse_roots[i] < mod ? inverse_roots[i] : inverse_roots[i] - mod) * 998244351u;
			roots[i] = out(roots[i]);
			inverse_roots[i] = out(inverse_roots[i]);
			if (composition) twists[i] = out(twists[i]), inverse_twists[i] = out(inverse_twists[i]);
		}
	}

#if defined(__ARM_NEON)
	using I = uint32x4_t;
	static I mul4(I a, I b) {
		uint64x2_t x = vmull_u32(vget_low_u32(a), vget_low_u32(b));
		uint64x2_t y = vmull_u32(vget_high_u32(a), vget_high_u32(b));
		auto reduce2 = [](uint64x2_t z) {
			uint32x2_t q = vmul_u32(vmovn_u64(z), vdup_n_u32(998244351));
			return vshrn_n_u64(vaddq_u64(z, vmull_u32(q, vdup_n_u32(mod))), 32);
		};
		return vcombine_u32(reduce2(x), reduce2(y));
	}
	static I add4(I a, I b) {
		I c = vaddq_u32(a, b);
		return vminq_u32(c, vsubq_u32(c, vdupq_n_u32(twice)));
	}
	static I sub4(I a, I b) {
		I c = vsubq_u32(a, b);
		return vminq_u32(c, vaddq_u32(c, vdupq_n_u32(twice)));
	}
	static I fixed_mul4(I a, I r, I q) {
		uint64x2_t x = vmull_u32(vget_low_u32(a), vget_low_u32(q));
		uint64x2_t y = vmull_u32(vget_high_u32(a), vget_high_u32(q));
		I t = vcombine_u32(vshrn_n_u64(x, 32), vshrn_n_u64(y, 32));
		return vmlsq_u32(vmulq_u32(a, r), t, vdupq_n_u32(mod));
	}
	template <bool inverse>
	static void fixed_butterfly4(I &x, I &y, I r, I q, bool unit = false) {
		if constexpr (!inverse)
			if (!unit) y = fixed_mul4(y, r, q);
		I z = sub4(x, y);
		x = add4(x, y);
		y = inverse && !unit ? fixed_mul4(z, r, q) : z;
	}
#endif
#if defined(__x86_64__) || defined(__i386__)
	static bool vectorized() {
#if defined(__AVX2__)
		return true;
#else
		static const bool available = __builtin_cpu_supports("avx2");
		return available;
#endif
	}
	using I8 = __m256i;
	__attribute__((target("avx2"))) static I8 load8(const U *a) {
		return _mm256_loadu_si256(reinterpret_cast<const I8 *>(a));
	}
	__attribute__((target("avx2"))) static void store8(U *a, I8 x) {
		_mm256_storeu_si256(reinterpret_cast<I8 *>(a), x);
	}
	__attribute__((target("avx2"))) static I8 mul8(I8 a, I8 b) {
		I8 x = _mm256_mul_epu32(a, b);
		I8 y = _mm256_mul_epu32(_mm256_srli_epi64(a, 32), _mm256_srli_epi64(b, 32));
		I8 r = _mm256_set1_epi32(998244351), p = _mm256_set1_epi32(mod);
		x = _mm256_add_epi64(x, _mm256_mul_epu32(_mm256_mul_epu32(x, r), p));
		y = _mm256_add_epi64(y, _mm256_mul_epu32(_mm256_mul_epu32(y, r), p));
		return _mm256_blend_epi32(_mm256_srli_epi64(x, 32), y, 0xaa);
	}
	__attribute__((target("avx2"))) static I8 fixed_mul8(I8 a, I8 r, I8 q) {
		I8 x = _mm256_srli_epi64(_mm256_mul_epu32(a, q), 32);
		I8 y = _mm256_mul_epu32(_mm256_srli_epi64(a, 32), q);
		I8 t = _mm256_blend_epi32(x, y, 0xaa);
		return _mm256_sub_epi32(_mm256_mullo_epi32(a, r),
								_mm256_mullo_epi32(t, _mm256_set1_epi32(mod)));
	}
	__attribute__((target("avx2"))) static I8 fixed_even8(I8 a, I8 r, I8 q) {
		I8 t = _mm256_srli_epi64(_mm256_mul_epu32(a, q), 32);
		return _mm256_sub_epi64(_mm256_mul_epu32(a, r),
								_mm256_mul_epu32(t, _mm256_set1_epi32(mod)));
	}
	__attribute__((target("avx2"))) static I8 add8(I8 a, I8 b) {
		I8 x = _mm256_add_epi32(a, b);
		return _mm256_min_epu32(x, _mm256_sub_epi32(x, _mm256_set1_epi32(twice)));
	}
	__attribute__((target("avx2"))) static I8 sub8(I8 a, I8 b) {
		I8 x = _mm256_sub_epi32(a, b);
		return _mm256_min_epu32(x, _mm256_add_epi32(x, _mm256_set1_epi32(twice)));
	}

	template <int stride>
	__attribute__((target("avx2"))) static I8 rates8(const V &values, int k) {
		if constexpr (stride == 1)
			return _mm256_cvtepu32_epi64(
				_mm_loadu_si128(reinterpret_cast<const __m128i *>(values.data() + 4 * k)));
		if constexpr (stride == 2)
			return _mm256_permute4x64_epi64(
				_mm256_cvtepu32_epi64(
					_mm_loadl_epi64(reinterpret_cast<const __m128i *>(values.data() + 2 * k))),
				0x50);
		if constexpr (stride == 4) return _mm256_set1_epi32(values[k]);
	}
	template <int stride, bool inverse>
	__attribute__((target("avx2"))) I8 small8(I8 z, int k) const {
		if constexpr (stride == 4) {
			if (k == 0) {
				I8 x = _mm256_permute2x128_si256(z, z, 0x00),
				   y = _mm256_permute2x128_si256(z, z, 0x11);
				return _mm256_blend_epi32(add8(x, y), sub8(x, y), 0xf0);
			}
		}
		const V &w = inverse ? inverse_roots : roots;
		const V &quot = inverse ? inverse_root_quotients : root_quotients;
		I8 r = rates8<stride>(w, k), q = rates8<stride>(quot, k);
		if constexpr (inverse) {
			if constexpr (stride == 2) z = _mm256_shuffle_epi32(z, 0xd8);
			if constexpr (stride == 4)
				z = _mm256_permutevar8x32_epi32(z, _mm256_setr_epi32(0, 4, 1, 5, 2, 6, 3, 7));
			I8 y = _mm256_srli_epi64(z, 32), x = add8(z, y);
			I8 d = fixed_even8(sub8(z, y), r, q);
			z = _mm256_blend_epi32(x, _mm256_slli_epi64(d, 32), 0xaa);
		} else {
			if constexpr (stride == 1) z = _mm256_shuffle_epi32(z, 0xb1);
			if constexpr (stride == 2) z = _mm256_shuffle_epi32(z, 0x72);
			if constexpr (stride == 4)
				z = _mm256_permutevar8x32_epi32(z, _mm256_setr_epi32(4, 0, 5, 1, 6, 2, 7, 3));
			I8 y = fixed_even8(z, r, q);
			I8 x = _mm256_blend_epi32(_mm256_srli_epi64(z, 32), z, 0xaa);
			y = _mm256_blend_epi32(
				y, _mm256_sub_epi32(_mm256_set1_epi32(twice), _mm256_slli_epi64(y, 32)), 0xaa);
			z = add8(x, y);
		}
		if constexpr (stride == 2) z = _mm256_shuffle_epi32(z, 0xd8);
		if constexpr (stride == 4)
			z = _mm256_permutevar8x32_epi32(z, _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7));
		return z;
	}
	template <bool inverse>
	__attribute__((target("avx2"))) void fft8(U *a, int length, int limit = 0) const {
		if (!limit) limit = length;
		if constexpr (!inverse) {
			int h = limit / 2;
			for (; h >= 16; h >>= 2) radix4_8<false>(a, h / 2, length);
			if (h == 8) {
				if (limit >= 32)
					for (int i = 0; i < length; i += 8) {
						I8 x = load8(a + i), p = _mm256_set1_epi32(twice);
						store8(a + i, _mm256_min_epu32(x, _mm256_sub_epi32(x, p)));
					}
				stage8<false>(a, h, length);
			}
		}
		for (int j = 0; j < length; j += 8) {
			I8 z = load8(a + j);
			if constexpr (!inverse)
				if ((limit >= 32) && (__builtin_ctz((unsigned int)limit) & 1)) {
					I8 p = _mm256_set1_epi32(twice);
					z = _mm256_min_epu32(z, _mm256_sub_epi32(z, p));
				}
			if constexpr (inverse)
				z = small8<4, true>(small8<2, true>(small8<1, true>(z, j / 8), j / 8), j / 8);
			else
				z = small8<1, false>(small8<2, false>(small8<4, false>(z, j / 8), j / 8),
									 j / 8);
			store8(a + j, z);
		}
		if constexpr (inverse) {
			int h = 8;
			for (; 4 * h <= limit; h <<= 2) radix4_8<true>(a, h, length);
			if (h < limit) stage8<true>(a, h, length);
		}
	}
	template <bool inverse>
	__attribute__((target("avx2"))) void radix4_8(U *a, int stride, int length) const {
		const V &w = inverse ? inverse_roots : roots;
		const V &q = inverse ? inverse_root_quotients : root_quotients;
		prepare_cubes(length / (4 * stride));
		const V &w3 = inverse ? inverse_cube_roots : cube_roots;
		const V &q3s = inverse ? inverse_cube_quotients : cube_quotients;
		I8 p = _mm256_set1_epi32(twice);
		I8 im = _mm256_set1_epi32(w[1]), iq = _mm256_set1_epi32(q[1]);
		for (int j = 0, k = 0; j < length; j += 4 * stride, ++k) {
			I8 r = _mm256_set1_epi32(w[k]), rq = _mm256_set1_epi32(q[k]);
			I8 r0 = _mm256_set1_epi32(w[2 * k]), q0 = _mm256_set1_epi32(q[2 * k]);
			I8 r3 = _mm256_set1_epi32(w3[k]), q3 = _mm256_set1_epi32(q3s[k]);
			for (int i = 0; i < stride; i += 8) {
				I8 x0 = load8(a + j + i), x1 = load8(a + j + stride + i);
				I8 x2 = load8(a + j + 2 * stride + i), x3 = load8(a + j + 3 * stride + i);
				if constexpr (inverse) {
					I8 u = add8(x0, x1), v = sub8(x0, x1);
					I8 s = add8(x2, x3);
					I8 t = fixed_mul8(_mm256_sub_epi32(_mm256_add_epi32(x2, p), x3), im, iq);
					x0 = add8(u, s);
					x1 = _mm256_add_epi32(v, t);
					x2 = _mm256_sub_epi32(_mm256_add_epi32(u, p), s);
					x3 = _mm256_sub_epi32(_mm256_add_epi32(v, p), t);
					if (k) {
						x1 = fixed_mul8(x1, r0, q0);
						x2 = fixed_mul8(x2, r, rq);
						x3 = fixed_mul8(x3, r3, q3);
					} else {
						x1 = _mm256_min_epu32(x1, _mm256_sub_epi32(x1, p));
						x2 = _mm256_min_epu32(x2, _mm256_sub_epi32(x2, p));
						x3 = _mm256_min_epu32(x3, _mm256_sub_epi32(x3, p));
					}
				} else {
					x0 = _mm256_min_epu32(x0, _mm256_sub_epi32(x0, p));
					if (!k) {
						x1 = _mm256_min_epu32(x1, _mm256_sub_epi32(x1, p));
						x2 = _mm256_min_epu32(x2, _mm256_sub_epi32(x2, p));
						x3 = _mm256_min_epu32(x3, _mm256_sub_epi32(x3, p));
					}
					if (k) {
						x1 = fixed_mul8(x1, r0, q0);
						x2 = fixed_mul8(x2, r, rq);
						x3 = fixed_mul8(x3, r3, q3);
					}
					I8 u = add8(x0, x2), v = sub8(x0, x2), s = add8(x1, x3);
					I8 t = fixed_mul8(_mm256_sub_epi32(_mm256_add_epi32(x1, p), x3), im, iq);
					x0 = _mm256_add_epi32(u, s);
					x1 = _mm256_sub_epi32(_mm256_add_epi32(u, p), s);
					x2 = _mm256_add_epi32(v, t);
					x3 = _mm256_sub_epi32(_mm256_add_epi32(v, p), t);
				}
				store8(a + j + i, x0);
				store8(a + j + stride + i, x1);
				store8(a + j + 2 * stride + i, x2);
				store8(a + j + 3 * stride + i, x3);
			}
		}
	}
	template <bool inverse>
	__attribute__((target("avx2"))) void stage8(U *a, int stride, int length, int first = 0) const {
		const V &w = inverse ? inverse_roots : roots;
		const V &quot = inverse ? inverse_root_quotients : root_quotients;
		if (stride < 8) {
			for (int j = 0; j < length; j += 8) {
				I8 z = load8(a + j);
				if (stride == 1)
					z = small8<1, inverse>(z, j / 8);
				else if (stride == 2)
					z = small8<2, inverse>(z, j / 8);
				else
					z = small8<4, inverse>(z, j / 8);
				store8(a + j, z);
			}
		} else {
			for (int j = 0, k = first; j < length; j += 2 * stride, ++k) {
				I8 r = _mm256_set1_epi32(w[k]), q = _mm256_set1_epi32(quot[k]);
				for (int i = 0; i < stride; i += 8) {
					I8 x = load8(a + j + i), y = load8(a + j + stride + i);
					if constexpr (!inverse)
						if (k) y = fixed_mul8(y, r, q);
					I8 z = sub8(x, y);
					x = add8(x, y);
					if constexpr (inverse)
						if (k) z = fixed_mul8(z, r, q);
					store8(a + j + i, x);
					store8(a + j + stride + i, z);
				}
			}
		}
	}

	template <bool inverse, bool lazy = false>
	__attribute__((target("avx2"), always_inline)) static void fixed_butterfly8(I8 &x, I8 &y, I8 r, I8 q, bool unit = false) {
		if constexpr (lazy && !inverse) {
			I8 p = _mm256_set1_epi32(twice);
			x = _mm256_min_epu32(x, _mm256_sub_epi32(x, p));
			y = unit ? _mm256_min_epu32(y, _mm256_sub_epi32(y, p)) : fixed_mul8(y, r, q);
			I8 z = _mm256_sub_epi32(_mm256_add_epi32(x, p), y);
			x = _mm256_add_epi32(x, y);
			y = z;
		} else {
			if constexpr (!inverse)
				if (!unit) y = fixed_mul8(y, r, q);
			I8 z = sub8(x, y);
			x = add8(x, y);
			y = inverse && !unit ? fixed_mul8(z, r, q) : z;
		}
	}
	template <bool inverse, bool lazy = false>
	__attribute__((target("avx2"))) void stage_pair8(U *a, int stride, int length, int first = 0) const {
		const V &w = inverse ? inverse_roots : roots;
		const V &q = inverse ? inverse_root_quotients : root_quotients;
		for (int j = 0, k = first; j < length; j += 4 * stride, ++k) {
			I8 r = _mm256_set1_epi32(w[k]), rq = _mm256_set1_epi32(q[k]);
			I8 r0 = _mm256_set1_epi32(w[2 * k]), q0 = _mm256_set1_epi32(q[2 * k]);
			I8 r1 = _mm256_set1_epi32(w[2 * k + 1]), q1 = _mm256_set1_epi32(q[2 * k + 1]);
			for (int i = 0; i < stride; i += 8) {
				I8 x0 = load8(a + j + i), x1 = load8(a + j + stride + i);
				I8 x2 = load8(a + j + 2 * stride + i), x3 = load8(a + j + 3 * stride + i);
				if constexpr (inverse) {
					fixed_butterfly8<inverse, lazy>(x0, x1, r0, q0, k == 0);
					fixed_butterfly8<inverse, lazy>(x2, x3, r1, q1, false);
					fixed_butterfly8<inverse, lazy>(x0, x2, r, rq, k == 0);
					fixed_butterfly8<inverse, lazy>(x1, x3, r, rq, k == 0);
				} else {
					fixed_butterfly8<inverse, lazy>(x0, x2, r, rq, k == 0);
					fixed_butterfly8<inverse, lazy>(x1, x3, r, rq, k == 0);
					fixed_butterfly8<inverse, lazy>(x0, x1, r0, q0, k == 0);
					fixed_butterfly8<inverse, lazy>(x2, x3, r1, q1, false);
				}
				store8(a + j + i, x0);
				store8(a + j + stride + i, x1);
				store8(a + j + 2 * stride + i, x2);
				store8(a + j + 3 * stride + i, x3);
			}
		}
	}
	__attribute__((target("avx2"))) static void fold8(U *a, int columns) {
		for (int i = 0; i < columns; i += 8) {
			I8 x = add8(load8(a + i), load8(a + columns + i));
			store8(a + i, x);
			store8(a + columns + i, x);
		}
	}
	__attribute__((target("avx2"))) void project8_8(U *a, int rows) const {
		for (int j = 0; j < rows; ++j) {
			I8 x0 = load8(a + 16 * j + 0);
			x0 = small8<4, true>(small8<2, true>(small8<1, true>(x0, 0), 0), 0);
			I8 x1 = load8(a + 16 * j + 8);
			x1 = small8<4, true>(small8<2, true>(small8<1, true>(x1, 1), 1), 1);
			x0 = add8(x0, x1);
			x1 = x0;
			x0 = small8<1, false>(small8<2, false>(small8<4, false>(x0, 0), 0), 0);
			store8(a + 16 * j + 0, x0);
			x1 = small8<1, false>(small8<2, false>(small8<4, false>(x1, 1), 1), 1);
			store8(a + 16 * j + 8, x1);
		}
	}
	__attribute__((target("avx2"))) void project16_8(U *a, int rows) const {
		for (int j = 0; j < rows; ++j) {
			I8 x0 = load8(a + 32 * j + 0);
			x0 = small8<4, true>(small8<2, true>(small8<1, true>(x0, 0), 0), 0);
			I8 x1 = load8(a + 32 * j + 8);
			x1 = small8<4, true>(small8<2, true>(small8<1, true>(x1, 1), 1), 1);
			I8 x2 = load8(a + 32 * j + 16);
			x2 = small8<4, true>(small8<2, true>(small8<1, true>(x2, 2), 2), 2);
			I8 x3 = load8(a + 32 * j + 24);
			x3 = small8<4, true>(small8<2, true>(small8<1, true>(x3, 3), 3), 3);
			fixed_butterfly8<true>(x0, x1, _mm256_set1_epi32(inverse_roots[0]), _mm256_set1_epi32(inverse_root_quotients[0]), true);
			fixed_butterfly8<true>(x2, x3, _mm256_set1_epi32(inverse_roots[1]), _mm256_set1_epi32(inverse_root_quotients[1]), false);
			x0 = add8(x0, x2);
			x2 = x0;
			x1 = add8(x1, x3);
			x3 = x1;
			fixed_butterfly8<false>(x0, x1, _mm256_set1_epi32(roots[0]), _mm256_set1_epi32(root_quotients[0]), true);
			fixed_butterfly8<false>(x2, x3, _mm256_set1_epi32(roots[1]), _mm256_set1_epi32(root_quotients[1]), false);
			x0 = small8<1, false>(small8<2, false>(small8<4, false>(x0, 0), 0), 0);
			store8(a + 32 * j + 0, x0);
			x1 = small8<1, false>(small8<2, false>(small8<4, false>(x1, 1), 1), 1);
			store8(a + 32 * j + 8, x1);
			x2 = small8<1, false>(small8<2, false>(small8<4, false>(x2, 2), 2), 2);
			store8(a + 32 * j + 16, x2);
			x3 = small8<1, false>(small8<2, false>(small8<4, false>(x3, 3), 3), 3);
			store8(a + 32 * j + 24, x3);
		}
	}
	__attribute__((target("avx2"))) static void scale_rows8(U *a, int rows, int columns, const U *factors) {
		if (columns == 2) {
			int j = 0;
			for (; j + 1 < rows; j += 2) {
				U w0 = factors[j], w1 = factors[j + 1], q0 = quotient(w0), q1 = quotient(w1);
				I8 w = _mm256_setr_epi32(w0, w0, w0, w0, w1, w1, w1, w1);
				I8 q = _mm256_setr_epi32(q0, q0, q0, q0, q1, q1, q1, q1);
				store8(a + 4 * j, fixed_mul8(load8(a + 4 * j), w, q));
			}
			if (j < rows)
				for (int i = 0; i < 4; ++i) a[4 * j + i] = fixed_mul(a[4 * j + i], factors[j], quotient(factors[j]));
		} else {
			for (int j = 0; j < rows; ++j) {
				U w = factors[j];
				I8 v = _mm256_set1_epi32(w), q = _mm256_set1_epi32(quotient(w));
				for (int i = 0; i < 2 * columns; i += 8) store8(a + 2 * columns * j + i, fixed_mul8(load8(a + 2 * columns * j + i), v, q));
			}
		}
	}
	__attribute__((target("avx2"))) void project4_8(U *a, int rows) const {
		for (int j = 0; j < rows; ++j) {
			I8 x = small8<2, true>(small8<1, true>(load8(a + 8 * j), 0), 0);
			x = add8(x, _mm256_permute2x128_si256(x, x, 1));
			x = small8<1, false>(small8<2, false>(x, 0), 0);
			store8(a + 8 * j, x);
		}
	}
	__attribute__((target("avx2"))) static void load_pairs(const U *a, I8 &x, I8 &y) {
		I8 index = _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7);
		I8 u = _mm256_permutevar8x32_epi32(load8(a), index),
		   v = _mm256_permutevar8x32_epi32(load8(a + 8), index);
		x = _mm256_permute2x128_si256(u, v, 0x20);
		y = _mm256_permute2x128_si256(u, v, 0x31);
	}
	__attribute__((target("avx2"))) static void store_pairs(U *a, I8 x, I8 y) {
		I8 index = _mm256_setr_epi32(0, 4, 1, 5, 2, 6, 3, 7);
		store8(a, _mm256_permutevar8x32_epi32(_mm256_permute2x128_si256(x, y, 0x20), index));
		store8(a + 8,
			   _mm256_permutevar8x32_epi32(_mm256_permute2x128_si256(x, y, 0x31), index));
	}
	__attribute__((target("avx2"))) static void scale_copy8(U *dst, const U *src, int length, U r) {
		U quotient_r = quotient(r);
		I8 v = _mm256_set1_epi32(r), q = _mm256_set1_epi32(quotient_r);
		int i = 0;
		for (; i + 8 <= length; i += 8) store8(dst + i, fixed_mul8(load8(src + i), v, q));
		for (; i < length; ++i) dst[i] = fixed_mul(src[i], r, quotient_r);
	}
	__attribute__((target("avx2"))) static void scale_normalize8(U *a, int n, U r) {
		U qr = quotient(r);
		I8 v = _mm256_set1_epi32(r), q = _mm256_set1_epi32(qr), p = _mm256_set1_epi32(mod);
		int i = 0;
		for (; i + 8 <= n; i += 8) {
			I8 x = fixed_mul8(load8(a + i), v, q);
			store8(a + i, _mm256_min_epu32(x, _mm256_sub_epi32(x, p)));
		}
		for (; i < n; ++i) {
			U x = fixed_mul(a[i], r, qr);
			a[i] = x < mod ? x : x - mod;
		}
	}
	__attribute__((target("avx2"))) static void square_fold8(U *a, const U *q, int n, bool negative) {
		for (int i = 0; i < n; i += 8) {
			I8 x = load8(a + i), y = load8(a + n + i);
			I8 u = mul8(load8(q + i), load8(q + i));
			I8 v = mul8(load8(q + n + i), load8(q + n + i));
			x = add8(x, x);
			store8(a + i, negative ? sub8(x, u) : add8(x, u));
			store8(a + n + i, add8(add8(y, y), v));
		}
	}
	__attribute__((target("avx2"))) static void scale8(U *a, int length, U r) {
		I8 v = _mm256_set1_epi32(r), q = _mm256_set1_epi32(quotient(r));
		int i = 0;
		for (; i + 8 <= length; i += 8) store8(a + i, fixed_mul8(load8(a + i), v, q));
		for (; i < length; ++i) a[i] = fixed_mul(a[i], r, quotient(r));
	}
	__attribute__((target("avx2"))) static void pair_product8(U *a, const U *b, int length) {
		I8 ni = _mm256_set1_epi32(998244351), modulus = _mm256_set1_epi32(mod), p = _mm256_set1_epi32(twice);
		I8 index = _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7);
		for (int i = 0; i < length; i += 8) {
			I8 x = load8(b + 2 * i), y = load8(b + 2 * i + 8);
			x = _mm256_min_epu32(x, _mm256_sub_epi32(x, p));
			y = _mm256_min_epu32(y, _mm256_sub_epi32(y, p));
			x = _mm256_mul_epu32(x, _mm256_srli_epi64(x, 32));
			y = _mm256_mul_epu32(y, _mm256_srli_epi64(y, 32));
			x = _mm256_add_epi64(x, _mm256_mul_epu32(_mm256_mul_epu32(x, ni), modulus));
			y = _mm256_add_epi64(y, _mm256_mul_epu32(_mm256_mul_epu32(y, ni), modulus));
			x = _mm256_blend_epi32(_mm256_srli_epi64(x, 32), y, 0xaa);
			store8(a + i, _mm256_permutevar8x32_epi32(x, index));
		}
	}
	__attribute__((target("avx2"))) static void pair_multiply8(U *a, const U *b, int length) {
		I8 index0 = _mm256_setr_epi32(0, 0, 1, 1, 2, 2, 3, 3), index1 = _mm256_setr_epi32(4, 4, 5, 5, 6, 6, 7, 7);
		I8 p = _mm256_set1_epi32(twice);
		for (int i = 0; i < length; i += 8) {
			I8 r = load8(b + i), x = load8(a + 2 * i), y = load8(a + 2 * i + 8);
			x = _mm256_min_epu32(x, _mm256_sub_epi32(x, p));
			y = _mm256_min_epu32(y, _mm256_sub_epi32(y, p));
			x = _mm256_shuffle_epi32(x, 0xb1);
			y = _mm256_shuffle_epi32(y, 0xb1);
			store8(a + 2 * i, mul8(x, _mm256_permutevar8x32_epi32(r, index0)));
			store8(a + 2 * i + 8, mul8(y, _mm256_permutevar8x32_epi32(r, index1)));
		}
	}
	__attribute__((target("avx2"))) static void difference8(U *a, const U *b, const U *c,
															int length) {
		for (int i = 0; i < length; i += 8) store8(a + i, sub8(load8(b + i), load8(c + i)));
	}
	__attribute__((target("avx2"))) static void project1_8(U *a, int rows) {
		for (int j = 0; j < rows; j += 8) {
			I8 x, y;
			load_pairs(a + 2 * j, x, y);
			x = add8(x, y);
			store_pairs(a + 2 * j, x, x);
		}
	}
	__attribute__((target("avx2"))) void project2_8(U *a, int rows) const {
		const I8 index = _mm256_setr_epi32(0, 4, 1, 5, 2, 6, 3, 7);
		for (int j = 0; j < rows; j += 8) {
			U *p = a + 4 * j;
			I8 u = load8(p), v = load8(p + 8), w = load8(p + 16), z = load8(p + 24);
			I8 uv0 = _mm256_unpacklo_epi32(u, v), uv1 = _mm256_unpackhi_epi32(u, v);
			I8 wz0 = _mm256_unpacklo_epi32(w, z), wz1 = _mm256_unpackhi_epi32(w, z);
			I8 a0 = _mm256_permutevar8x32_epi32(_mm256_unpacklo_epi64(uv0, wz0), index);
			I8 a1 = _mm256_permutevar8x32_epi32(_mm256_unpackhi_epi64(uv0, wz0), index);
			I8 a2 = _mm256_permutevar8x32_epi32(_mm256_unpacklo_epi64(uv1, wz1), index);
			I8 a3 = _mm256_permutevar8x32_epi32(_mm256_unpackhi_epi64(uv1, wz1), index);
			I8 x = add8(add8(a0, a1), add8(a2, a3));
			I8 y =
				add8(sub8(a0, a1), fixed_mul8(sub8(a2, a3), _mm256_set1_epi32(inverse_roots[1]),
											  _mm256_set1_epi32(inverse_root_quotients[1])));
			a0 = add8(x, y);
			a1 = sub8(x, y);
			y = fixed_mul8(y, _mm256_set1_epi32(roots[1]),
						   _mm256_set1_epi32(root_quotients[1]));
			a2 = add8(x, y);
			a3 = sub8(x, y);
			uv0 = _mm256_unpacklo_epi32(a0, a1);
			uv1 = _mm256_unpackhi_epi32(a0, a1);
			wz0 = _mm256_unpacklo_epi32(a2, a3);
			wz1 = _mm256_unpackhi_epi32(a2, a3);
			u = _mm256_unpacklo_epi64(uv0, wz0);
			v = _mm256_unpackhi_epi64(uv0, wz0);
			w = _mm256_unpacklo_epi64(uv1, wz1);
			z = _mm256_unpackhi_epi64(uv1, wz1);
			store8(p, _mm256_permute2x128_si256(u, v, 0x20));
			store8(p + 8, _mm256_permute2x128_si256(w, z, 0x20));
			store8(p + 16, _mm256_permute2x128_si256(u, v, 0x31));
			store8(p + 24, _mm256_permute2x128_si256(w, z, 0x31));
		}
	}
#endif

#if defined(__ARM_NEON)
	template <int stride, bool inverse>
	void small4(I &a, I &b, int k) const {
		const V &w = inverse ? inverse_roots : roots;
		const V &quot = inverse ? inverse_root_quotients : root_quotients;
		if constexpr (stride == 4) {
			if (k == 0) {
				I x = sub4(a, b);
				a = add4(a, b);
				b = x;
			} else
				fixed_butterfly4<inverse>(a, b, vdupq_n_u32(w[k]), vdupq_n_u32(quot[k]));
		}
		if constexpr (stride == 2) {
			I x = vcombine_u32(vget_low_u32(a), vget_low_u32(b));
			I y = vcombine_u32(vget_high_u32(a), vget_high_u32(b));
			fixed_butterfly4<inverse>(
				x, y, vcombine_u32(vdup_n_u32(w[2 * k]), vdup_n_u32(w[2 * k + 1])),
				vcombine_u32(vdup_n_u32(quot[2 * k]), vdup_n_u32(quot[2 * k + 1])));
			a = vcombine_u32(vget_low_u32(x), vget_low_u32(y));
			b = vcombine_u32(vget_high_u32(x), vget_high_u32(y));
		}
		if constexpr (stride == 1) {
			uint32x4x2_t z = vuzpq_u32(a, b);
			fixed_butterfly4<inverse>(z.val[0], z.val[1], vld1q_u32(w.data() + 4 * k),
									  vld1q_u32(quot.data() + 4 * k));
			z = vzipq_u32(z.val[0], z.val[1]);
			a = z.val[0];
			b = z.val[1];
		}
	}
	template <bool inverse>
	void fft4(U *a, int length, int limit = 0) const {
		if (!limit) limit = length;
		if constexpr (!inverse) {
			int h = limit / 2;
			for (; h >= 16; h >>= 2) stage_pair<false>(a, h / 2, length);
			if (h == 8) stage<false>(a, h, length);
		}
		for (int j = 0; j < length; j += 8) {
			I x = vld1q_u32(a + j), y = vld1q_u32(a + j + 4);
			if constexpr (inverse) {
				small4<1, true>(x, y, j / 8);
				small4<2, true>(x, y, j / 8);
				small4<4, true>(x, y, j / 8);
			} else {
				small4<4, false>(x, y, j / 8);
				small4<2, false>(x, y, j / 8);
				small4<1, false>(x, y, j / 8);
			}
			vst1q_u32(a + j, x);
			vst1q_u32(a + j + 4, y);
		}
		if constexpr (inverse) {
			int h = 8;
			for (; 4 * h <= limit; h <<= 2) stage_pair<true>(a, h, length);
			if (h < limit) stage<true>(a, h, length);
		}
	}
#endif
	template <bool inverse>
	void stage(U *a, int stride, int length, int first = 0) const {
		const V &w = inverse ? inverse_roots : roots;
		const V &quot = inverse ? inverse_root_quotients : root_quotients;
#if defined(__x86_64__) || defined(__i386__)
		if (length >= 8 && vectorized()) {
			stage8<inverse>(a, stride, length, first);
			return;
		}
#endif
		for (int j = 0, k = first; j < length; j += 2 * stride, ++k) {
			U r = w[k], q = quot[k];
			int i = 0;
#if defined(__ARM_NEON)
			for (; i + 4 <= stride; i += 4) {
				I x = vld1q_u32(a + j + i), y = vld1q_u32(a + j + stride + i);
				if constexpr (!inverse)
					if (k) y = fixed_mul4(y, vdupq_n_u32(r), vdupq_n_u32(q));
				I z = sub4(x, y);
				x = add4(x, y);
				if constexpr (inverse)
					if (k) z = fixed_mul4(z, vdupq_n_u32(r), vdupq_n_u32(q));
				y = z;
				vst1q_u32(a + j + i, x);
				vst1q_u32(a + j + stride + i, y);
			}
#endif
			for (; i < stride; ++i) {
				U x = a[j + i], y = a[j + stride + i];
				if constexpr (!inverse)
					if (k) y = fixed_mul(y, r, q);
				a[j + i] = add(x, y);
				a[j + stride + i] = inverse && k ? fixed_mul(sub(x, y), r, q) : sub(x, y);
			}
		}
	}

	template <bool inverse>
	void stage_pair(U *a, int stride, int length, int first = 0) const {
#if defined(__x86_64__) || defined(__i386__)
		if (stride >= 8 && vectorized()) {
			stage_pair8<inverse>(a, stride, length, first);
			return;
		}
#endif
#if defined(__ARM_NEON)
		if (stride >= 4) {
			const V &w = inverse ? inverse_roots : roots;
			const V &q = inverse ? inverse_root_quotients : root_quotients;
			for (int j = 0, k = first; j < length; j += 4 * stride, ++k) {
				I r = vdupq_n_u32(w[k]), rq = vdupq_n_u32(q[k]);
				I r0 = vdupq_n_u32(w[2 * k]), q0 = vdupq_n_u32(q[2 * k]);
				I r1 = vdupq_n_u32(w[2 * k + 1]), q1 = vdupq_n_u32(q[2 * k + 1]);
				for (int i = 0; i < stride; i += 4) {
					I x0 = vld1q_u32(a + j + i), x1 = vld1q_u32(a + j + stride + i);
					I x2 = vld1q_u32(a + j + 2 * stride + i),
					  x3 = vld1q_u32(a + j + 3 * stride + i);
					if constexpr (inverse) {
						fixed_butterfly4<true>(x0, x1, r0, q0, k == 0);
						fixed_butterfly4<true>(x2, x3, r1, q1);
						fixed_butterfly4<true>(x0, x2, r, rq, k == 0);
						fixed_butterfly4<true>(x1, x3, r, rq, k == 0);
					} else {
						fixed_butterfly4<false>(x0, x2, r, rq, k == 0);
						fixed_butterfly4<false>(x1, x3, r, rq, k == 0);
						fixed_butterfly4<false>(x0, x1, r0, q0, k == 0);
						fixed_butterfly4<false>(x2, x3, r1, q1);
					}
					vst1q_u32(a + j + i, x0);
					vst1q_u32(a + j + stride + i, x1);
					vst1q_u32(a + j + 2 * stride + i, x2);
					vst1q_u32(a + j + 3 * stride + i, x3);
				}
			}
			return;
		}
#endif
		if constexpr (inverse)
			stage<true>(a, stride, length, 2 * first), stage<true>(a, 2 * stride, length, first);
		else
			stage<false>(a, 2 * stride, length, first), stage<false>(a, stride, length, 2 * first);
	}
	template <bool inverse>
	void fft(U *a, int length, int limit = 0) const {
		if (!limit) limit = length;
#if defined(__ARM_NEON)
		if (limit >= 8) {
			fft4<inverse>(a, length, limit);
			return;
		}
#elif defined(__x86_64__) || defined(__i386__)
		if (limit >= 8 && vectorized()) {
			fft8<inverse>(a, length, limit);
			return;
		}
#endif
		if constexpr (inverse) {
			for (int h = 1; h < limit; h <<= 1) stage<true>(a, h, length);
		} else {
			for (int h = limit / 2; h; h >>= 1) stage<false>(a, h, length);
		}
	}
#if defined(__x86_64__) || defined(__i386__)
	__attribute__((target("avx2"))) static void product8(U *dst, const U *a, const U *b, int n) {
		int i = 0;
		for (; i + 8 <= n; i += 8) {
			I8 x = mul8(load8(a + i), load8(b + i));
			store8(dst + i, x);
		}
		for (; i < n; ++i) dst[i] = mul(a[i], b[i]);
	}
	__attribute__((target("avx2"))) static void weighted8(U *dst, const U *src, int n, int first) {
		alignas(32) U factors[8];
		for (int j = 0; j < 8; ++j) factors[j] = in(first + j);
		I8 r = load8(factors), step = _mm256_set1_epi32(in(8));
		int i = 0;
		for (; i + 8 <= n; i += 8) {
			store8(dst + i, mul8(load8(src + i), r));
			r = add8(r, step);
		}
		for (; i < n; ++i) dst[i] = mul(src[i], in(first + i));
	}
	__attribute__((target("avx2"), always_inline)) static I8 inverse8(I8 a) {
		I8 a2 = mul8(a, a), a3 = mul8(a2, a);
		I8 a7 = mul8(mul8(a3, a3), a), a31 = a7;
		for (int repeat_index = 0; repeat_index < 2; ++repeat_index) a31 = mul8(a31, a31);
		a31 = mul8(a31, a3);
		I8 r = a7;
		for (int repeat_index = 0; repeat_index < 3; ++repeat_index) r = mul8(r, r);
		r = mul8(r, a3);
		r = mul8(r, r);
		for (int j = 0; j < 4; ++j) {
			for (int repeat_index = 0; repeat_index < 5; ++repeat_index) r = mul8(r, r);
			r = mul8(r, a31);
		}
		for (int repeat_index = 0; repeat_index < 3; ++repeat_index) r = mul8(r, r);
		return mul8(r, a7);
	}
	__attribute__((target("avx2"))) static void batch_inverse8(U *dst, const U *src, int n) {
		if (n < 16) {
			for (int i = 0; i < n; ++i) dst[i] = T::word(src[i]).inv().v;
			return;
		}
		int count = n / 16 * 16;
		I8 p = load8(src), q = load8(src + 8);
		store8(dst, p);
		store8(dst + 8, q);
		for (int i = 16; i < count; i += 16) {
			p = mul8(p, load8(src + i));
			q = mul8(q, load8(src + i + 8));
			store8(dst + i, p);
			store8(dst + i + 8, q);
		}
		alignas(32) U lanes[16];
		store8(lanes, p);
		store8(lanes + 8, q);
		for (auto x : lanes) assert(out(x) != 0);
		I8 inverse = inverse8(p), other = inverse8(q);
		for (int i = count - 16; i >= 16; i -= 16) {
			store8(dst + i, mul8(inverse, load8(dst + i - 16)));
			store8(dst + i + 8, mul8(other, load8(dst + i - 8)));
			inverse = mul8(inverse, load8(src + i));
			other = mul8(other, load8(src + i + 8));
		}
		store8(dst, inverse);
		store8(dst + 8, other);
		for (int i = count; i < n; ++i) dst[i] = T::word(src[i]).inv().v;
	}
	__attribute__((target("avx2"))) static void evaluate_leaf8(U *dst, const U *cur, const U *tree, int degree, const U *points, int count) {
		I8 p = _mm256_set1_epi32(mod), pp = _mm256_set1_epi32(twice);
		I8 x = load8(cur);
		x = _mm256_min_epu32(x, _mm256_sub_epi32(x, pp));
		x = _mm256_min_epu32(x, _mm256_sub_epi32(x, p));
		alignas(32) U window[16] = {}, coefficients[8] = {}, polynomial[8], tail[8] = {};
		store8(window + 8, x);
		std::copy_n(tree, std::min(degree, 8), coefficients);
		I8 y = load8(coefficients);
		y = _mm256_min_epu32(y, _mm256_sub_epi32(y, pp));
		store8(coefficients, _mm256_min_epu32(y, _mm256_sub_epi32(y, p)));
		I8 e = _mm256_setzero_si256(), o = e;
		for (int i = 0; i < 8; ++i) {
			I8 z = load8(window + 8 - i), r = _mm256_set1_epi32(coefficients[i]);
			e = _mm256_add_epi64(e, _mm256_mul_epu32(z, r));
			o = _mm256_add_epi64(o, _mm256_mul_epu32(_mm256_srli_epi64(z, 32), r));
		}
		I8 ni = _mm256_set1_epi32(998244351);
		e = _mm256_add_epi64(e, _mm256_mul_epu32(_mm256_mul_epu32(e, ni), p));
		o = _mm256_add_epi64(o, _mm256_mul_epu32(_mm256_mul_epu32(o, ni), p));
		x = _mm256_blend_epi32(_mm256_srli_epi64(e, 32), o, 0xaa);
		x = _mm256_min_epu32(x, _mm256_sub_epi32(x, pp));
		store8(polynomial, _mm256_min_epu32(x, _mm256_sub_epi32(x, p)));
		if (count == 8) x = load8(points);
		else {
			std::copy_n(points, count, tail);
			x = load8(tail);
		}
		I8 value = _mm256_setzero_si256();
		for (int i = 0; i < 8; ++i) value = add8(mul8(value, x), _mm256_set1_epi32(polynomial[i]));
		if (count == 8) store8(dst, value);
		else {
			store8(tail, value);
			std::copy_n(tail, count, dst);
		}
	}
	__attribute__((target("avx2"), always_inline)) static U sum_product8(I8 a, I8 b) {
		I8 p = _mm256_set1_epi32(mod);
		a = _mm256_min_epu32(a, _mm256_sub_epi32(a, p));
		b = _mm256_min_epu32(b, _mm256_sub_epi32(b, p));
		I8 sum = _mm256_add_epi64(_mm256_mul_epu32(a, b), _mm256_mul_epu32(_mm256_srli_epi64(a, 32), _mm256_srli_epi64(b, 32)));
		sum = _mm256_add_epi64(sum, _mm256_permute2x128_si256(sum, sum, 1));
		sum = _mm256_add_epi64(sum, _mm256_shuffle_epi32(sum, 0x4e));
		U result = reduce((unsigned long long)_mm_cvtsi128_si64(_mm256_castsi256_si128(sum)));
		return result < twice ? result : result - twice;
	}
	__attribute__((target("avx2"))) static void interpolate_leaf8(U *dst, const U *q, const U *points, const U *weights) {
		I8 x = load8(points), w = load8(weights), c = _mm256_set1_epi32(one);
		I8 sum = add8(w, _mm256_permute2x128_si256(w, w, 1));
		sum = add8(sum, _mm256_shuffle_epi32(sum, 0x4e));
		sum = add8(sum, _mm256_shuffle_epi32(sum, 0xb1));
		dst[7] = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum));
		for (int k = 7 - 1; k >= 0; --k) {
			c = add8(mul8(c, x), _mm256_set1_epi32(q[k + 1]));
			dst[k] = sum_product8(c, w);
		}
	}
	__attribute__((target("avx2"))) static void evaluate_middle8(U *left, U *right, const U *cur, const U *ql, int dl, const U *qr, int dr) {
		I8 p = _mm256_set1_epi32(mod), pp = _mm256_set1_epi32(twice), ni = _mm256_set1_epi32(998244351);
		alignas(32) U window[16], coefficients[2][9] = {};
		for (int i = 0; i < 16; i += 8) {
			I8 x = load8(cur + i);
			x = _mm256_min_epu32(x, _mm256_sub_epi32(x, pp));
			store8(window + i, _mm256_min_epu32(x, _mm256_sub_epi32(x, p)));
		}
		for (int i = 0; i < std::min(dl, 9); ++i) {
			U x = ql[i];
			if (x >= twice) x -= twice;
			if (x >= mod) x -= mod;
			coefficients[0][i] = x;
		}
		for (int i = 0; i < std::min(dr, 9); ++i) {
			U x = qr[i];
			if (x >= twice) x -= twice;
			if (x >= mod) x -= mod;
			coefficients[1][i] = x;
		}
		I8 e[2] = {}, o[2] = {};
		for (int i = 0; i < 9; ++i) {
			I8 z = load8(window + 8 - i), high = _mm256_srli_epi64(z, 32);
			for (int k = 0; k < 2; ++k) {
				I8 r = _mm256_set1_epi32(coefficients[k][i]);
				e[k] = _mm256_add_epi64(e[k], _mm256_mul_epu32(z, r));
				o[k] = _mm256_add_epi64(o[k], _mm256_mul_epu32(high, r));
			}
		}
		U *dst[2] = {left, right};
		for (int k = 0; k < 2; ++k) {
			e[k] = _mm256_add_epi64(e[k], _mm256_mul_epu32(_mm256_mul_epu32(e[k], ni), p));
			o[k] = _mm256_add_epi64(o[k], _mm256_mul_epu32(_mm256_mul_epu32(o[k], ni), p));
			I8 x = _mm256_blend_epi32(_mm256_srli_epi64(e[k], 32), o[k], 0xaa);
			store8(dst[k], _mm256_min_epu32(x, _mm256_sub_epi32(x, pp)));
		}
	}
	__attribute__((target("avx2"))) static void inverse_sequence8(U *dst, int first, int n) {

		alignas(32) U initial[8];
		for (int j = 0; j < 8; ++j) initial[j] = in(first + j);
		I8 factor = load8(initial), step = _mm256_set1_epi32(in(8));
		I8 p = factor;
		store8(dst, p);
		for (int i = 8; i < n; i += 8) {
			factor = add8(factor, step);
			p = mul8(p, factor);
			store8(dst + i, p);
		}
		I8 inverse = _mm256_set1_epi32(one);
		for (U e = mod - 2; e; e >>= 1, p = mul8(p, p))
			if (e & 1) inverse = mul8(inverse, p);
		for (int i = n - 8; i >= 8; i -= 8) {
			store8(dst + i, mul8(inverse, load8(dst + i - 8)));
			inverse = mul8(inverse, factor);
			factor = sub8(factor, step);
		}
		store8(dst, inverse);
	}
	__attribute__((target("avx2"))) static void normalize8(U *a, int n) {
		I8 p = _mm256_set1_epi32(mod);
		int i = 0;
		for (; i + 8 <= n; i += 8) {
			I8 x = load8(a + i);
			store8(a + i, _mm256_min_epu32(x, _mm256_sub_epi32(x, p)));
		}
		for (; i < n; ++i)
			if (a[i] >= mod) a[i] -= mod;
	}
#endif
	static void weighted(U *dst, const U *src, int n, int first) {
#if defined(__x86_64__) || defined(__i386__)
		if (n >= 8 && vectorized()) {
			weighted8(dst, src, n, first);
			return;
		}
#endif
		for (int i = 0; i < n; ++i) dst[i] = mul(src[i], in(first + i));
	}
	static void scale_normalize(U *a, int n, U r) {
#if defined(__x86_64__) || defined(__i386__)
		if (n >= 8 && vectorized()) {
			scale_normalize8(a, n, r);
			return;
		}
#endif
		U qr = quotient(r);
		int i = 0;
#if defined(__ARM_NEON)
		I v = vdupq_n_u32(r), q = vdupq_n_u32(qr), p = vdupq_n_u32(mod);
		for (; i + 4 <= n; i += 4) {
			I x = fixed_mul4(vld1q_u32(a + i), v, q);
			vst1q_u32(a + i, vminq_u32(x, vsubq_u32(x, p)));
		}
#endif
		for (; i < n; ++i) {
			U x = fixed_mul(a[i], r, qr);
			a[i] = x < mod ? x : x - mod;
		}
	}
	static void square_fold(U *a, const U *q, int n, bool negative) {
#if defined(__x86_64__) || defined(__i386__)
		if (n >= 8 && vectorized()) {
			square_fold8(a, q, n, negative);
			return;
		}
#endif
		int i = 0;
#if defined(__ARM_NEON)
		for (; i + 4 <= n; i += 4) {
			I x = vld1q_u32(a + i), y = vld1q_u32(a + n + i);
			I u = mul4(vld1q_u32(q + i), vld1q_u32(q + i));
			I v = mul4(vld1q_u32(q + n + i), vld1q_u32(q + n + i));
			x = add4(x, x);
			vst1q_u32(a + i, negative ? sub4(x, u) : add4(x, u));
			vst1q_u32(a + n + i, add4(add4(y, y), v));
		}
#endif
		for (; i < n; ++i) {
			U u = mul(q[i], q[i]), v = mul(q[n + i], q[n + i]);
			U x = add(a[i], a[i]);
			a[i] = negative ? sub(x, u) : add(x, u);
			a[n + i] = add(add(a[n + i], a[n + i]), v);
		}
	}
	static void normalize(U *a, int n) {
#if defined(__x86_64__) || defined(__i386__)
		if (n >= 8 && vectorized()) {
			normalize8(a, n);
			return;
		}
#endif
		int i = 0;
#if defined(__ARM_NEON)
		for (; i + 4 <= n; i += 4) {
			I x = vld1q_u32(a + i);
			vst1q_u32(a + i, vminq_u32(x, vsubq_u32(x, vdupq_n_u32(mod))));
		}
#endif
		for (; i < n; ++i)
			if (a[i] >= mod) a[i] -= mod;
	}
	static void product(U *dst, const U *a, const U *b, int n) {
#if defined(__x86_64__) || defined(__i386__)
		if (n >= 8 && vectorized()) {
			product8(dst, a, b, n);
			return;
		}
#endif
		int i = 0;
#if defined(__ARM_NEON)
		for (; i + 4 <= n; i += 4) {
			I x = mul4(vld1q_u32(a + i), vld1q_u32(b + i));
			vst1q_u32(dst + i, x);
		}
#endif
		for (; i < n; ++i) dst[i] = mul(a[i], b[i]);
	}
#if defined(__x86_64__) || defined(__i386__)
	template <bool inverse>
	__attribute__((target("avx2"), always_inline)) void convolution_stage8(U *a, int stride, int k) const {
		const V &w = inverse ? inverse_roots : roots;
		const V &q = inverse ? inverse_root_quotients : root_quotients;
		I8 r = _mm256_set1_epi32(w[k]), rq = _mm256_set1_epi32(q[k]);
		for (int i = 0; i < stride; i += 8) {
			I8 x = load8(a + i), y = load8(a + stride + i);
			fixed_butterfly8<inverse, true>(x, y, r, rq, k == 0);
			store8(a + i, x);
			store8(a + stride + i, y);
		}
	}
	__attribute__((target("avx2"), always_inline)) void leaf8(U *a, const U *b, int k) const {
		U w = roots[k / 2], q = root_quotients[k / 2];
		if (k & 1) w = mod - w, q = ~q;
		I8 p = _mm256_set1_epi32(mod), pp = _mm256_set1_epi32(twice);
		I8 x = load8(a), y = load8(b);
		x = _mm256_min_epu32(x, _mm256_sub_epi32(x, pp));
		y = _mm256_min_epu32(y, _mm256_sub_epi32(y, pp));
		x = _mm256_min_epu32(x, _mm256_sub_epi32(x, p));
		y = _mm256_min_epu32(y, _mm256_sub_epi32(y, p));
		alignas(32) U window[16], coefficients[8];
		I8 low = fixed_mul8(x, _mm256_set1_epi32(w), _mm256_set1_epi32(q));
		store8(window, _mm256_min_epu32(low, _mm256_sub_epi32(low, p)));
		store8(window + 8, x);
		store8(coefficients, y);
		I8 e = _mm256_setzero_si256(), o = e;
		for (int i = 0; i < 8; ++i) {
			I8 z = load8(window + 8 - i), r = _mm256_set1_epi32(coefficients[i]);
			e = _mm256_add_epi64(e, _mm256_mul_epu32(z, r));
			o = _mm256_add_epi64(o, _mm256_mul_epu32(_mm256_srli_epi64(z, 32), r));
		}
		I8 ni = _mm256_set1_epi32(998244351);
		e = _mm256_add_epi64(e, _mm256_mul_epu32(_mm256_mul_epu32(e, ni), p));
		o = _mm256_add_epi64(o, _mm256_mul_epu32(_mm256_mul_epu32(o, ni), p));
		x = _mm256_blend_epi32(_mm256_srli_epi64(e, 32), o, 0xaa);
		x = _mm256_min_epu32(x, _mm256_sub_epi32(x, pp));
		store8(a, _mm256_min_epu32(x, _mm256_sub_epi32(x, p)));
	}
	template <int count, bool second = false>
	__attribute__((target("avx2"), always_inline)) void inverse_leaf_products8(U *a, U *b, const U *ws, const U *qs) const {
		I8 p = _mm256_set1_epi32(mod), pp = _mm256_set1_epi32(twice);
		alignas(32) U window[count][16], coefficients[count][8];
		I8 e[count], o[count];
		for (int k = 0; k < count; ++k) {
			I8 x = load8(a + 8 * k), y = load8(b + 8 * k);
			if constexpr (!second) x = _mm256_min_epu32(x, _mm256_sub_epi32(x, pp));
			if constexpr (!second) y = _mm256_min_epu32(y, _mm256_sub_epi32(y, pp));
			x = _mm256_min_epu32(x, _mm256_sub_epi32(x, p));
			if constexpr (!second) y = _mm256_min_epu32(y, _mm256_sub_epi32(y, p));
			I8 low = fixed_mul8(x, _mm256_set1_epi32(ws[k]), _mm256_set1_epi32(qs[k]));
			store8(window[k], _mm256_min_epu32(low, _mm256_sub_epi32(low, p)));
			store8(window[k] + 8, x);
			store8(coefficients[k], y);
			if constexpr (!second) store8(b + 8 * k, y);
			e[k] = o[k] = _mm256_setzero_si256();
		}
		for (int i = 0; i < 8; ++i)
			for (int k = 0; k < count; ++k) {
				I8 z = load8(window[k] + 8 - i), r = _mm256_set1_epi32(coefficients[k][i]);
				e[k] = _mm256_add_epi64(e[k], _mm256_mul_epu32(z, r));
				o[k] = _mm256_add_epi64(o[k], _mm256_mul_epu32(_mm256_srli_epi64(z, 32), r));
			}
		I8 ni = _mm256_set1_epi32(998244351);
		for (int k = 0; k < count; ++k) {
			I8 x = e[k], y = o[k];
			x = _mm256_add_epi64(x, _mm256_mul_epu32(_mm256_mul_epu32(x, ni), p));
			y = _mm256_add_epi64(y, _mm256_mul_epu32(_mm256_mul_epu32(y, ni), p));
			x = _mm256_blend_epi32(_mm256_srli_epi64(x, 32), y, 0xaa);
			store8(a + 8 * k, _mm256_min_epu32(x, _mm256_sub_epi32(x, pp)));
		}
	}
	__attribute__((target("avx2"))) static U reverse_dot8(const U *a, const U *b, int n) {
		I8 p = _mm256_set1_epi32(mod), pp = _mm256_set1_epi32(twice), ni = _mm256_set1_epi32(998244351);
		I8 reversed = _mm256_setr_epi32(7, 6, 5, 4, 3, 2, 1, 0), total = _mm256_setzero_si256();
		int i = 0;
		for (; i + 8 <= n;) {
			I8 even = _mm256_setzero_si256(), odd = even;
			int end = std::min(n / 8 * 8, i + 64);
			for (; i < end; i += 8) {
				I8 x = load8(a + i), y = _mm256_permutevar8x32_epi32(load8(b - i - 7), reversed);
				x = _mm256_min_epu32(x, _mm256_sub_epi32(x, pp));
				y = _mm256_min_epu32(y, _mm256_sub_epi32(y, pp));
				x = _mm256_min_epu32(x, _mm256_sub_epi32(x, p));
				y = _mm256_min_epu32(y, _mm256_sub_epi32(y, p));
				even = _mm256_add_epi64(even, _mm256_mul_epu32(x, y));
				odd = _mm256_add_epi64(odd, _mm256_mul_epu32(_mm256_srli_epi64(x, 32), _mm256_srli_epi64(y, 32)));
			}
			even = _mm256_add_epi64(even, _mm256_mul_epu32(_mm256_mul_epu32(even, ni), p));
			odd = _mm256_add_epi64(odd, _mm256_mul_epu32(_mm256_mul_epu32(odd, ni), p));
			I8 value = _mm256_blend_epi32(_mm256_srli_epi64(even, 32), odd, 0xaa);
			value = _mm256_min_epu32(value, _mm256_sub_epi32(value, pp));
			total = add8(total, value);
		}
		total = add8(total, _mm256_permute2x128_si256(total, total, 1));
		total = add8(total, _mm256_shuffle_epi32(total, 0x4e));
		total = add8(total, _mm256_shuffle_epi32(total, 0xb1));
		U result = _mm_cvtsi128_si32(_mm256_castsi256_si128(total));
		for (; i < n; ++i) {
			U x = a[i], y = b[-i];
			if (x >= twice) x -= twice;
			if (y >= twice) y -= twice;
			if (x >= mod) x -= mod;
			if (y >= mod) y -= mod;
			result = add(result, mul(x, y));
		}
		return result;
	}

	__attribute__((target("avx2"))) static void inverse_fold8(U *a, const U *f, int n) {
		int i = 0;
		for (; i + 8 <= n; i += 8) store8(a + i, add8(load8(a + i), load8(f + i)));
		for (; i < n; ++i) {
			a[i] += f[i];
			if (a[i] >= twice) a[i] -= twice;
		}
	}
	__attribute__((target("avx2"))) static void inverse_recover8(U *g, const U *a, int m) {
		int n = 2 * m;
		U ratio = T(3).pow((mod - 1) / (2 * (n / 8))).inv().val();
		U first = mod - (mod - 1) / (2 * (n / 8));
		U phase = T(3).pow((mod - 1) / 4).inv().val();
		U powers[8];
		powers[0] = first;
		for (int j = 1; j < 8; ++j) powers[j] = (unsigned long long)powers[j - 1] * ratio % mod;
		U step = T(ratio).pow(8).val();
		I8 factors = load8(powers), s = _mm256_set1_epi32(step), sq = _mm256_set1_epi32(quotient(step));
		I8 vphase = _mm256_set1_epi32(phase), qphase = _mm256_set1_epi32(quotient(phase)), p = _mm256_set1_epi32(mod), z = _mm256_setzero_si256();
		for (int i = 0; i < m; i += 64) {
			store8(powers, factors);
			for (int j = 0; j < 8; ++j) {
				int at = i + 8 * j;
				U w = powers[j];
				I8 v = _mm256_set1_epi32(w), q = _mm256_set1_epi32(quotient(w));
				I8 r = _mm256_add_epi32(load8(a + at), fixed_mul8(load8(a + m + at), vphase, qphase));
				r = fixed_mul8(r, v, q);
				store8(g + m + at, sub8(z, add8(r, load8(g + m + at))));
			}
			factors = fixed_mul8(factors, s, sq);
			factors = _mm256_min_epu32(factors, _mm256_sub_epi32(factors, p));
		}
	}
	template <bool inverse>
	__attribute__((target("avx2"))) void inverse_partial8(U *a, int n, int limit = 0) const {
		if (!limit) limit = n;
		if constexpr (!inverse) {
			int stride = limit / 2;
			for (; stride >= 16; stride >>= 2) radix4_8<false>(a, stride / 2, n);
			if (stride == 8)
				for (int i = 0; i < n; i += 16) convolution_stage8<false>(a + i, 8, i / 16);
		} else {
			int stride = 8;
			for (; 4 * stride <= n; stride <<= 2) radix4_8<true>(a, stride, n);
			if (stride < n) stage8<true>(a, stride, n);
		}
	}
	__attribute__((target("avx2"))) static void inverse_twist8(U *a, U *b, int n, U ratio, U first = one) {
		U powers[8];
		powers[0] = out(first);
		U r = out(ratio);
		for (int j = 1; j < 8; ++j) powers[j] = (unsigned long long)powers[j - 1] * r % mod;
		U step = T::word(ratio).pow(8).val();
		I8 factors = load8(powers), s = _mm256_set1_epi32(step), sq = _mm256_set1_epi32(quotient(step)), p = _mm256_set1_epi32(mod);
		for (int i = 0; i < n; i += 64) {
			store8(powers, factors);
			for (int j = 0; j < 8 && i + 8 * j < n; ++j) {
				U w = powers[j];
				I8 v = _mm256_set1_epi32(w), q = _mm256_set1_epi32(quotient(w));
				store8(a + i + 8 * j, fixed_mul8(load8(a + i + 8 * j), v, q));
				if (b && i + 8 * j < n / 2) store8(b + i + 8 * j, fixed_mul8(load8(b + i + 8 * j), v, q));
			}
			factors = fixed_mul8(factors, s, sq);
			factors = _mm256_min_epu32(factors, _mm256_sub_epi32(factors, p));
		}
	}
	template <bool negative>
	__attribute__((target("avx2"))) void inverse_cubic8(U *a, U *b, int n) const {
		U shift = one;
		if constexpr (negative) {
			T psi = T(3).pow((mod - 1) / (2 * (n / 8)));
			shift = psi.v;
			if (shift >= mod) shift -= mod;
			inverse_twist8(a, b, n, shift);
		}
		inverse_partial8<false>(a, n);
		if constexpr (negative) {
			std::copy_n(b, n / 2, b + n / 2);
			inverse_partial8<false>(b, n, n / 2);
		} else inverse_partial8<false>(b, n);
		for (int i = 0; i < n; i += 8 * 8) {
			U ws[8], qs[8];
			for (int j = 0; j < 8; ++j) {
				int k = i / 8 + j;
				U w = roots[k / 2];
				if (k & 1) w = mod - w;
				if constexpr (negative) {
					w = mul(w, shift);
					if (w >= mod) w -= mod;
					qs[j] = quotient(w);
				} else qs[j] = (k & 1) ? ~root_quotients[k / 2] : root_quotients[k / 2];
				ws[j] = w;
			}
			inverse_leaf_products8<8>(a + i, b + i, ws, qs);
			inverse_leaf_products8<8, true>(a + i, b + i, ws, qs);
			for (int j = 0; j < 8; ++j) store8(a + i + 8 * j, sub8(load8(a + i + 8 * j), load8(b + i + 8 * j)));
		}
		inverse_partial8<true>(a, n);
		U normal = mod - (mod - 1) / (2 * (n / 8));
		if constexpr (!negative) scale8(a, n, normal);
	}
	template <int n>
	__attribute__((target("avx2"), always_inline)) void convolution_tile8(U *a, U *b, int k) const {
		if constexpr (n == 8) leaf8(a, b, k);
		else if constexpr (n == 16) {
			convolution_stage8<false>(a, n / 2, k);
			convolution_stage8<false>(b, n / 2, k);
			convolution_tile8<n / 2>(a, b, 2 * k);
			convolution_tile8<n / 2>(a + n / 2, b + n / 2, 2 * k + 1);
			convolution_stage8<true>(a, n / 2, k);
		} else {
			stage_pair8<false, true>(a, n / 4, n, k);
			stage_pair8<false, true>(b, n / 4, n, k);
			for (int i = 0; i < 4; ++i) convolution_tile8<n / 4>(a + i * (n / 4), b + i * (n / 4), 4 * k + i);
			stage_pair8<true>(a, n / 4, n, k);
		}
	}
	__attribute__((target("avx2"))) void convolution8(U *a, U *b, int n, int k = 0) const {
		if (n == 256) {
			convolution_tile8<256>(a, b, k);
			return;
		}
		if (n == 128) {
			convolution_tile8<128>(a, b, k);
			return;
		}
		if (n == 64) {
			convolution_tile8<64>(a, b, k);
			return;
		}
		stage_pair8<false, true>(a, n / 4, n, k);
		stage_pair8<false, true>(b, n / 4, n, k);
		for (int i = 0; i < 4; ++i) convolution8(a + i * (n / 4), b + i * (n / 4), n / 4, 4 * k + i);
		stage_pair8<true>(a, n / 4, n, k);
	}
#endif
#if defined(__x86_64__) || defined(__i386__)
	template <bool odd>
	__attribute__((target("avx2"))) static void reduce_pairs8(U *a, U *b, const U *twist, int n) {
		for (int i = 0; i < n; i += 8) {
			I8 x, y, u, v;
			load_pairs(a + 2 * i, x, y);
			load_pairs(b + 2 * i, u, v);
			x = mul8(x, v);
			y = mul8(y, u);
			if constexpr (odd) x = mul8(sub8(x, y), load8(twist + i));
			else x = add8(x, y);
			store8(a + i, x);
			store8(b + i, mul8(u, v));
		}
	}
	__attribute__((target("avx2"))) static void combine8(U *a, const U *b, const U *x, const U *y, int n) {
		for (int i = 0; i < n; i += 8) store8(a + i, add8(mul8(load8(a + i), load8(y + i)), mul8(load8(b + i), load8(x + i))));
	}
#endif

	template <bool odd>
	static void reduce_pairs(U *a, U *b, const U *twist, int n) {
#if defined(__x86_64__) || defined(__i386__)
		if (n >= 8 && vectorized()) {
			reduce_pairs8<odd>(a, b, twist, n);
			return;
		}
#endif
		int i = 0;
#if defined(__ARM_NEON)
		for (; i + 4 <= n; i += 4) {
			auto p = vld2q_u32(a + 2 * i), q = vld2q_u32(b + 2 * i);
			I x = mul4(p.val[0], q.val[1]), y = mul4(p.val[1], q.val[0]);
			if constexpr (odd) x = mul4(sub4(x, y), vld1q_u32(twist + i));
			else x = add4(x, y);
			vst1q_u32(a + i, x);
			vst1q_u32(b + i, mul4(q.val[0], q.val[1]));
		}
#endif
		for (; i < n; ++i) {
			U u = mul(a[2 * i], b[2 * i + 1]), v = mul(a[2 * i + 1], b[2 * i]);
			a[i] = odd ? mul(sub(u, v), twist[i]) : add(u, v);
			b[i] = mul(b[2 * i], b[2 * i + 1]);
		}
	}
	static void combine(U *a, const U *b, const U *x, const U *y, int n) {
#if defined(__x86_64__) || defined(__i386__)
		if (n >= 8 && vectorized()) {
			combine8(a, b, x, y, n);
			return;
		}
#endif
		int i = 0;
#if defined(__ARM_NEON)
		for (; i + 4 <= n; i += 4) vst1q_u32(a + i, add4(mul4(vld1q_u32(a + i), vld1q_u32(y + i)), mul4(vld1q_u32(b + i), vld1q_u32(x + i))));
#endif
		for (; i < n; ++i) a[i] = add(mul(a[i], y[i]), mul(b[i], x[i]));
	}
#if defined(__x86_64__) || defined(__i386__)
	__attribute__((target("avx2"))) static void geometric8(U *a, int n, U first, U step) {
		alignas(32) U factors[8];
		U power = one;
		for (int i = 0; i < 8; ++i) {
			factors[i] = mul(first, power);
			power = mul(power, step);
		}
		I8 r = load8(factors), ratio = _mm256_set1_epi32(power);
		int i = 0;
		for (; i + 8 <= n; i += 8) {
			store8(a + i, mul8(load8(a + i), r));
			r = mul8(r, ratio);
		}
		store8(factors, r);
		for (int j = 0; i < n; ++i, ++j) a[i] = mul(a[i], factors[j]);
	}
#endif
	static void geometric(U *a, int n, U first, U step) {
#if defined(__x86_64__) || defined(__i386__)
		if (n >= 8 && vectorized()) {
			geometric8(a, n, first, step);
			return;
		}
#endif
		int i = 0;
#if defined(__ARM_NEON)
		U factors[4], power = one;
		for (int j = 0; j < 4; ++j) {
			factors[j] = mul(first, power);
			power = mul(power, step);
		}
		I r = vld1q_u32(factors), ratio = vdupq_n_u32(power);
		for (; i + 4 <= n; i += 4) {
			vst1q_u32(a + i, mul4(vld1q_u32(a + i), r));
			r = mul4(r, ratio);
		}
		vst1q_u32(factors, r);
		first = factors[0];
#endif
		for (; i < n; ++i) {
			a[i] = mul(a[i], first);
			first = mul(first, step);
		}
	}
#if defined(__x86_64__) || defined(__i386__)
	__attribute__((target("avx2"))) static void add_scaled8(U *a, const U *b, int n, U r) {
		U q = quotient(r);
		I8 v = _mm256_set1_epi32(r), vq = _mm256_set1_epi32(q);
		int i = 0;
		for (; i + 8 <= n; i += 8) store8(a + i, add8(load8(a + i), fixed_mul8(load8(b + i), v, vq)));
		for (; i < n; ++i) a[i] = add(a[i], fixed_mul(b[i], r, q));
	}
#endif
	static void add_scaled(U *a, const U *b, int n, U r) {
#if defined(__x86_64__) || defined(__i386__)
		if (n >= 8 && vectorized()) {
			add_scaled8(a, b, n, r);
			return;
		}
#endif
		int i = 0;
		U q = quotient(r);
#if defined(__ARM_NEON)
		for (; i + 4 <= n; i += 4) vst1q_u32(a + i, add4(vld1q_u32(a + i), fixed_mul4(vld1q_u32(b + i), vdupq_n_u32(r), vdupq_n_u32(q))));
#endif
		for (; i < n; ++i) a[i] = add(a[i], fixed_mul(b[i], r, q));
	}
	void leaf(U *a, const U *b, int k) const {
		U w = roots[k / 2], q = root_quotients[k / 2];
		if (k & 1) w = mod - w, q = ~q;
		U window[16], coefficients[8];
		for (int i = 0; i < 8; ++i) {
			U x = a[i] < mod ? a[i] : a[i] - mod;
			U z = fixed_mul(x, w, q);
			window[i] = z < mod ? z : z - mod;
			window[8 + i] = x;
			coefficients[i] = b[i] < mod ? b[i] : b[i] - mod;
		}
#if defined(__ARM_NEON)
		for (int j = 0; j < 8; j += 4) {
			uint64x2_t lo = vdupq_n_u64(0), hi = lo;
			for (int i = 0; i < 8; ++i) {
				I z = vld1q_u32(window + 8 - i + j);
				lo = vmlal_n_u32(lo, vget_low_u32(z), coefficients[i]);
				hi = vmlal_n_u32(hi, vget_high_u32(z), coefficients[i]);
			}
			uint32x2_t x = vmul_n_u32(vmovn_u64(lo), 998244351u), y = vmul_n_u32(vmovn_u64(hi), 998244351u);
			I z = vcombine_u32(vshrn_n_u64(vmlal_n_u32(lo, x, mod), 32), vshrn_n_u64(vmlal_n_u32(hi, y, mod), 32));
			z = vminq_u32(z, vsubq_u32(z, vdupq_n_u32(twice)));
			vst1q_u32(a + j, vminq_u32(z, vsubq_u32(z, vdupq_n_u32(mod))));
		}
#else
		for (int j = 0; j < 8; ++j) {
			unsigned long long sum = 0;
			for (int i = 0; i < 8; ++i) sum += (unsigned long long)window[8 - i + j] * coefficients[i];
			U x = reduce(sum);
			if (x >= twice) x -= twice;
			a[j] = x < mod ? x : x - mod;
		}
#endif
	}
	template <int n>
	void convolution_tile(U *a, U *b, int k) const {
		if constexpr (n == 8) leaf(a, b, k);
		else if constexpr (n == 16) {
			stage<false>(a, n / 2, n, k);
			stage<false>(b, n / 2, n, k);
			convolution_tile<n / 2>(a, b, 2 * k);
			convolution_tile<n / 2>(a + n / 2, b + n / 2, 2 * k + 1);
			stage<true>(a, n / 2, n, k);
		} else {
			stage_pair<false>(a, n / 4, n, k);
			stage_pair<false>(b, n / 4, n, k);
			for (int i = 0; i < 4; ++i) convolution_tile<n / 4>(a + i * (n / 4), b + i * (n / 4), 4 * k + i);
			stage_pair<true>(a, n / 4, n, k);
		}
	}
	void convolution(U *a, U *b, int n, int k = 0) const {
#if defined(__x86_64__) || defined(__i386__)
		if (vectorized()) {
			convolution8(a, b, n, k);
			return;
		}
#endif
		if (n == 256) {
			convolution_tile<256>(a, b, k);
			return;
		}
		if (n == 128) {
			convolution_tile<128>(a, b, k);
			return;
		}
		if (n == 64) {
			convolution_tile<64>(a, b, k);
			return;
		}
		stage_pair<false>(a, n / 4, n, k);
		stage_pair<false>(b, n / 4, n, k);
		for (int i = 0; i < 4; ++i) convolution(a + i * (n / 4), b + i * (n / 4), n / 4, 4 * k + i);
		stage_pair<true>(a, n / 4, n, k);
	}
#if defined(__x86_64__) || defined(__i386__)
	__attribute__((target("avx2"))) static void accumulate8(U *dst, const U *a, const U *b, const U *c, int n, bool minus) {
		for (int i = 0; i < n; i += 8) {
			I8 x = load8(b + i), y = load8(c + i);
			x = minus ? _mm256_sub_epi32(_mm256_add_epi32(x, _mm256_set1_epi32(twice)), y)
					  : _mm256_add_epi32(x, y);
			store8(dst + i, add8(load8(dst + i), mul8(load8(a + i), x)));
		}
	}
#endif
#if defined(__x86_64__) || defined(__i386__)
	template <int count>
	__attribute__((target("avx2"))) static void accumulate_batch8(U *dst, const U *const *a, const U *const *b,
																  const U *const *c, int n, bool minus) {
		static_assert(count >= 1 && count <= 3);
		I8 p = _mm256_set1_epi32(mod), pp = _mm256_set1_epi32(twice), ni = _mm256_set1_epi32(998244351);
		int offset = minus ? n : 0;
		for (int i = 0; i < n; i += 8) {
			I8 even = _mm256_setzero_si256(), odd = even;
#pragma GCC unroll 3
			for (int k = 0; k < count; ++k) {
				I8 x = load8(b[k] + offset + i), y = load8(c[k] + offset + i);
				x = minus ? _mm256_sub_epi32(_mm256_add_epi32(x, pp), y) : _mm256_add_epi32(x, y);
				I8 z = load8(a[k] + offset + i);
				even = _mm256_add_epi64(even, _mm256_mul_epu32(z, x));
				odd = _mm256_add_epi64(odd, _mm256_mul_epu32(_mm256_srli_epi64(z, 32), _mm256_srli_epi64(x, 32)));
			}
			even = _mm256_add_epi64(even, _mm256_mul_epu32(_mm256_mul_epu32(even, ni), p));
			odd = _mm256_add_epi64(odd, _mm256_mul_epu32(_mm256_mul_epu32(odd, ni), p));
			I8 value = _mm256_blend_epi32(_mm256_srli_epi64(even, 32), odd, 0xaa);
			value = _mm256_min_epu32(value, _mm256_sub_epi32(value, pp));
			store8(dst + i, add8(load8(dst + i), value));
		}
	}
#endif
#if defined(__x86_64__) || defined(__i386__)
	template <int count, bool first, bool canonical>
	__attribute__((target("avx2"))) static void accumulate_wide8(U *dst, const U *const *a, const U *const *b, const U *const *c, int n, bool minus) {
		static_assert(count >= 1 && count <= 7);
		I8 p = _mm256_set1_epi32(mod), pp = _mm256_set1_epi32(twice), ni = _mm256_set1_epi32(998244351);
		int offset = minus ? n : 0;
		for (int i = 0; i < n; i += 8) {
			I8 even = _mm256_setzero_si256(), odd = even;
#pragma GCC unroll 7
			for (int k = 0; k < count; ++k) {
				I8 x = load8(b[k] + offset + i), y = load8(c[k] + offset + i);
				x = minus ? _mm256_sub_epi32(_mm256_add_epi32(x, canonical ? p : pp), y) : _mm256_add_epi32(x, y);
				if constexpr (!canonical) x = _mm256_min_epu32(x, _mm256_sub_epi32(x, pp));
				I8 z = load8(a[k] + offset + i);
				even = _mm256_add_epi64(even, _mm256_mul_epu32(z, x));
				odd = _mm256_add_epi64(odd, _mm256_mul_epu32(_mm256_srli_epi64(z, 32), _mm256_srli_epi64(x, 32)));
			}
			even = _mm256_add_epi64(even, _mm256_mul_epu32(_mm256_mul_epu32(even, ni), p));
			odd = _mm256_add_epi64(odd, _mm256_mul_epu32(_mm256_mul_epu32(odd, ni), p));
			I8 value = _mm256_blend_epi32(_mm256_srli_epi64(even, 32), odd, 0xaa);
			value = _mm256_min_epu32(value, _mm256_sub_epi32(value, pp));
			if constexpr (count == 7) value = _mm256_min_epu32(value, _mm256_sub_epi32(value, pp));
			if constexpr (!first) value = add8(load8(dst + i), value);
			store8(dst + i, value);
		}
	}
#endif
	template <int count, bool first, bool canonical>
	static void accumulate_wide(U *dst, const U *const *a, const U *const *b, const U *const *c, int n, bool minus) {
#if defined(__x86_64__) || defined(__i386__)
		if (n >= 8 && vectorized()) {
			accumulate_wide8<count, first, canonical>(dst, a, b, c, n, minus);
			return;
		}
#endif
		if constexpr (first) std::fill_n(dst, n, 0u);
		int offset = minus ? n : 0;
		for (int k = 0; k < count; ++k) accumulate(dst, a[k] + offset, b[k] + offset, c[k] + offset, n, minus);
	}
	template <int count>
	static void accumulate_batch(U *dst, const U *const *a, const U *const *b, const U *const *c, int n, bool minus) {
#if defined(__x86_64__) || defined(__i386__)
		if (n >= 8 && vectorized()) {
			accumulate_batch8<count>(dst, a, b, c, n, minus);
			return;
		}
#endif
		int offset = minus ? n : 0;
		for (int k = 0; k < count; ++k) accumulate(dst, a[k] + offset, b[k] + offset, c[k] + offset, n, minus);
	}
	static void accumulate(U *dst, const U *a, const U *b, const U *c, int n, bool minus) {
#if defined(__x86_64__) || defined(__i386__)
		if (n >= 8 && vectorized()) {
			accumulate8(dst, a, b, c, n, minus);
			return;
		}
#endif
		int i = 0;
#if defined(__ARM_NEON)
		for (; i + 4 <= n; i += 4) {
			I x = vld1q_u32(b + i), y = vld1q_u32(c + i);
			x = minus ? sub4(x, y) : add4(x, y);
			vst1q_u32(dst + i, add4(vld1q_u32(dst + i), mul4(vld1q_u32(a + i), x)));
		}
#endif
		for (; i < n; ++i) dst[i] = add(dst[i], mul(a[i], minus ? sub(b[i], c[i]) : add(b[i], c[i])));
	}
	static void scale_copy(U *dst, const U *src, int length, U r) {
#if defined(__x86_64__) || defined(__i386__)
		if (length >= 8 && vectorized()) {
			scale_copy8(dst, src, length, r);
			return;
		}
#endif
		U q = quotient(r);
		int i = 0;
#if defined(__ARM_NEON)
		for (; i + 4 <= length; i += 4) vst1q_u32(dst + i, fixed_mul4(vld1q_u32(src + i), vdupq_n_u32(r), vdupq_n_u32(q)));
#endif
		for (; i < length; ++i) dst[i] = fixed_mul(src[i], r, q);
	}
	void scale(U *a, int length, U r) const {
#if defined(__x86_64__) || defined(__i386__)
		if (length >= 8 && vectorized()) {
			scale8(a, length, r);
			return;
		}
#endif
		int i = 0;
		U q = quotient(r);
#if defined(__ARM_NEON)
		for (; i + 4 <= length; i += 4)
			vst1q_u32(a + i, fixed_mul4(vld1q_u32(a + i), vdupq_n_u32(r), vdupq_n_u32(q)));
#endif
		for (; i < length; ++i) a[i] = fixed_mul(a[i], r, q);
	}

	static void pair_product(U *a, const U *b, int length) {
#if defined(__x86_64__) || defined(__i386__)
		if (length >= 8 && vectorized()) {
			pair_product8(a, b, length);
			return;
		}
#endif
		int i = 0;
#if defined(__ARM_NEON)
		for (; i + 4 <= length; i += 4) {
			uint32x4x2_t v = vld2q_u32(b + 2 * i);
			vst1q_u32(a + i, mul4(v.val[0], v.val[1]));
		}
#endif
		for (; i < length; ++i) a[i] = mul(b[2 * i], b[2 * i + 1]);
	}
	static void pair_multiply(U *a, const U *b, int length) {
#if defined(__x86_64__) || defined(__i386__)
		if (length >= 8 && vectorized()) {
			pair_multiply8(a, b, length);
			return;
		}
#endif
		int i = 0;
#if defined(__ARM_NEON)
		for (; i + 4 <= length; i += 4) {
			uint32x4x2_t v = vld2q_u32(a + 2 * i);
			I x = mul4(v.val[1], vld1q_u32(b + i));
			v.val[1] = mul4(v.val[0], vld1q_u32(b + i));
			v.val[0] = x;
			vst2q_u32(a + 2 * i, v);
		}
#endif
		for (; i < length; ++i) {
			U x = a[2 * i];
			a[2 * i] = mul(b[i], a[2 * i + 1]);
			a[2 * i + 1] = mul(b[i], x);
		}
	}
	static void difference(U *a, const U *b, const U *c, int length) {
#if defined(__x86_64__) || defined(__i386__)
		if (length >= 8 && vectorized()) {
			difference8(a, b, c, length);
			return;
		}
#endif
		int i = 0;
#if defined(__ARM_NEON)
		for (; i + 4 <= length; i += 4)
			vst1q_u32(a + i, sub4(vld1q_u32(b + i), vld1q_u32(c + i)));
#endif
		for (; i < length; ++i) a[i] = sub(b[i], c[i]);
	}
	static void fold(U *a, int columns) {
#if defined(__x86_64__) || defined(__i386__)
		if (columns >= 8 && vectorized()) {
			fold8(a, columns);
			return;
		}
#endif
		int i = 0;
#if defined(__ARM_NEON)
		for (; i + 4 <= columns; i += 4) {
			I x = add4(vld1q_u32(a + i), vld1q_u32(a + columns + i));
			vst1q_u32(a + i, x);
			vst1q_u32(a + columns + i, x);
		}
#endif
		for (; i < columns; ++i) a[i] = a[columns + i] = add(a[i], a[columns + i]);
	}
	void project(U *a, int rows, int columns) const {
		if (columns == 1) {
#if defined(__x86_64__) || defined(__i386__)
			if (rows >= 8 && vectorized()) {
				project1_8(a, rows);
				return;
			}
#endif
			int j = 0;
#if defined(__ARM_NEON)
			for (; j + 4 <= rows; j += 4) {
				uint32x4x2_t v = vld2q_u32(a + 2 * j);
				v.val[0] = v.val[1] = add4(v.val[0], v.val[1]);
				vst2q_u32(a + 2 * j, v);
			}
#endif
			for (; j < rows; ++j) a[2 * j] = a[2 * j + 1] = add(a[2 * j], a[2 * j + 1]);
			return;
		}
		if (columns == 2) {
#if defined(__x86_64__) || defined(__i386__)
			if (rows >= 8 && vectorized()) {
				project2_8(a, rows);
				return;
			}
#endif
			int j = 0;
#if defined(__ARM_NEON)
			for (; j + 4 <= rows; j += 4) {
				uint32x4x4_t v = vld4q_u32(a + 4 * j);
				I x = add4(add4(v.val[0], v.val[1]), add4(v.val[2], v.val[3]));
				I y = add4(sub4(v.val[0], v.val[1]),
						   fixed_mul4(sub4(v.val[2], v.val[3]), vdupq_n_u32(inverse_roots[1]),
									  vdupq_n_u32(inverse_root_quotients[1])));
				v.val[0] = add4(x, y);
				v.val[1] = sub4(x, y);
				y = fixed_mul4(y, vdupq_n_u32(roots[1]), vdupq_n_u32(root_quotients[1]));
				v.val[2] = add4(x, y);
				v.val[3] = sub4(x, y);
				vst4q_u32(a + 4 * j, v);
			}
#endif
			for (; j < rows; ++j) {
				U *r = a + 4 * j;
				U x = add(add(r[0], r[1]), add(r[2], r[3]));
				U y = add(sub(r[0], r[1]), fixed_mul(sub(r[2], r[3]), inverse_roots[1],
													 inverse_root_quotients[1]));
				r[0] = add(x, y);
				r[1] = sub(x, y);
				y = fixed_mul(y, roots[1], root_quotients[1]);
				r[2] = add(x, y);
				r[3] = sub(x, y);
			}
			return;
		}
		if (columns == 4) {
#if defined(__x86_64__) || defined(__i386__)
			if (vectorized()) {
				project4_8(a, rows);
				return;
			}
#endif
#if defined(__ARM_NEON)
			for (int j = 0; j < rows; ++j) {
				I x = vld1q_u32(a + 8 * j), y = vld1q_u32(a + 8 * j + 4);
				small4<1, true>(x, y, 0);
				small4<2, true>(x, y, 0);
				x = y = add4(x, y);
				small4<2, false>(x, y, 0);
				small4<1, false>(x, y, 0);
				vst1q_u32(a + 8 * j, x);
				vst1q_u32(a + 8 * j + 4, y);
			}
			return;
#endif
		}
#if defined(__x86_64__) || defined(__i386__)
		if (columns == 8 && vectorized()) {
			project8_8(a, rows);
			return;
		}
#endif
#if defined(__x86_64__) || defined(__i386__)
		if (columns == 16 && vectorized()) {
			project16_8(a, rows);
			return;
		}
#endif
		for (int j = 0; j < rows; ++j) {
			U *row = a + 2 * columns * j;
			fft<true>(row, 2 * columns, columns);
			fold(row, columns);
			fft<false>(row, 2 * columns, columns);
		}
	}
	void composition_difference(U *a, const U *b, const U *c, int n) const {
#if defined(__x86_64__) || defined(__i386__)
		if (vectorized()) {
			composition_difference8(a, b, c, n);
			return;
		}
#endif
		difference(a, b, c, n);
	}
#if defined(__x86_64__) || defined(__i386__)
	__attribute__((target("avx2"))) static void composition_difference8(U *a, const U *b, const U *c, int n) {
		I8 p = _mm256_set1_epi32(twice);
		for (int i = 0; i < n; i += 8) {
			I8 x = load8(c + i);
			x = _mm256_min_epu32(x, _mm256_sub_epi32(x, p));
			store8(a + i, sub8(load8(b + i), x));
		}
	}
#endif
	void composition_finish(U *p, const U *v, int columns, U factor) const {
#if defined(__x86_64__) || defined(__i386__)
		if (columns >= 8 && vectorized()) {
			composition_finish8(p, v, columns, out(factor));
			return;
		}
#endif
		for (int i = 0; i < columns; ++i) p[i] = mul(sub(add(v[i], v[columns + i]), add(v[2 * columns + i], v[3 * columns + i])), factor);
	}
#if defined(__x86_64__) || defined(__i386__)
	__attribute__((target("avx2"))) static void composition_finish8(U *dst, const U *v, int n, U factor) {
		I8 r = _mm256_set1_epi32(factor), q = _mm256_set1_epi32(quotient(factor));
		for (int i = 0; i < n; i += 8) {
			I8 x = add8(load8(v + i), load8(v + n + i)), y = add8(load8(v + 2 * n + i), load8(v + 3 * n + i));
			store8(dst + i, fixed_mul8(sub8(x, y), r, q));
		}
	}
#endif
	struct CompositionDelete {
		void operator()(unsigned int *p) const { ::operator delete[](p, std::align_val_t(64)); }
	};
	using CompositionWords = std::unique_ptr<unsigned int[], CompositionDelete>;
	void composition_forward_y(U *a, int size, int columns) const {
#if defined(__x86_64__) || defined(__i386__)
		if (vectorized()) {
			int h = size;
			bool lazy = false;
			for (; h / 2 > columns && h / 2 >= 8; h >>= 2) {
				radix4_8<false>(a, h / 2, 2 * size);
				lazy = true;
			}
			if (lazy && h > columns) composition_shrink8(a, 2 * size);
			for (; h / 2 > columns; h >>= 2) stage_pair<false>(a, h / 2, 2 * size);
			if (h > columns) stage<false>(a, h, 2 * size);
			return;
		}
#endif
		int h = size;
		for (; h / 2 > columns; h >>= 2) stage_pair<false>(a, h / 2, 2 * size);
		if (h > columns) stage<false>(a, h, 2 * size);
	}
#if defined(__x86_64__) || defined(__i386__)
	__attribute__((target("avx2"))) static void composition_shrink8(U *a, int n) {
		I8 p = _mm256_set1_epi32(twice);
		for (int i = 0; i < n; i += 8) {
			I8 x = load8(a + i);
			store8(a + i, _mm256_min_epu32(x, _mm256_sub_epi32(x, p)));
		}
	}
#endif
	void double_y(U *a, int rows, int columns, bool inverse) const {
		const int size = rows * columns;
		int h = 2 * columns;
		for (; h <= size / 2; h <<= 2) {
#if defined(__x86_64__) || defined(__i386__)
			if (h >= 8 && vectorized()) radix4_8<true>(a, h, 2 * size);
			else
#endif
				stage_pair<true>(a, h, 2 * size);
		}
		if (h <= size) stage<true>(a, h, 2 * size);
#if defined(__x86_64__) || defined(__i386__)
		if (columns >= 2 && vectorized()) {
			scale_rows8(a, rows, columns, (inverse ? inverse_twists : twists).data() + rows);
			return;
		}
#endif
		for (int j = 0; j < rows; ++j)
			scale(a + 2 * columns * j, 2 * columns,
				  (inverse ? inverse_twists : twists)[rows + j]);
	}
#if defined(__x86_64__) || defined(__i386__)
	__attribute__((target("avx2"))) static void projection_dot8(U *p, const U *q, int n) {
		I8 ni = _mm256_set1_epi32(998244351), modulus = _mm256_set1_epi32(mod), pp = _mm256_set1_epi32(twice);
		I8 index = _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7);
		for (int i = 0; i < n; i += 8) {
			I8 x = load8(p + 2 * i), y = load8(p + 2 * i + 8), a = load8(q + 2 * i), b = load8(q + 2 * i + 8);
			x = _mm256_min_epu32(x, _mm256_sub_epi32(x, pp));
			y = _mm256_min_epu32(y, _mm256_sub_epi32(y, pp));
			a = _mm256_min_epu32(a, _mm256_sub_epi32(a, pp));
			b = _mm256_min_epu32(b, _mm256_sub_epi32(b, pp));
			I8 lo = _mm256_add_epi64(_mm256_mul_epu32(x, _mm256_srli_epi64(a, 32)), _mm256_mul_epu32(_mm256_srli_epi64(x, 32), a));
			I8 hi = _mm256_add_epi64(_mm256_mul_epu32(y, _mm256_srli_epi64(b, 32)), _mm256_mul_epu32(_mm256_srli_epi64(y, 32), b));
			lo = _mm256_add_epi64(lo, _mm256_mul_epu32(_mm256_mul_epu32(lo, ni), modulus));
			hi = _mm256_add_epi64(hi, _mm256_mul_epu32(_mm256_mul_epu32(hi, ni), modulus));
			I8 r = _mm256_permutevar8x32_epi32(_mm256_blend_epi32(_mm256_srli_epi64(lo, 32), hi, 0xaa), index);
			store8(p + i, _mm256_min_epu32(r, _mm256_sub_epi32(r, pp)));
		}
	}
#endif

	static void projection_dot(U *p, const U *q, int n) {
#if defined(__x86_64__) || defined(__i386__)
		if (n >= 8 && vectorized()) {
			projection_dot8(p, q, n);
			return;
		}
#endif
		for (int i = 0; i < n; ++i) {
			U x = p[2 * i], y = p[2 * i + 1], a = q[2 * i], b = q[2 * i + 1];
			if (x >= twice) x -= twice;
			if (y >= twice) y -= twice;
			if (a >= twice) a -= twice;
			if (b >= twice) b -= twice;
			U r = reduce((unsigned long long)x * b + (unsigned long long)y * a);
			p[i] = r < twice ? r : r - twice;
		}
	}
	void transpose_roots() {
		roots.swap(inverse_roots);
		root_quotients.swap(inverse_root_quotients);
		cube_roots.swap(inverse_cube_roots);
		cube_quotients.swap(inverse_cube_quotients);
	}
	U reversion_projection(U *p, U *q, int size) {
		int levels = __builtin_ctz((unsigned int)size);
		V saved(levels, one);
		U q_one = one, p_one = one;
		for (int i = 1, columns = size / 2; columns > 1; ++i, columns >>= 1) {
			q_one = mul(mul(q_one, q_one), in(2 * columns));
			saved[i] = q_one;
		}
		for (int i = levels - 1, columns = 2; i >= 0; --i, columns <<= 1) {
			p_one = mul(p_one, saved[i]);
			if (i) p_one = mul(p_one, in(4 * columns));
		}
		U factor = T::word(mul(p_one, in(4 * size))).inv().v;
		fft<false>(q, 2 * size);
		for (int i = 0; i < 2 * size; ++i) {
			U x = q[i];
			q[i] = add(x, one);
			q[2 * size + i] = sub(x, one);
		}

		for (int i = 0; i < size; ++i) {
			U x = fixed_mul(sub(q[2 * i + 1], q[2 * i]), roots[i], root_quotients[i]);
			p[i] = i < size / 2 ? x : sub(0, x);
			p[size + i] = sub(0, p[i]);
		}
		q_one = one;
		for (int columns = size / 2, rows = 2; columns > 1; columns >>= 1, rows <<= 1) {
			pair_product(q, q, 2 * size);
			q_one = mul(q_one, q_one);
			project(q, rows, columns);
			std::copy_n(q, 2 * size, q + 2 * size);
			double_y(q + 2 * size, rows, columns, false);
			q_one = mul(q_one, in(2 * columns));
			U correction = add(q_one, q_one);
			for (int i = 0; i < 2 * columns; ++i) q[2 * size + i] = sub(q[2 * size + i], correction);
			composition_forward_y(q + 2 * size, size, columns);
			transpose_roots();
			project(p, rows, columns);
			for (int i = 0; i < 2 * size; ++i) p[2 * size + i] = sub(0, p[i]);
			double_y(p + 2 * size, rows, columns, true);
			composition_forward_y(p + 2 * size, size, columns);
			transpose_roots();
			projection_dot(p, q, 2 * size);
		}
		for (int i = 0; i < size; ++i) p[i] = add(p[2 * i], p[2 * i + 1]);
		transpose_roots();
		fft<true>(p, size);
		transpose_roots();
		return factor;
	}

	void compose(U *p, const U *q, int rows, int columns, U &p_one, U &q_one, U *workspace = nullptr, int level = 0) const {
		const int size = rows * columns;
		if (columns == 1) {
			fft<false>(p, rows);
			for (int i = rows - 1; i >= 0; --i) p[2 * i] = p[2 * i + 1] = p[i];
			return;
		}
		CompositionWords storage;
		if (!workspace) {
			int levels = __builtin_ctz((unsigned int)columns);
			storage.reset(static_cast<unsigned int *>(::operator new[]((4 * size + 16) * levels * sizeof(unsigned int), std::align_val_t(64))));
			workspace = storage.get();
		}
		U *v = workspace + (4 * size + 16) * level;
		if (rows == 1) {
			std::fill_n(v + columns, columns, U(0));
			for (int i = 0; i < columns; ++i) v[i] = sub(0, q[i]);
			fft<false>(v, 2 * columns);
			for (int i = 0; i < 2 * columns; ++i) {
				U x = v[i];
				v[i] = add(x, one);
				v[i + 2 * columns] = sub(x, one);
			}
		} else {
			pair_product(v, q, 2 * size);
			q_one = mul(q_one, q_one);
			project(v, rows, columns);
			std::copy_n(v, 2 * size, v + 2 * size);
			U *upper = v + 2 * size;
			double_y(upper, rows, columns, false);
			q_one = mul(q_one, in(2 * columns));
			U correction = add(q_one, q_one);
			for (int i = 0; i < 2 * columns; ++i) upper[i] = sub(upper[i], correction);
			composition_forward_y(upper, size, columns);
		}
		U saved = q_one;
		compose(p, v, 2 * rows, columns / 2, p_one, q_one, workspace, level + 1);
		pair_multiply(v, p, 2 * size);
		p_one = mul(p_one, saved);
		if (rows == 1) {
			fft<true>(v, 2 * columns, columns);
			fft<true>(v + 2 * columns, 2 * columns, columns);
			U factor = in(T(out(mul(p_one, in(4 * size)))).inv().val());
			composition_finish(p, v, columns, factor);
		} else {
			U *upper = v + 2 * size;
			double_y(upper, rows, columns, true);
			composition_forward_y(upper, size, columns);
			p_one = add(p_one, p_one);
			composition_difference(p, v, upper, 2 * size);
			project(p, rows, columns);
			p_one = mul(p_one, in(2 * columns));
		}
	}
};
inline NTT &ntt_plan(int size) {
	static NTT t(1, false);
	if (static_cast<int>(t.roots.size()) < size) t = NTT(size, false);
	return t;
}

}

inline void ntt(std::vector<fps_mint> &a, bool inverse = false) {
	assert(!a.empty() && !(a.size() & (a.size() - 1)) && a.size() <= (1U << 23));
	int n = static_cast<int>(a.size());
	auto &t = fps_detail::ntt_plan(std::max(1, n / 2));
	auto reorder = [&] {
		for (int i = 1, j = 0; i < n; ++i) {
			int bit = n >> 1;
			for (; j & bit; bit >>= 1) j ^= bit;
			j ^= bit;
			if (i < j) std::swap(a[i], a[j]);
		}
	};
	static_assert(sizeof(fps_mint) == sizeof(fps_detail::NTT::U) && std::is_trivially_copyable_v<fps_mint>);
	auto *words = reinterpret_cast<fps_detail::NTT::U *>(a.data());
	if (inverse) {
		reorder();
		t.fft<true>(words, n);
		t.scale(words, n, fps_mint(n).inv().val());
	} else {
		t.fft<false>(words, n);
		reorder();
	}
}

#if defined(__GNUC__) && !defined(__clang__) && (defined(__x86_64__) || defined(__i386__))
#pragma GCC pop_options
#endif

#endif
