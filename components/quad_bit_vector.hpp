#ifndef QUAD_BIT_VECTOR_HPP
#define QUAD_BIT_VECTOR_HPP

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>
#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

struct QuadBitVector {
	struct Block {
		std::uint64_t low = 0, high = 0;
		std::array<int, 4> prefix{};
	};
	int n;
	std::vector<Block> blocks;

#if defined(__x86_64__) || defined(__i386__)
	__attribute__((target("avx2"))) static void pack_avx2(const std::uint32_t *a, int n, int shift, Block *blocks) {
		auto low_shift = _mm_cvtsi32_si128(31 - shift), high_shift = _mm_cvtsi32_si128(30 - shift);
		for (int w = 0; w < (n >> 6) + bool(n & 63); ++w) {
			int start = w * 64;
			std::uint64_t low = 0, high = 0;
			int length = std::min(64, n - start), j = 0;
			for (; j + 8 <= length; j += 8) {
				auto x = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(a + start + j));
				low |= std::uint64_t(_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_sll_epi32(x, low_shift)))) << j;
				high |= std::uint64_t(_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_sll_epi32(x, high_shift)))) << j;
			}
			for (; j < length; ++j) {
				low |= std::uint64_t((a[start + j] >> shift) & 1) << j;
				high |= std::uint64_t((a[start + j] >> (shift + 1)) & 1) << j;
			}
			blocks[w].low = low;
			blocks[w].high = high;
		}
		blocks[n >> 6].low = n & 63 ? blocks[n >> 6].low : 0;
		blocks[n >> 6].high = n & 63 ? blocks[n >> 6].high : 0;
	}
#endif
	explicit QuadBitVector(int n = 0) : n(n), blocks((n >> 6) + 1) {}
	void build(const std::vector<std::uint32_t> &a, int shift) {
		assert(a.size() == std::size_t(n) && 0 <= shift && shift <= 30);
		[[maybe_unused]] bool vectorized = false;
#if defined(__x86_64__) || defined(__i386__)
		vectorized = n >= 64 && __builtin_cpu_supports("avx2");
		if (vectorized) pack_avx2(a.data(), n, shift, blocks.data());
#endif
		std::array<int, 4> counts{};
		for (int w = 0; w < static_cast<int>(blocks.size()); ++w) {
			std::uint64_t low = 0, high = 0;
			int start = w * 64, length = std::min(64, n - start);
			if (vectorized) low = blocks[w].low, high = blocks[w].high;
			else
				for (int j = 0; j < length; ++j) {
					low |= std::uint64_t((a[start + j] >> shift) & 1) << j;
					high |= std::uint64_t((a[start + j] >> (shift + 1)) & 1) << j;
				}
			blocks[w] = {low, high, counts};
			int c1 = std::popcount(low & ~high), c2 = std::popcount(~low & high), c3 = std::popcount(low & high);
			counts[0] += length - c1 - c2 - c3;
			counts[1] += c1;
			counts[2] += c2;
			counts[3] += c3;
		}
	}
	std::uint64_t mask(const Block &b, int digit) const {
		return ((digit & 1) ? b.low : ~b.low) & ((digit & 2) ? b.high : ~b.high);
	}

	int rank(int r, int digit) const {
		const auto &b = blocks[r >> 6];
		return b.prefix[digit] + std::popcount(mask(b, digit) & ((std::uint64_t(1) << (r & 63)) - 1));
	}
	int below(int r, int digit) const {
		if (!digit) return 0;
		if (digit == 1) return rank(r, 0);
		if (digit == 3) return r - rank(r, 3);
		const auto &b = blocks[r >> 6];
		return b.prefix[0] + b.prefix[1] + std::popcount(~b.high & ((std::uint64_t(1) << (r & 63)) - 1));
	}
	std::array<int, 4> ranks(int r) const {
		const auto &b = blocks[r >> 6];
		std::uint64_t m = (std::uint64_t(1) << (r & 63)) - 1;
		std::array<int, 4> a;
		a[1] = b.prefix[1] + std::popcount(b.low & ~b.high & m);
		a[2] = b.prefix[2] + std::popcount(~b.low & b.high & m);
		a[3] = b.prefix[3] + std::popcount(b.low & b.high & m);
		a[0] = r - a[1] - a[2] - a[3];
		return a;
	}
};

#endif
