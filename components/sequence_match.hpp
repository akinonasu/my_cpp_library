#ifndef SEQUENCE_MATCH_HPP
#define SEQUENCE_MATCH_HPP

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

namespace sequence_match_detail {
#if defined(__x86_64__)
__attribute__((target("avx2"))) inline std::size_t forward_avx2(const unsigned char *a, const unsigned char *b, std::size_t size) {
	std::size_t i = 0;
	for (; size - i >= 32; i += 32) {
		auto x = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(a + i));
		auto y = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(b + i));
		auto difference = ~std::uint32_t(_mm256_movemask_epi8(_mm256_cmpeq_epi8(x, y)));
		if (difference) return i + std::countr_zero(difference);
	}
	while (i < size && a[i] == b[i]) ++i;
	return i;
}

__attribute__((target("avx2"))) inline std::size_t backward_avx2(const unsigned char *a, const unsigned char *b, std::size_t size) {
	std::size_t i = 0;
	for (; size - i >= 32; i += 32) {
		auto x = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(a - i - 32));
		auto y = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(b - i - 32));
		auto difference = ~std::uint32_t(_mm256_movemask_epi8(_mm256_cmpeq_epi8(x, y)));
		if (difference) return i + std::countl_zero(difference);
	}
	while (i < size && a[-std::ptrdiff_t(i) - 1] == b[-std::ptrdiff_t(i) - 1]) ++i;
	return i;
}
#endif

template <class T, bool Vectorize = true>
class Matcher {
	const T *data;
	std::vector<std::array<int, 2>> runs;
	static constexpr bool bytes = std::is_integral_v<T> && std::has_unique_object_representations_v<T>;
public:
	explicit Matcher(std::span<const T> values) : data(values.data()) {}
	void index_runs(int n) {
		runs.resize(n);
		for (int first = 0, last; first < n; first = last) {
			last = first + 1;
			while (last < n && data[first] == data[last]) ++last;
			for (int i = first; i < last; ++i) runs[i] = {first, last};
		}
	}
	std::pair<int, int> forward_counted(int i, int j, int limit) const {
		if (runs.empty() || !limit || data[i] != data[j]) {
			int length = forward(i, j, limit);
			return {length, length};
		}
		int a = runs[i][1] - i, b = runs[j][1] - j, skip = std::min({limit, a, b});
		if (a != b || skip == limit) return {skip, 1};
		int length = forward(i + skip, j + skip, limit - skip);
		return {skip + length, length + 1};
	}
	std::pair<int, int> backward_counted(int i, int j, int limit) const {
		if (runs.empty() || !limit || data[i - 1] != data[j - 1]) {
			int length = backward(i, j, limit);
			return {length, length};
		}
		int a = i - runs[i - 1][0], b = j - runs[j - 1][0], skip = std::min({limit, a, b});
		if (a != b || skip == limit) return {skip, 1};
		int length = backward(i - skip, j - skip, limit - skip);
		return {skip + length, length + 1};
	}
	int forward(int i, int j, int limit) const {
		if (!limit) return 0;
		int matched = 0;
		if constexpr (bytes) {
			const auto *a = reinterpret_cast<const unsigned char *>(data + i);
			const auto *b = reinterpret_cast<const unsigned char *>(data + j);
			std::size_t offset = 0, size = std::size_t(limit) * sizeof(T);
			while (size - offset >= 8) {
				std::uint64_t x, y;
				std::memcpy(&x, a + offset, 8);
				std::memcpy(&y, b + offset, 8);
				auto difference = x ^ y;
				if (difference) {
					int bits = std::endian::native == std::endian::little ? std::countr_zero(difference) : std::countl_zero(difference);
					return static_cast<int>((offset + bits / 8) / sizeof(T));
				}
				offset += 8;
#if defined(__x86_64__)
				if constexpr (Vectorize) if (size - offset >= 64 && __builtin_cpu_supports("avx2"))
					return static_cast<int>((offset + forward_avx2(a + offset, b + offset, size - offset)) / sizeof(T));
#endif
			}
			matched = static_cast<int>(offset / sizeof(T));
		}
		while (matched < limit && data[i + matched] == data[j + matched]) ++matched;
		return matched;
	}
	int backward(int i, int j, int limit) const {
		if (!limit) return 0;
		int matched = 0;
		if constexpr (bytes) {
			const auto *a = reinterpret_cast<const unsigned char *>(data + i);
			const auto *b = reinterpret_cast<const unsigned char *>(data + j);
			std::size_t offset = 0, size = std::size_t(limit) * sizeof(T);
			while (size - offset >= 8) {
				std::uint64_t x, y;
				std::memcpy(&x, a - offset - 8, 8);
				std::memcpy(&y, b - offset - 8, 8);
				auto difference = x ^ y;
				if (difference) {
					int bits = std::endian::native == std::endian::little ? std::countl_zero(difference) : std::countr_zero(difference);
					return static_cast<int>((offset + bits / 8) / sizeof(T));
				}
				offset += 8;
#if defined(__x86_64__)
				if constexpr (Vectorize) if (size - offset >= 64 && __builtin_cpu_supports("avx2"))
					return static_cast<int>((offset + backward_avx2(a - offset, b - offset, size - offset)) / sizeof(T));
#endif
			}
			matched = static_cast<int>(offset / sizeof(T));
		}
		while (matched < limit && data[i - matched - 1] == data[j - matched - 1]) ++matched;
		return matched;
	}
};

template <class F>
void z_table(int n, int *z, F match) {
	if (!n) return;
	z[0] = n;
	for (int i = 1, base = 0, end = 0; i < n; ++i) {
		int length = 0;
		if (i < end) {
			length = z[i - base];
			if (length < end - i) {
				z[i] = length;
				continue;
			}
			length = end - i;
		}
		length += match(length, i + length, n - i - length);
		z[i] = length;
		if (i + length > end) base = i, end = i + length;
	}
}

template <class F, class Emit>
void prefix_matches(int n, int m, const int *z, F match, Emit emit) {
	for (int i = 0, base = 0, end = 0; i < n; ++i) {
		int length = 0;
		if (i < end) {
			length = z[i - base];
			if (length < end - i) {
				emit(i, length);
				continue;
			}
			length = end - i;
		}
		length += match(length, i + length, m - length);
		if (i + length > end) base = i, end = i + length;
		emit(i, length);
	}
}
}

#endif
