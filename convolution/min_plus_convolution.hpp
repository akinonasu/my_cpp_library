#ifndef MIN_PLUS_CONVOLUTION_HPP
#define MIN_PLUS_CONVOLUTION_HPP

#include "../components/smawk.hpp"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <vector>

namespace min_plus_convolution_detail {

template <class T, class Compare>
std::vector<T> multiply(const std::vector<T> &a, const std::vector<T> &b, Compare compare) {
	if (a.empty() || b.empty()) return {};
	assert(a.size() + b.size() <= std::size_t(std::numeric_limits<int>::max()));
	int n = a.size(), m = b.size(), length = n + m - 1;
	std::vector<int> columns(m);
	std::iota(columns.begin(), columns.end(), 0);
	auto better = [&](int row, int x, int y) {
		if (row < y) return false;
		if (row - x >= n) return true;
		return compare(a[row - y] + b[y], a[row - x] + b[x]);
	};
	auto best = smawk(length, columns, better);
	std::vector<T> result;
	result.reserve(length);
	for (int i = 0; i < length; ++i) result.push_back(a[i - best[i]] + b[best[i]]);
	return result;
}
template <class T, class Compare>
std::vector<T> opposite(const std::vector<T> &a, const std::vector<T> &b, Compare compare) {
	if (a.empty() || b.empty()) return {};
	assert(a.size() + b.size() <= std::size_t(std::numeric_limits<int>::max()));
	int n = a.size(), m = b.size(), length = n + m - 1;
	std::vector<T> result(length, a[0] + b[0]);
	std::vector<std::uint8_t> done(length);
	auto solve = [&](auto &&self, int l, int r, int parent_l, int parent_r) -> void {
		int low = std::max(0, r - n), high = std::min(m, l + 1);
		std::vector<int> columns;
		for (int i = high - 1; i >= low;) {
			if (parent_l <= i && i < parent_r) i = parent_l - 1;
			else columns.push_back(i--);
		}
		if (!columns.empty()) {
			auto best = smawk(r - l, columns, [&](int row, int x, int y) {
				return compare(a[l + row - y] + b[y], a[l + row - x] + b[x]);
			});
			for (int i = l; i < r; ++i) {
				T value = a[i - best[i - l]] + b[best[i - l]];
				if (!done[i] || compare(value, result[i])) result[i] = value;
				done[i] = 1;
			}
		}
		if (r - l > 1) {
			int mid = l + (r - l) / 2;
			self(self, l, mid, low, high);
			self(self, mid, r, low, high);
		}
	};
	solve(solve, 0, length, 0, 0);
	return result;
}

}

template <class T>
std::vector<T> min_plus_convolution_convex(const std::vector<T> &a, const std::vector<T> &b) {
	return min_plus_convolution_detail::multiply(a, b, std::less<T>{});
}
template <class T>
std::vector<T> max_plus_convolution_concave(const std::vector<T> &a, const std::vector<T> &b) {
	return min_plus_convolution_detail::multiply(a, b, std::greater<T>{});
}

template <class T>
std::vector<T> min_plus_convolution_concave(const std::vector<T> &a, const std::vector<T> &b) {
	return min_plus_convolution_detail::opposite(a, b, std::less<T>{});
}
template <class T>
std::vector<T> max_plus_convolution_convex(const std::vector<T> &a, const std::vector<T> &b) {
	return min_plus_convolution_detail::opposite(a, b, std::greater<T>{});
}

#endif
