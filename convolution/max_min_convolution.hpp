#ifndef MAX_MIN_CONVOLUTION_HPP
#define MAX_MIN_CONVOLUTION_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <numeric>
#include <vector>

#include "../math/fps.hpp"

namespace max_min_convolution_detail {
template <class T, class Compare>
std::vector<T> multiply(const std::vector<T> &a, const std::vector<T> &b, Compare compare) {
	if (a.empty() || b.empty()) return {};
	assert(a.size() + b.size() - 1 <= std::size_t(1) << 23);
	int n = a.size(), m = b.size(), length = n + m - 1;
	std::vector<T> result(length, a[0]);
	std::vector<std::uint8_t> done(length);
	if (std::min(n, m) <= 32) {
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j) {
				T value = compare(a[i], b[j]) ? b[j] : a[i];
				if (!done[i + j] || compare(value, result[i + j])) result[i + j] = value;
				done[i + j] = 1;
			}
		return result;
	}
	std::vector<int> order(n + m), pending;
	std::iota(order.begin(), order.end(), 0);
	auto value = [&](int i) -> decltype(auto) { return i < n ? a[i] : b[i - n]; };
	std::sort(order.begin(), order.end(), [&](int i, int j) { return compare(value(i), value(j)); });
	int block = std::max(1, int(std::sqrt(double(n + m) * std::log2(n + m))));
	fps active_a(n), active_b(m);
	int remaining = length;
	for (int first = 0; first < n + m && remaining; first += block) {
		int last = std::min(n + m, first + block);
		fps x = active_a, y = active_b;
		for (int k = first; k < last; ++k) {
			int i = order[k];
			if (i < n) x[i] = 1;
			else y[i - n] = 1;
		}
		auto counts = x * y;
		pending.clear();
		for (int k = 0; k < length; ++k)
			if (!done[k] && counts[k] != mint(0)) pending.push_back(k);
		for (int k = first; k < last; ++k) {
			int i = order[k];
			for (int pos : pending) {
				if (done[pos]) continue;
				int j = pos - (i < n ? i : i - n);
				bool exists = i < n ? 0 <= j && j < m && active_b[j] != mint(0) : 0 <= j && j < n && active_a[j] != mint(0);
				if (exists) result[pos] = value(i), done[pos] = 1, --remaining;
			}
			if (i < n) active_a[i] = 1;
			else active_b[i - n] = 1;
		}
	}
	return result;
}
}

template <class T>
std::vector<T> max_min_convolution(const std::vector<T> &a, const std::vector<T> &b) {
	return max_min_convolution_detail::multiply(a, b, std::greater<T>{});
}
template <class T>
std::vector<T> min_max_convolution(const std::vector<T> &a, const std::vector<T> &b) {
	return max_min_convolution_detail::multiply(a, b, std::less<T>{});
}

#endif
