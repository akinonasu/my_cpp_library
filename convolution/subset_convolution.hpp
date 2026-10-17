#ifndef SUBSET_CONVOLUTION_HPP
#define SUBSET_CONVOLUTION_HPP

#include <algorithm>
#include <bit>
#include <cstddef>
#include <vector>

template <class T>
std::vector<T> subset_convolution(const std::vector<T> &a, const std::vector<T> &b) {
	if (a.empty() || b.empty()) return {};
	std::size_t n = std::bit_ceil(std::max(a.size(), b.size()));
	int levels = std::bit_width(n) - 1, width = levels + 1;
	std::vector<T> x(n * width, T(0)), y(n * width, T(0)), row(width, T(0)), result(n, T(0));
	for (std::size_t i = 0; i < a.size(); ++i) x[i * width + std::popcount(i)] = a[i];
	for (std::size_t i = 0; i < b.size(); ++i) y[i * width + std::popcount(i)] = b[i];
	for (std::size_t half = 1; half < n; half *= 2)
		for (std::size_t start = 0; start < n; start += 2 * half)
			for (std::size_t j = start; j < start + half; ++j)
				for (int k = 0, end = std::popcount(j); k <= end; ++k) {
					x[(j + half) * width + k] += x[j * width + k];
					y[(j + half) * width + k] += y[j * width + k];
				}
	for (std::size_t mask = 0; mask < n; ++mask) {
		int count = std::popcount(mask);
		std::fill(row.begin(), row.end(), T(0));
		for (int i = 0; i <= count; ++i)
			for (int j = std::max(0, count - i); j <= count && i + j <= levels; ++j)
				row[i + j] += x[mask * width + i] * y[mask * width + j];
		std::copy(row.begin(), row.end(), x.begin() + mask * width);
	}
	for (std::size_t half = 1; half < n; half *= 2)
		for (std::size_t start = 0; start < n; start += 2 * half)
			for (std::size_t j = start; j < start + half; ++j)
				for (int k = std::popcount(j + half); k <= levels; ++k)
					x[(j + half) * width + k] -= x[j * width + k];
	for (std::size_t mask = 0; mask < n; ++mask) result[mask] = x[mask * width + std::popcount(mask)];
	return result;
}

#endif
