#ifndef CONVOLUTION_RING_HPP
#define CONVOLUTION_RING_HPP

#include <algorithm>
#include <cstddef>
#include <span>
#include <vector>

namespace convolution_ring_detail {
template <class T>
std::vector<T> multiply(std::span<const T> a, std::span<const T> b) {
	if (a.empty() || b.empty()) return {};
	std::vector<T> result(a.size() + b.size() - 1, T(0));
	if (std::min(a.size(), b.size()) <= 32) {
		for (std::size_t i = 0; i < a.size(); ++i)
			for (std::size_t j = 0; j < b.size(); ++j) result[i + j] += a[i] * b[j];
		return result;
	}
	std::size_t half = std::max(a.size(), b.size()) / 2;
	auto al = a.first(std::min(a.size(), half)), ar = a.subspan(al.size());
	auto bl = b.first(std::min(b.size(), half)), br = b.subspan(bl.size());
	auto low = multiply(al, bl), high = multiply(ar, br);
	std::vector<T> x(std::max(al.size(), ar.size()), T(0)), y(std::max(bl.size(), br.size()), T(0));
	for (std::size_t i = 0; i < al.size(); ++i) x[i] += al[i];
	for (std::size_t i = 0; i < ar.size(); ++i) x[i] += ar[i];
	for (std::size_t i = 0; i < bl.size(); ++i) y[i] += bl[i];
	for (std::size_t i = 0; i < br.size(); ++i) y[i] += br[i];
	auto middle = multiply<T>(x, y);
	for (std::size_t i = 0; i < low.size(); ++i) result[i] += low[i], middle[i] -= low[i];
	for (std::size_t i = 0; i < high.size(); ++i) result[i + 2 * half] += high[i], middle[i] -= high[i];
	for (std::size_t i = 0; i < middle.size() && i + half < result.size(); ++i) result[i + half] += middle[i];
	return result;
}
}

template <class T>
std::vector<T> convolution_ring(const std::vector<T> &a, const std::vector<T> &b) {
	return convolution_ring_detail::multiply<T>(a, b);
}

#endif
