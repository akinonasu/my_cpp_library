#ifndef CYCLIC_CONVOLUTION_HPP
#define CYCLIC_CONVOLUTION_HPP

#include <cassert>
#include <cstddef>
#include <vector>

#include "convolution.hpp"

template <bool negative = false, class T>
std::vector<T> cyclic_convolution(const std::vector<T> &a, const std::vector<T> &b, std::size_t n) {
	assert(n > 0);
	std::vector<T> x(n, T(0)), y(n, T(0));
	for (std::size_t i = 0; i < a.size(); ++i) {
		if constexpr (negative) {
			if ((i / n) & 1) x[i % n] -= a[i];
			else x[i % n] += a[i];
		} else x[i % n] += a[i];
	}
	for (std::size_t i = 0; i < b.size(); ++i) {
		if constexpr (negative) {
			if ((i / n) & 1) y[i % n] -= b[i];
			else y[i % n] += b[i];
		} else y[i % n] += b[i];
	}
	auto product = convolution(x, y);
	x.assign(n, T(0));
	for (std::size_t i = 0; i < product.size(); ++i) {
		if constexpr (negative) {
			if (i >= n) x[i - n] -= product[i];
			else x[i] += product[i];
		} else x[i % n] += product[i];
	}
	return x;
}

#endif
