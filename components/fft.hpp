#ifndef FFT_HPP
#define FFT_HPP

#include <algorithm>
#include <bit>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <vector>

template <class Real>
void fft(std::vector<std::complex<Real>> &a, bool inverse = false) {
	std::size_t n = a.size();
	assert(!n || std::has_single_bit(n));
	for (std::size_t i = 1, j = 0; i < n; ++i) {
		std::size_t bit = n >> 1;
		for (; j & bit; bit >>= 1) j ^= bit;
		j ^= bit;
		if (i < j) std::swap(a[i], a[j]);
	}
	for (std::size_t width = 2; width <= n; width *= 2) {
		Real angle = (inverse ? -2 : 2) * std::numbers::pi_v<Real> / Real(width);
		std::complex<Real> step(std::cos(angle), std::sin(angle));
		for (std::size_t start = 0; start < n; start += width) {
			std::complex<Real> w(1);
			for (std::size_t j = 0; j < width / 2; ++j, w *= step) {
				auto x = a[start + j], y = a[start + j + width / 2] * w;
				a[start + j] = x + y;
				a[start + j + width / 2] = x - y;
			}
		}
	}
	if (inverse)
		for (auto &x : a) x /= Real(n);
}

#endif
