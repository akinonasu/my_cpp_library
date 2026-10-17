#ifndef CONVOLUTION_FFT_HPP
#define CONVOLUTION_FFT_HPP

#include "../components/fft.hpp"
#include <algorithm>
#include <bit>
#include <complex>
#include <cstddef>
#include <type_traits>
#include <vector>

template <class Real>
std::vector<std::complex<Real>> convolution_fft(const std::vector<std::complex<Real>> &a, const std::vector<std::complex<Real>> &b) {
	if (a.empty() || b.empty()) return {};
	std::size_t length = a.size() + b.size() - 1, n = std::bit_ceil(length);
	std::vector<std::complex<Real>> x(a.begin(), a.end()), y(b.begin(), b.end());
	x.resize(n);
	y.resize(n);
	fft(x, false);
	if (&a == &b) y = x;
	else fft(y, false);
	for (std::size_t i = 0; i < n; ++i) x[i] *= y[i];
	fft(x, true);
	x.resize(length);
	return x;
}
template <class Real>
	requires std::is_floating_point_v<Real>
std::vector<Real> convolution_fft(const std::vector<Real> &a, const std::vector<Real> &b) {
	std::vector<std::complex<Real>> x(a.begin(), a.end()), y(b.begin(), b.end());
	auto values = &a == &b ? convolution_fft(x, x) : convolution_fft(x, y);
	std::vector<Real> result;
	result.reserve(values.size());
	for (const auto &x : values) result.push_back(x.real());
	return result;
}

#endif
