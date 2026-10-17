#ifndef MULTIVARIATE_CONVOLUTION_HPP
#define MULTIVARIATE_CONVOLUTION_HPP

#include <bit>
#include <cassert>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

#include "convolution.hpp"
#include "../components/ntt.hpp"

template <class T>
std::vector<T> multivariate_convolution(const std::vector<T> &a, const std::vector<T> &b, const std::vector<int> &dimensions) {
	std::size_t size = 1;
	std::vector<std::size_t> bases;
	for (int d : dimensions) {
		assert(d > 0 && size <= std::size_t(-1) / d);
		if (d == 1) continue;
		size *= d;
		bases.push_back(size);
	}
	assert(a.size() == size && b.size() == size);
	if (bases.empty()) return {a[0] * b[0]};
	int colors = bases.size();
	std::vector<int> rank(size);
	for (std::size_t i = 0; i < size; ++i) {
		std::size_t value = 0;
		for (std::size_t base : bases) value += i / base;
		rank[i] = value % colors;
	}
	if constexpr (requires { T::mod(); T::raw(0); a[0].val(); }) {
		if constexpr (T::mod() == 998244353) {
			assert(size <= std::size_t(1) << 22);
			std::size_t length = std::bit_ceil(2 * size - 1);
			using NTT = fps_detail::NTT;
			NTT transform(length / 2, false);
			std::vector<std::vector<std::uint32_t>> x(colors, std::vector<std::uint32_t>(length)), y = x;
			for (std::size_t i = 0; i < size; ++i) x[rank[i]][i] = NTT::in(a[i].val()), y[rank[i]][i] = NTT::in(b[i].val());
			for (int i = 0; i < colors; ++i) {
				transform.fft<false>(x[i].data(), length);
				transform.fft<false>(y[i].data(), length);
			}
			std::vector<std::uint32_t> row(colors);
			for (std::size_t k = 0; k < length; ++k) {
				for (int c = 0; c < colors; ++c) {
					std::uint32_t value = 0;
					for (int i = 0; i < colors; ++i) {
						value = NTT::add(value, NTT::mul(x[i][k], y[(c - i + colors) % colors][k]));
					}
					row[c] = value;
				}
				for (int c = 0; c < colors; ++c) x[c][k] = row[c];
			}
			for (auto &v : x) {
				transform.fft<true>(v.data(), length);
				transform.scale(v.data(), length, 998244353 - 998244352 / length);
			}
			std::vector<T> result;
			result.reserve(size);
			for (std::size_t i = 0; i < size; ++i) result.push_back(T::raw(NTT::out(x[rank[i]][i])));
			return result;
		}
	}
	std::size_t stride = 2 * colors - 1;
	assert(size <= std::size_t(-1) / stride / 2);
	std::vector<T> x(size * stride, T(0)), y(size * stride, T(0));
	for (std::size_t i = 0; i < size; ++i) x[i * stride + rank[i]] = a[i], y[i * stride + rank[i]] = b[i];
	auto product = convolution(x, y);
	std::vector<T> result(size, T(0));
	for (std::size_t i = 0; i < size; ++i) {
		result[i] = product[i * stride + rank[i]];
		if (rank[i] + colors < int(stride)) result[i] += product[i * stride + rank[i] + colors];
	}
	return result;
}

#endif
