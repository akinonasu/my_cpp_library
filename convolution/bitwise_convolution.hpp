#ifndef BITWISE_CONVOLUTION_HPP
#define BITWISE_CONVOLUTION_HPP

#include "../components/bitwise_zeta.hpp"
#include "../components/walsh_hadamard.hpp"
#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace bitwise_convolution_detail {
template <int kind, class T>
std::vector<T> multiply(std::vector<T> a, std::vector<T> b) {
	if (a.empty() || b.empty()) return {};
	std::size_t n = std::bit_ceil(std::max(a.size(), b.size()));
	a.resize(n, T(0));
	b.resize(n, T(0));
	if constexpr (kind == 2) {
		walsh_hadamard(a);
		walsh_hadamard(b);
	} else {
		bitwise_zeta<false, kind == 0>(a);
		bitwise_zeta<false, kind == 0>(b);
	}
	for (std::size_t i = 0; i < n; ++i) a[i] *= b[i];
	if constexpr (kind == 2) {
		walsh_hadamard(a);
		if constexpr (std::numeric_limits<T>::is_integer)
			for (T &x : a) x /= std::int64_t(n);
		else {
			T inverse = T(1) / T(n);
			for (T &x : a) x *= inverse;
		}
	} else bitwise_zeta<true, kind == 0>(a);
	return a;
}
}

template <class T>
std::vector<T> and_convolution(std::vector<T> a, std::vector<T> b) {
	return bitwise_convolution_detail::multiply<0>(std::move(a), std::move(b));
}
template <class T>
std::vector<T> or_convolution(std::vector<T> a, std::vector<T> b) {
	return bitwise_convolution_detail::multiply<1>(std::move(a), std::move(b));
}
template <class T>
std::vector<T> xor_convolution(std::vector<T> a, std::vector<T> b) {
	return bitwise_convolution_detail::multiply<2>(std::move(a), std::move(b));
}

#endif
