#ifndef BOOLEAN_CONVOLUTION_HPP
#define BOOLEAN_CONVOLUTION_HPP

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include "bitwise_convolution.hpp"
#include "gcd_lcm_convolution.hpp"
#include "subset_convolution.hpp"

namespace boolean_convolution_detail {
template <int index, bool dual = false, class T>
std::vector<T> multiply(const std::vector<T> &a, const std::vector<T> &b) {
	static_assert(std::is_integral_v<T> && sizeof(T) <= 8 && 0 <= index && index <= 5);
	if (a.empty() || b.empty()) return {};
	using U = std::make_unsigned_t<std::conditional_t<std::is_same_v<T, bool>, unsigned char, T>>;
	std::size_t n = std::max(a.size(), b.size());
	if constexpr (index < 3 || index == 5) n = std::bit_ceil(n);
	assert(n <= std::size_t(1) << 30);
	constexpr int width = std::is_same_v<T, bool> ? 1 : std::numeric_limits<U>::digits;
	std::vector<T> result(n, dual ? T(~U(0)) : T(0));
	for (int bit = 0; bit < width; ++bit) {
		std::vector<__int128_t> x(n), y(n);
		for (std::size_t i = index == 3 || index == 4; i < a.size(); ++i) x[i] = ((U(a[i]) >> bit) & 1) ^ dual;
		for (std::size_t i = index == 3 || index == 4; i < b.size(); ++i) y[i] = ((U(b[i]) >> bit) & 1) ^ dual;
		std::vector<__int128_t> counts;
		if constexpr (index < 3) counts = bitwise_convolution_detail::multiply<index>(std::move(x), std::move(y));
		else if constexpr (index == 5) counts = subset_convolution(std::move(x), std::move(y));
		else counts = gcd_lcm_convolution_detail::multiply<index == 3>(std::move(x), std::move(y));
		for (std::size_t i = 0; i < n; ++i)
			if (counts[i]) {
				U flag = U(1) << bit;
				result[i] = T(dual ? U(result[i]) & ~flag : U(result[i]) | flag);
			}
	}
	if constexpr (index == 3 || index == 4) result[0] = T(0);
	return result;
}
}

template <class T>
std::vector<T> and_or_convolution(const std::vector<T> &a, const std::vector<T> &b) { return boolean_convolution_detail::multiply<0>(a, b); }
template <class T>
std::vector<T> or_or_convolution(const std::vector<T> &a, const std::vector<T> &b) { return boolean_convolution_detail::multiply<1>(a, b); }
template <class T>
std::vector<T> xor_or_convolution(const std::vector<T> &a, const std::vector<T> &b) { return boolean_convolution_detail::multiply<2>(a, b); }
template <class T>
std::vector<T> gcd_or_convolution(const std::vector<T> &a, const std::vector<T> &b) { return boolean_convolution_detail::multiply<3>(a, b); }
template <class T>
std::vector<T> lcm_or_convolution(const std::vector<T> &a, const std::vector<T> &b) { return boolean_convolution_detail::multiply<4>(a, b); }
template <class T>
std::vector<T> and_and_convolution(const std::vector<T> &a, const std::vector<T> &b) { return boolean_convolution_detail::multiply<0, true>(a, b); }
template <class T>
std::vector<T> or_and_convolution(const std::vector<T> &a, const std::vector<T> &b) { return boolean_convolution_detail::multiply<1, true>(a, b); }
template <class T>
std::vector<T> xor_and_convolution(const std::vector<T> &a, const std::vector<T> &b) { return boolean_convolution_detail::multiply<2, true>(a, b); }
template <class T>
std::vector<T> gcd_and_convolution(const std::vector<T> &a, const std::vector<T> &b) { return boolean_convolution_detail::multiply<3, true>(a, b); }
template <class T>
std::vector<T> lcm_and_convolution(const std::vector<T> &a, const std::vector<T> &b) { return boolean_convolution_detail::multiply<4, true>(a, b); }

template <class T>
std::vector<T> subset_or_convolution(const std::vector<T> &a, const std::vector<T> &b) { return boolean_convolution_detail::multiply<5>(a, b); }
template <class T>
std::vector<T> subset_and_convolution(const std::vector<T> &a, const std::vector<T> &b) { return boolean_convolution_detail::multiply<5, true>(a, b); }

#endif
