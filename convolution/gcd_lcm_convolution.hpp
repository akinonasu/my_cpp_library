#ifndef GCD_LCM_CONVOLUTION_HPP
#define GCD_LCM_CONVOLUTION_HPP

#include "../components/divisor_zeta.hpp"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

namespace gcd_lcm_convolution_detail {

template <bool gcd, class T>
std::vector<T> multiply(std::vector<T> a, std::vector<T> b) {
	if (a.empty() || b.empty()) return {};
	assert(a[0] == T(0) && b[0] == T(0));
	std::size_t n = std::max(a.size(), b.size());
	a.resize(n, T(0));
	b.resize(n, T(0));
	auto ps = prime_sieve(n);
	divisor_zeta<false, gcd>(a, ps);
	divisor_zeta<false, gcd>(b, ps);
	for (std::size_t i = 1; i < n; ++i) a[i] *= b[i];
	divisor_zeta<true, gcd>(a, ps);
	return a;
}
}

template <class T>
std::vector<T> gcd_convolution(std::vector<T> a, std::vector<T> b) {
	return gcd_lcm_convolution_detail::multiply<true>(std::move(a), std::move(b));
}
template <class T>
std::vector<T> lcm_convolution(std::vector<T> a, std::vector<T> b) {
	return gcd_lcm_convolution_detail::multiply<false>(std::move(a), std::move(b));
}

#endif
