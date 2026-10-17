#ifndef DIVISOR_ZETA_HPP
#define DIVISOR_ZETA_HPP

#include "../math/prime_sieve.hpp"
#include <cstddef>
#include <vector>

template <bool inverse = false, bool multiples = false, class T>
void divisor_zeta(std::vector<T> &a, const std::vector<std::size_t> &ps) {
	if (a.empty()) return;
	for (std::size_t p : ps) {
		std::size_t end = (a.size() - 1) / p;
		for (std::size_t step = 1; step <= end; ++step) {
			std::size_t i = inverse == multiples ? step : end + 1 - step;
			T &to = a[multiples ? i : i * p];
			const T &from = a[multiples ? i * p : i];
			if constexpr (inverse) to -= from;
			else to += from;
		}
	}
}

template <bool inverse = false, bool multiples = false, class T>
void divisor_zeta(std::vector<T> &a) { divisor_zeta<inverse, multiples>(a, prime_sieve(a.size())); }

#endif
