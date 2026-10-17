#ifndef WALSH_HADAMARD_HPP
#define WALSH_HADAMARD_HPP

#include <bit>
#include <cassert>
#include <cstddef>
#include <vector>

template <class T>
void walsh_hadamard(std::vector<T> &a) {
	assert(a.empty() || std::has_single_bit(a.size()));
	for (std::size_t half = 1; half < a.size(); half *= 2)
		for (std::size_t start = 0; start < a.size(); start += 2 * half)
			for (std::size_t j = 0; j < half; ++j) {
				T x = a[start + j], y = a[start + j + half];
				a[start + j] = x + y;
				a[start + j + half] = x - y;
			}
}

#endif
