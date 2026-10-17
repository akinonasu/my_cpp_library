#ifndef ORDERED_CONVOLUTION_HPP
#define ORDERED_CONVOLUTION_HPP

#include <algorithm>
#include <cstddef>
#include <vector>

#include "convolution.hpp"

template <class T>
std::vector<T> ordered_convolution(const std::vector<T> &a, const std::vector<T> &b) {
	if (a.empty() || b.empty()) return {};
	std::vector<T> result(a.size() + b.size() - 1, T(0));
	auto solve = [&](auto &&self, std::size_t l, std::size_t r) -> void {
		if (r - l <= 1 || l >= a.size() || l >= b.size()) return;
		std::size_t mid = (l + r) / 2;
		if (mid < b.size()) {
			std::vector<T> x(a.begin() + l, a.begin() + std::min(mid, a.size()));
			std::vector<T> y(b.begin() + mid, b.begin() + std::min(r, b.size()));
			auto product = convolution(x, y);
			for (std::size_t i = 0; i < product.size(); ++i) result[l + mid + i] += product[i];
		}
		self(self, l, mid);
		self(self, mid, r);
	};
	solve(solve, 0, std::max(a.size(), b.size()));
	return result;
}

#endif
