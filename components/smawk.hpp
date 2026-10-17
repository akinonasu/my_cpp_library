#ifndef SMAWK_HPP
#define SMAWK_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <numeric>
#include <vector>

template <class Better>
std::vector<int> smawk(int length, const std::vector<int> &columns, Better better) {
	assert(length >= 0 && (!length || !columns.empty()));
	if (!length) return {};
	std::vector<int> best(length);
	auto solve = [&](auto &&self, std::size_t first, std::size_t stride, const std::vector<int> &cols) -> void {
		int rows = (length - 1 - first) / stride + 1;
		std::vector<int> reduced;
		reduced.reserve(std::min<std::size_t>(rows, cols.size()));
		for (int col : cols) {
			while (!reduced.empty() && better(first + (int(reduced.size()) - 1) * stride, reduced.back(), col)) reduced.pop_back();
			if (int(reduced.size()) < rows) reduced.push_back(col);
		}
		if (rows > 1) self(self, first + stride, stride * 2, reduced);
		std::size_t left = 0;
		for (std::size_t row = first; row < std::size_t(length); row += 2 * stride) {
			std::size_t right = left;
			if (row + stride < std::size_t(length))
				while (reduced[right] != best[row + stride]) ++right;
			else right = reduced.size() - 1;
			best[row] = reduced[left];
			for (std::size_t i = left + 1; i <= right; ++i)
				if (better(row, best[row], reduced[i])) best[row] = reduced[i];
			left = right;
		}
	};
	solve(solve, 0, 1, columns);
	return best;
}

template <class Better>
std::vector<int> smawk(int rows, int columns, Better better) {
	assert(columns >= 0);
	std::vector<int> cols(columns);
	std::iota(cols.begin(), cols.end(), 0);
	return smawk(rows, cols, better);
}

#endif
