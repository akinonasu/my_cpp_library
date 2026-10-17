#ifndef D_ARY_HEAP_HPP
#define D_ARY_HEAP_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

template <class T, class Compare = std::less<T>, std::size_t D = 2>
class DaryHeap {
	static_assert(D >= 2);
	std::vector<T> data;
	[[no_unique_address]] Compare compare{};
	void up(std::size_t i, T value) {
		while (i) {
			std::size_t p = (i - 1) / D;
			if (!compare(data[p], value)) break;
			data[i] = std::move(data[p]);
			i = p;
		}
		data[i] = std::move(value);
	}
	void down(std::size_t i, T value) {
		std::size_t n = data.size();
		while (i < n / D + (n % D > 1)) {
			std::size_t best = i * D + 1, end = n;
			if constexpr (D <= std::size_t(-1) / 2) end = std::min(best + D, n);
			for (std::size_t j = best + 1; j < end; ++j)
				if (compare(data[best], data[j])) best = j;
			if (!compare(value, data[best])) break;
			data[i] = std::move(data[best]);
			i = best;
		}
		data[i] = std::move(value);
	}

  public:
	DaryHeap() = default;
	explicit DaryHeap(Compare compare) : compare(std::move(compare)) {}
	explicit DaryHeap(std::vector<T> values, Compare compare = {}) : data(std::move(values)), compare(std::move(compare)) {
		for (std::size_t i = data.size() / D + (data.size() % D > 1); i--;) down(i, std::move(data[i]));
	}
	bool empty() const { return data.empty(); }
	std::size_t size() const { return data.size(); }
	void reserve(std::size_t n) { data.reserve(n); }
	void clear() { data.clear(); }
	const std::vector<T> &values() const { return data; }
	decltype(auto) top() const {
		assert(!empty());
		return data.front();
	}
	void push(T value) {
		data.push_back(std::move(value));
		up(data.size() - 1, std::move(data.back()));
	}
	T replace_top(T value) {
		assert(!empty());
		T result = std::move(data.front());
		down(0, std::move(value));
		return result;
	}
	T push_pop(T value) {
		if (empty() || !compare(value, data.front())) return value;
		return replace_top(std::move(value));
	}
	T pop_value() {
		assert(!empty());
		T result = std::move(data.front());
		if (data.size() == 1) data.pop_back();
		else {
			T value = std::move(data.back());
			data.pop_back();
			down(0, std::move(value));
		}
		return result;
	}
	void pop() { pop_value(); }
	void swap(DaryHeap &other) noexcept(noexcept(std::swap(compare, other.compare))) {
		data.swap(other.data);
		std::swap(compare, other.compare);
	}
};

#endif
