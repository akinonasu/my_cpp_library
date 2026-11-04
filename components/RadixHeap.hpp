#ifndef RADIX_HEAP_HPP
#define RADIX_HEAP_HPP

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

template <class Key = std::uint64_t, class T = int>
class RadixHeap {
	static_assert(std::is_integral_v<Key> && !std::is_same_v<Key, bool> && sizeof(Key) <= 8);
	using U = std::make_unsigned_t<Key>;
	using State = std::pair<Key, T>;
	std::vector<State> bucket[std::numeric_limits<U>::digits + 1];
	Key lower;
	std::size_t count = 0;
	int index(Key key) const { return std::bit_width(U(U(key) ^ U(lower))); }
	void refill() {
		assert(count);
		if (!bucket[0].empty()) return;
		int i = 1;
		while (bucket[i].empty()) ++i;
		lower = bucket[i][0].first;
		for (const auto &s : bucket[i]) lower = std::min(lower, s.first);
		for (auto &s : bucket[i]) bucket[index(s.first)].push_back(std::move(s));
		bucket[i].clear();
	}

  public:
	explicit RadixHeap(Key minimum = 0) : lower(minimum) {}
	RadixHeap(const RadixHeap &) = default;
	RadixHeap &operator=(const RadixHeap &) = default;
	RadixHeap(RadixHeap &&other) noexcept : lower(other.lower), count(std::exchange(other.count, 0)) {
		for (int i = 0; i <= std::numeric_limits<U>::digits; ++i) bucket[i] = std::move(other.bucket[i]);
	}
	RadixHeap &operator=(RadixHeap &&other) noexcept {
		if (this != &other) {
			lower = other.lower;
			count = std::exchange(other.count, 0);
			for (int i = 0; i <= std::numeric_limits<U>::digits; ++i) bucket[i] = std::move(other.bucket[i]);
		}
		return *this;
	}

	bool empty() const { return count == 0; }
	std::size_t size() const { return count; }
	Key last() const { return lower; }
	void clear(Key minimum = 0) {
		for (auto &b : bucket) b.clear();
		lower = minimum;
		count = 0;
	}
	void push(State s) {
		assert(s.first >= lower);
		bucket[index(s.first)].push_back(std::move(s));
		++count;
	}
	const State &top() {
		refill();
		return bucket[0].back();
	}
	void pop() {
		refill();
		bucket[0].pop_back();
		--count;
	}
	State pop_value() {
		refill();
		State s = std::move(bucket[0].back());
		bucket[0].pop_back();
		--count;
		return s;
	}
};

#endif
