#ifndef SMALL_VECTOR_HPP
#define SMALL_VECTOR_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <memory>
#include <type_traits>
#include <utility>

template <class T, std::size_t N = 2>
class SmallVector {
	static_assert(N && std::is_trivially_copyable_v<T> && std::is_trivially_destructible_v<T>);
	static constexpr std::size_t flag = std::size_t(1) << (std::numeric_limits<std::size_t>::digits - 1);
	struct Heap { T *data; std::size_t capacity; };
	union Storage {
		T local[N];
		Heap heap;
		Storage() : local{} {}
		~Storage() {}
	} storage;
	std::size_t length = 0;
	bool large() const { return length & flag; }
	void set_size(std::size_t n) { length = (length & flag) | n; }

public:
	SmallVector() = default;
	explicit SmallVector(std::size_t n) { resize(n); }
	SmallVector(std::initializer_list<T> values) { assign(values.begin(), values.end()); }
	template <std::input_iterator It>
	SmallVector(It first, It last) { assign(first, last); }
	SmallVector(const SmallVector &b) { assign(b.begin(), b.end()); }
	SmallVector(SmallVector &&b) noexcept : length(b.length) {
		if (large()) storage.heap = b.storage.heap;
		else std::memcpy(storage.local, b.storage.local, size() * sizeof(T));
		b.length = 0;
		b.storage.local[0] = T{};
	}
	SmallVector &operator=(const SmallVector &b) {
		if (this != &b) assign(b.begin(), b.end());
		return *this;
	}
	SmallVector &operator=(SmallVector &&b) noexcept {
		if (this == &b) return *this;
		if (large()) delete[] storage.heap.data;
		length = b.length;
		if (large()) storage.heap = b.storage.heap;
		else std::memcpy(storage.local, b.storage.local, size() * sizeof(T));
		b.length = 0;
		b.storage.local[0] = T{};
		return *this;
	}
	~SmallVector() { if (large()) delete[] storage.heap.data; }
	std::size_t size() const { return length & ~flag; }
	std::size_t capacity() const { return large() ? storage.heap.capacity : N; }
	bool empty() const { return size() == 0; }
	T *data() { return large() ? storage.heap.data : storage.local; }
	const T *data() const { return large() ? storage.heap.data : storage.local; }
	T *begin() { return data(); }
	const T *begin() const { return data(); }
	T *end() { return data() + size(); }
	const T *end() const { return data() + size(); }
	T &operator[](std::size_t i) { return data()[i]; }
	const T &operator[](std::size_t i) const { return data()[i]; }
	T &back() { return data()[size() - 1]; }
	const T &back() const { return data()[size() - 1]; }
	void reserve(std::size_t n) {
		if (n <= capacity()) return;
		assert(n < flag / sizeof(T));
		std::size_t capacity = std::max(n, this->capacity() * 2);
		auto memory = std::make_unique_for_overwrite<T[]>(capacity);
		std::memcpy(memory.get(), data(), size() * sizeof(T));
		if (large()) delete[] storage.heap.data;
		storage.heap = {memory.release(), capacity};
		length |= flag;
	}
	void resize(std::size_t n) {
		std::size_t old = size();
		reserve(n);
		if (n > old) std::fill(data() + old, data() + n, T{});
		set_size(n);
	}
	void resize_for_overwrite(std::size_t n) { reserve(n); set_size(n); }
	void clear() { set_size(0); }
	void push_back(T value) {
		std::size_t n = size();
		if (n == capacity()) reserve(n + 1);
		data()[n] = value;
		set_size(n + 1);
	}
	void pop_back() { assert(!empty()); set_size(size() - 1); }
	void assign(std::size_t n, T value) {
		reserve(n);
		std::fill_n(n > N ? storage.heap.data : data(), n, value);
		set_size(n);
	}
	template <std::input_iterator It>
	void assign(It first, It last) {
		if constexpr (std::forward_iterator<It>) {
			std::size_t n = std::distance(first, last);
			reserve(n);
			std::copy(first, last, data());
			set_size(n);
		} else {
			clear();
			for (; first != last; ++first) push_back(*first);
		}
	}
	T *insert(const T *position, std::size_t count, T value) {
		std::size_t offset = position - data(), old = size();
		reserve(old + count);
		T *p = data();
		std::memmove(p + offset + count, p + offset, (old - offset) * sizeof(T));
		std::fill(p + offset, p + offset + count, value);
		set_size(old + count);
		return p + offset;
	}
	template <std::forward_iterator It>
	T *insert(const T *position, It first, It last) {
		std::size_t offset = position - data(), old = size(), count = std::distance(first, last);
		reserve(old + count);
		T *p = data();
		std::memmove(p + offset + count, p + offset, (old - offset) * sizeof(T));
		std::copy(first, last, p + offset);
		set_size(old + count);
		return p + offset;
	}
	T *erase(const T *first, const T *last) {
		std::size_t offset = first - data(), count = last - first;
		std::memmove(data() + offset, last, (end() - last) * sizeof(T));
		set_size(size() - count);
		return data() + offset;
	}
	friend bool operator==(const SmallVector &a, const SmallVector &b) {
		return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
	}
};

#endif
