#pragma once

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "components/default_init_allocator.hpp"

template <class S, auto op, auto e = nullptr>
class SparseTable {
    using Value = std::conditional_t<std::is_same_v<S, bool>, std::uint8_t, S>;
    using Reference = std::conditional_t<std::is_same_v<S, bool>, bool, const S&>;
    int n = 0;
    std::vector<std::vector<Value, DefaultInitAllocator<Value>>> data;

   public:
    SparseTable() = default;
    SparseTable(const SparseTable&) = default;
    SparseTable& operator=(const SparseTable&) = default;
    SparseTable(SparseTable&& other) noexcept : n(std::exchange(other.n, 0)), data(std::move(other.data)) {}
    SparseTable& operator=(SparseTable&& other) noexcept {
        if (this != &other) {
            n = std::exchange(other.n, 0);
            data = std::move(other.data);
        }
        return *this;
    }
    explicit SparseTable(const std::vector<S>& a) { build(a); }
    void clear() {
        data.clear();
        n = 0;
    }
    void build(const std::vector<S>& a) {
        assert(a.size() <= std::size_t(std::numeric_limits<int>::max()));
        build(static_cast<int>(a.size()), [&](int i) -> S { return a[i]; });
    }
    template <class F>
    void build(int size, F value) {
        assert(size >= 0);
        clear();
        n = size;
        if (!n) return;
        data.reserve(std::bit_width(unsigned(n)));
        data.emplace_back();
        if constexpr (std::is_default_constructible_v<Value>) {
            data[0].resize(n);
            for (int i = 0; i < n; ++i)
                data[0][i] = value(i);
        } else {
            data[0].reserve(n);
            for (int i = 0; i < n; ++i)
                data[0].push_back(value(i));
        }
        for (int k = 1; k < std::bit_width(unsigned(n)); ++k) {
            int length = n - (1U << k) + 1, half = 1U << (k - 1);
            data.emplace_back();
            if constexpr (std::is_default_constructible_v<Value>) {
                data[k].resize(length);
                for (int i = 0; i < length; ++i)
                    data[k][i] = op(S(data[k - 1][i]), S(data[k - 1][i + half]));
            } else {
                data[k].reserve(length);
                for (int i = 0; i < length; ++i)
                    data[k].push_back(op(data[k - 1][i], data[k - 1][i + half]));
            }
        }
    }
    int size() const { return n; }
    bool empty() const { return n == 0; }
    Reference operator[](int i) const {
        assert(0 <= i && i < n);
        return data[0][i];
    }
    [[gnu::always_inline]] S prod(int l, int r) const {
        assert(0 <= l && l <= r && r <= n);
        if constexpr (std::is_same_v<decltype(e), std::nullptr_t>)
            assert(l < r);
        else if (l == r)
            return e();
        return prod_unchecked(l, r);
    }
    [[gnu::always_inline]] S prod_unchecked(int l, int r) const {
        int k = std::bit_width(unsigned(r - l)) - 1;
        return op(S(data[k][l]), S(data[k][r - (1U << k)]));
    }
    S all_prod() const { return prod(0, n); }
};