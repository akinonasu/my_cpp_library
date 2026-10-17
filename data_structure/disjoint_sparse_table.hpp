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
class DisjointSparseTable {
    using Value = std::conditional_t<std::is_same_v<S, bool>, std::uint8_t, S>;
    using Reference = std::conditional_t<std::is_same_v<S, bool>, bool, const S&>;
    int n = 0;
    std::vector<std::vector<Value, DefaultInitAllocator<Value>>> data;
    static decltype(auto) combine(const Value& a, const Value& b) {
        if constexpr (std::is_same_v<S, bool>)
            return bool(op(bool(a), bool(b)));
        else
            return op(a, b);
    }

   public:
    DisjointSparseTable() = default;
    DisjointSparseTable(const DisjointSparseTable&) = default;
    DisjointSparseTable& operator=(const DisjointSparseTable&) = default;
    DisjointSparseTable(DisjointSparseTable&& other) noexcept : n(std::exchange(other.n, 0)), data(std::move(other.data)) {}
    DisjointSparseTable& operator=(DisjointSparseTable&& other) noexcept {
        if (this != &other) {
            n = std::exchange(other.n, 0);
            data = std::move(other.data);
        }
        return *this;
    }
    explicit DisjointSparseTable(const std::vector<S>& a) { build(a); }
    void build(const std::vector<S>& a) {
        assert(a.size() <= std::size_t(std::numeric_limits<int>::max()));
        n = static_cast<int>(a.size());
        int levels = std::max(1, int(std::bit_width(unsigned(std::max(0, n - 1)))));
        data.clear();
        data.resize(levels);
        for (int k = 0; k < levels; ++k) {
            std::size_t length = n, half = std::size_t(1) << k;
            if (k && (length & (half * 2 - 1)) <= half) length &= ~(half * 2 - 1);
            if constexpr (std::is_default_constructible_v<Value>)
                data[k].resize(length);
            else
                data[k].assign(a.begin(), a.begin() + length);
        }
        for (int i = 0; i < n; ++i)
            data[0][i] = a[i];
        for (int k = 1; k < levels; ++k) {
            Value* row = data[k].data();
            std::size_t half = std::size_t(1) << k, block = half * 2;
            for (std::size_t start = 0; start < data[k].size(); start += block) {
                std::size_t mid = start + half, end = std::min(start + block, data[k].size());
                row[mid - 1] = a[mid - 1];
                row[mid] = a[mid];
                for (std::size_t i = mid - 1; i > start; --i)
                    row[i - 1] = combine(a[i - 1], row[i]);
                for (std::size_t i = mid + 1; i < end; ++i)
                    row[i] = combine(row[i - 1], a[i]);
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
        if (l == --r) return data[0][l];
        const auto& row = data[std::bit_width(unsigned(l ^ r)) - 1];
        return combine(row[l], row[r]);
    }
    S all_prod() const { return prod(0, n); }
};
