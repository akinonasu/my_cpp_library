#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <type_traits>
#include <utility>

#include "components/d_ary_heap.hpp"

template <class T = std::int64_t, class Sum = T>
class SlopeTrick {
    static_assert((std::is_arithmetic_v<T> && std::is_signed_v<T>) || std::is_same_v<T, __int128_t>);
    static_assert((std::is_arithmetic_v<Sum> && std::is_signed_v<Sum>) || std::is_same_v<Sum, __int128_t>);
    DaryHeap<T> left;
    DaryHeap<T, std::greater<T>> right;
    T offset_left = 0, offset_right = 0;
    Sum minimum = 0;

   public:
    SlopeTrick() = default;
    SlopeTrick(const SlopeTrick&) = default;
    SlopeTrick& operator=(const SlopeTrick&) = default;
    SlopeTrick(SlopeTrick&& other) noexcept
        : left(std::move(other.left)),
          right(std::move(other.right)),
          offset_left(std::exchange(other.offset_left, 0)),
          offset_right(std::exchange(other.offset_right, 0)),
          minimum(std::exchange(other.minimum, Sum(0))) {}
    SlopeTrick& operator=(SlopeTrick&& other) noexcept {
        if (this != &other) {
            left = std::move(other.left);
            right = std::move(other.right);
            offset_left = std::exchange(other.offset_left, 0);
            offset_right = std::exchange(other.offset_right, 0);
            minimum = std::exchange(other.minimum, Sum(0));
        }
        return *this;
    }
    std::size_t size() const { return left.size() + right.size(); }
    bool empty() const { return size() == 0; }
    void reserve(std::size_t n) {
        left.reserve(n);
        right.reserve(n);
    }
    void clear_left() {
        left.clear();
        offset_left = 0;
    }
    void clear_right() {
        right.clear();
        offset_right = 0;
    }
    void clear() {
        clear_left();
        clear_right();
        minimum = 0;
    }
    Sum min() const { return minimum; }
    std::pair<std::optional<T>, std::optional<T>> argmin() const {
        return {left.empty() ? std::nullopt : std::optional<T>(left.top() + offset_left),
                right.empty() ? std::nullopt : std::optional<T>(right.top() + offset_right)};
    }
    void add_all(Sum a) { minimum += a; }
    void add_a_minus_x(T a) {
        if (!right.empty() && right.top() + offset_right < a) {
            T r = right.replace_top(a - offset_right) + offset_right;
            minimum += Sum(a) - Sum(r);
            left.push(r - offset_left);
        } else
            left.push(a - offset_left);
    }
    void add_x_minus_a(T a) {
        if (!left.empty() && a < left.top() + offset_left) {
            T l = left.replace_top(a - offset_left) + offset_left;
            minimum += Sum(l) - Sum(a);
            right.push(l - offset_right);
        } else
            right.push(a - offset_right);
    }
    void add_abs(T a) {
        if (!left.empty() && a < left.top() + offset_left) {
            T l = left.replace_top(a - offset_left) + offset_left;
            minimum += Sum(l) - Sum(a);
            left.push(a - offset_left);
            right.push(l - offset_right);
        } else if (!right.empty() && right.top() + offset_right < a) {
            T r = right.replace_top(a - offset_right) + offset_right;
            minimum += Sum(a) - Sum(r);
            right.push(a - offset_right);
            left.push(r - offset_left);
        } else
            left.push(a - offset_left), right.push(a - offset_right);
    }
    void shift(T a, T b) {
        assert(a <= b);
        offset_left += a;
        offset_right += b;
    }
    void shift(T a) { shift(a, a); }
    Sum eval(T x) const {
        Sum value = minimum;
        for (T a : left.values())
            if (x < a + offset_left) value += Sum(a + offset_left) - Sum(x);
        for (T a : right.values())
            if (a + offset_right < x) value += Sum(x) - Sum(a + offset_right);
        return value;
    }
    void swap(SlopeTrick& other) noexcept {
        left.swap(other.left);
        right.swap(other.right);
        std::swap(offset_left, other.offset_left);
        std::swap(offset_right, other.offset_right);
        std::swap(minimum, other.minimum);
    }
    void merge(SlopeTrick& other) {
        if (this == &other) {
            SlopeTrick copy(other);
            merge(copy);
            return;
        }
        if (size() < other.size()) swap(other);
        minimum += other.minimum;
        for (T a : other.left.values())
            add_a_minus_x(a + other.offset_left);
        for (T a : other.right.values())
            add_x_minus_a(a + other.offset_right);
        other = SlopeTrick{};
    }
    void merge(SlopeTrick&& other) { merge(other); }
};
