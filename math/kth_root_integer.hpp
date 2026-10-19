#pragma once
#include <bit>
#include <cassert>
#include <cstdint>

inline std::uint64_t kth_root_integer(std::uint64_t n, int k = 2) {
    assert(k >= 1);
    if (k == 1 || n <= 1) return n;
    if (k >= 64) return 1;
    auto fits = [n, k](std::uint64_t x) {
        std::uint64_t value = 1;
        for (int i = 0; i < k; ++i) {
            if (value > n / x) return false;
            value *= x;
        }
        return true;
    };
    std::uint64_t lo = 1, hi = std::uint64_t{1} << ((std::bit_width(n) + k - 1) / k);
    while (hi - lo > 1) {
        auto mid = lo + (hi - lo) / 2;
        if (fits(mid))
            lo = mid;
        else
            hi = mid;
    }
    return lo;
}