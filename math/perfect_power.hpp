#pragma once

#include <bit>
#include <cassert>
#include <cmath>
#include <utility>

namespace perfect_power_impl {

inline constexpr int PRIMES[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61};

// base^exp と target をオーバーフローなしで比較する (<: -1, ==: 0, >: 1)
inline int cmp_pow(unsigned long base, int exp, unsigned long target) {
    if (base == 0) return target == 0 ? 0 : -1;
    if (base == 1) return target == 1 ? 0 : (1 < target ? -1 : 1);
    unsigned long res = 1;
    unsigned long cur = base;
    while (exp > 0) {
        if (exp & 1) {
            if (cur > target / res) return 1;
            res *= cur;
        }
        exp >>= 1;
        if (exp > 0) {
            if (cur > target / cur) return 1;
            cur *= cur;
        }
    }
    if (res < target) return -1;
    if (res > target) return 1;
    return 0;
}

// p 乗数になり得るかの軽量な剰余フィルタ
inline bool pass_mod_filter(int p, unsigned long n) {
    if (p == 5) {
        unsigned long m = n % 11;
        return m <= 1 || m == 10;
    }
    if (p == 7) {
        unsigned long m = n % 29;
        return m <= 1 || m == 12 || m == 17 || m == 28;
    }
    if (p == 11) {
        unsigned long m = n % 23;
        return m <= 1 || m == 22;
    }
    if (p == 23) {
        unsigned long m = n % 47;
        return m <= 1 || m == 46;
    }
    if (p == 29) {
        unsigned long m = n % 59;
        return m <= 1 || m == 58;
    }
    return true;
}

// floor(n^(1/p)) を厳密に求める (p >= 1)
inline unsigned long integer_nth_root_u64(unsigned long n, int p) {
    if (p == 1 || n <= 1) return n;
    if (p >= 64) return 1;

    if (p >= 10) {
        int bw = std::bit_width(n);
        if (bw <= p) return 1;
        unsigned long lo = 1UL << ((bw - 1) / p);
        unsigned long hi = 1UL << ((bw + p - 1) / p);
        while (lo + 1 < hi) {
            unsigned long mid = lo + (hi - lo) / 2;
            if (cmp_pow(mid, p, n) <= 0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return lo;
    }

    double approx;
    if (p == 2) {
        approx = std::sqrt(static_cast<double>(n));
    } else if (p == 3) {
        approx = std::cbrt(static_cast<double>(n));
    } else {
        approx = std::pow(static_cast<double>(n), 1.0 / p);
    }

    unsigned long r = static_cast<unsigned long>(approx);
    while (cmp_pow(r + 1, p, n) <= 0)
        ++r;
    while (cmp_pow(r, p, n) > 0)
        --r;
    return r;
}

}  // namespace perfect_power_impl

// 非負整数 n に対して floor(n^(1/p)) を返す (n >= 0, p >= 1)
inline long integer_nth_root(long n, int p) {
    assert(n >= 0);
    assert(p >= 1);
    return static_cast<long>(perfect_power_impl::integer_nth_root_u64(static_cast<unsigned long>(n), p));
}

// n = a^b (b >= 2) を満たす整数 a, b が存在するか判定する (負数・0・1にも対応)
inline bool is_perfect_power(long n) {
    if (n >= -1 && n <= 1) return true;

    unsigned long abs_n = (n < 0) ? (0UL - static_cast<unsigned long>(n)) : static_cast<unsigned long>(n);

    int tz = std::countr_zero(abs_n);
    if (tz > 0) {
        if (tz == 1) return false;
        int rem = tz;
        if ((rem & 1) == 0) {
            if (n > 0) {
                unsigned long r = perfect_power_impl::integer_nth_root_u64(abs_n, 2);
                if (r * r == abs_n) return true;
            }
            while ((rem & 1) == 0)
                rem >>= 1;
        }
        for (int p = 3; p * p <= rem; p += 2) {
            if (rem % p == 0) {
                if (perfect_power_impl::pass_mod_filter(p, abs_n)) {
                    unsigned long r = perfect_power_impl::integer_nth_root_u64(abs_n, p);
                    if (perfect_power_impl::cmp_pow(r, p, abs_n) == 0) return true;
                }
                while (rem % p == 0)
                    rem /= p;
            }
        }
        if (rem > 1 && perfect_power_impl::pass_mod_filter(rem, abs_n)) {
            unsigned long r = perfect_power_impl::integer_nth_root_u64(abs_n, rem);
            if (perfect_power_impl::cmp_pow(r, rem, abs_n) == 0) return true;
        }
        return false;
    }

    if (n > 0 && (abs_n & 7) == 1) {
        unsigned long r = perfect_power_impl::integer_nth_root_u64(abs_n, 2);
        if (r * r == abs_n) return true;
    }

    for (int i = 1; perfect_power_impl::PRIMES[i] <= 37; ++i) {
        int p = perfect_power_impl::PRIMES[i];
        if (perfect_power_impl::cmp_pow(3, p, abs_n) > 0) break;
        if (!perfect_power_impl::pass_mod_filter(p, abs_n)) continue;
        unsigned long r = perfect_power_impl::integer_nth_root_u64(abs_n, p);
        if (perfect_power_impl::cmp_pow(r, p, abs_n) == 0) return true;
    }
    return false;
}

// n = a^b (b >= 2) に分解する。
// 累乗数なら {true, {a, b}} (b は最大)、そうでなければ {false, {n, 1}} を返す。
// 0, 1 は {true, {n, 2}}, -1 は {true, {-1, 3}} を返す。
inline std::pair<bool, std::pair<long, int>> decompose_perfect_power(long n) {
    if (n == 0) return {true, {0L, 2}};
    if (n == 1) return {true, {1L, 2}};
    if (n == -1) return {true, {-1L, 3}};

    unsigned long abs_n = (n < 0) ? (0UL - static_cast<unsigned long>(n)) : static_cast<unsigned long>(n);
    int b = 1;
    int tz = std::countr_zero(abs_n);

    auto try_factor = [&](int p) {
        while (perfect_power_impl::pass_mod_filter(p, abs_n)) {
            unsigned long r = perfect_power_impl::integer_nth_root_u64(abs_n, p);
            if (perfect_power_impl::cmp_pow(r, p, abs_n) != 0) break;
            abs_n = r;
            b *= p;
        }
    };

    if (tz > 0) {
        if (tz == 1) return {false, {n, 1}};
        int rem = tz;
        if ((rem & 1) == 0) {
            if (n > 0) try_factor(2);
            while ((rem & 1) == 0)
                rem >>= 1;
        }
        for (int p = 3; p * p <= rem; p += 2) {
            if (rem % p == 0) {
                try_factor(p);
                while (rem % p == 0)
                    rem /= p;
            }
        }
        if (rem > 1) {
            try_factor(rem);
        }
    } else {
        if (n > 0 && (abs_n & 7) == 1) {
            try_factor(2);
        }
        for (int i = 1; perfect_power_impl::PRIMES[i] <= 37; ++i) {
            int p = perfect_power_impl::PRIMES[i];
            if (perfect_power_impl::cmp_pow(3, p, abs_n) > 0) break;
            try_factor(p);
        }
    }

    if (b == 1) return {false, {n, 1}};
    long a = (n < 0) ? static_cast<long>(0UL - abs_n) : static_cast<long>(abs_n);
    return {true, {a, b}};
}
