#include <algorithm>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <numeric>
#include <print>
#include <ranges>
#include <utility>
#include <vector>

// ============================================================================
// 前計算あり・高速クエリ用（目安: N <= 10^7, int型）
// ============================================================================
class PrimeSieve {
   public:
    explicit PrimeSieve(int n) : n_(n) {
        assert(0 <= n && n < std::numeric_limits<int>::max());
        spf_.assign(n_ + 1, 0);
        for (int i = 2; i <= n_; ++i) {
            if (spf_[i] == 0) {
                spf_[i] = i;
                primes_.push_back(i);
            }
            for (int p : primes_) {
                if (p > spf_[i] || static_cast<long>(i) * p > n_) break;
                spf_[i * p] = p;
            }
        }
    }

    // 0 <= x <= n_ の素数判定: O(1)
    bool is_prime(int x) const {
        assert(0 <= x && x <= n_);
        return x >= 2 && spf_[x] == x;
    }

    // 0 <= x <= n_ の最小素因数（x < 2 の場合は -1）: O(1)
    int min_factor(int x) const {
        assert(0 <= x && x <= n_);
        if (x < 2) return -1;
        return spf_[x];
    }

    // n_ 以下の全素数（昇順）
    const std::vector<int>& primes() const { return primes_; }

    // 半開区間 [l, r) の全素数（昇順）: O(log(pi(N)) + 出力数)
    std::vector<int> primes(int l, int r) const {
        assert(0 <= l && l <= r && r <= n_ + 1);
        auto it_l = std::ranges::lower_bound(primes_, l);
        auto it_r = std::ranges::lower_bound(primes_, r);
        return std::vector<int>(it_l, it_r);
    }

    // 0 <= x <= n_ 以下の素数の個数: O(log(pi(N)))
    int count_primes(int x) const {
        assert(0 <= x && x <= n_);
        return static_cast<int>(std::ranges::upper_bound(primes_, x) - primes_.begin());
    }

    // 半開区間 [l, r) の素数の個数: O(log(pi(N)))
    int count_primes(int l, int r) const {
        assert(0 <= l && l <= r && r <= n_ + 1);
        return static_cast<int>(std::ranges::lower_bound(primes_, r) - std::ranges::lower_bound(primes_, l));
    }

    // 0-indexed で k 番目の素数（n_ 以下に存在しない場合は -1）: O(1)
    int kth_prime(int k) const {
        assert(k >= 0);
        if (k >= std::ssize(primes_)) return -1;
        return primes_[k];
    }

    // 1 <= x <= n_ の素因数分解 {素因数, 指数}（昇順）: O(log x)
    std::vector<std::pair<int, int>> factorize(int x) const {
        assert(1 <= x && x <= n_);
        std::vector<std::pair<int, int>> res;
        while (x > 1) {
            int p = spf_[x];
            int cnt = 0;
            while (spf_[x] == p) {
                x /= p;
                ++cnt;
            }
            res.emplace_back(p, cnt);
        }
        return res;
    }

    // 1 <= x <= n_ の全約数（昇順）: O(log x + d(x) log d(x))
    std::vector<int> divisors(int x) const {
        assert(1 <= x && x <= n_);
        std::vector<int> res = {1};
        for (auto [p, e] : factorize(x)) {
            int sz = static_cast<int>(res.size());
            res.reserve(sz * (e + 1));
            int mul = 1;
            for (int i = 0; i < e; ++i) {
                mul *= p;
                for (int j = 0; j < sz; ++j) {
                    res.push_back(res[j] * mul);
                }
            }
        }
        std::ranges::sort(res);
        return res;
    }

    // 1 <= x <= n_ のオイラーのトーシェント関数 phi(x): O(log x)
    int totient(int x) const {
        assert(1 <= x && x <= n_);
        int res = x;
        while (x > 1) {
            int p = spf_[x];
            res -= res / p;
            while (spf_[x] == p) {
                x /= p;
            }
        }
        return res;
    }

   private:
    int n_;
    std::vector<int> spf_;
    std::vector<int> primes_;
};

// ============================================================================
// 単発・64bit整数向け高速アルゴリズム（long型, 最大 9.22 * 10^18 まで対応）
// ============================================================================
namespace prime {

namespace detail {
__extension__ using u128 = unsigned __int128;

inline long mul_mod(long a, long b, long mod) {
    return static_cast<long>((static_cast<u128>(a) * b) % mod);
}

inline long isqrt(long n) {
    assert(n >= 0);
    if (n == 0) return 0;
    long x = static_cast<long>(std::sqrt(static_cast<double>(n)));
    while (x + 1 <= n / (x + 1))
        ++x;
    while (x > n / x)
        --x;
    return x;
}
}  // namespace detail

// (base^exp) % mod: O(log exp)
inline long pow_mod(long base, long exp, long mod) {
    assert(exp >= 0 && mod > 0);
    if (mod == 1) return 0;
    base %= mod;
    if (base < 0) base += mod;
    long res = 1;
    while (exp > 0) {
        if (exp & 1L) res = detail::mul_mod(res, base, mod);
        base = detail::mul_mod(base, base, mod);
        exp >>= 1;
    }
    return res;
}

// 決定的 Miller-Rabin 素数判定（64bit符号付き非負整数全域で正確）: O(log^3 n)
inline bool is_prime(long n) {
    if (n < 2) return false;
    for (long p : {2L, 3L, 5L, 7L, 11L, 13L, 17L, 19L, 23L, 29L, 31L, 37L}) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }
    if (n < 41L * 41L) return true;

    long d = n - 1;
    int s = std::countr_zero(static_cast<unsigned long>(d));
    d >>= s;

    auto check = [&](long a) -> bool {
        if (a % n == 0) return true;
        long x = pow_mod(a, d, n);
        if (x == 1 || x == n - 1) return true;
        for (int r = 1; r < s; ++r) {
            x = detail::mul_mod(x, x, n);
            if (x == n - 1) return true;
            if (x == 1) return false;
        }
        return false;
    };

    if (n <= 4294967295L) {
        for (long a : {2L, 7L, 61L}) {
            if (!check(a)) return false;
        }
        return true;
    }
    for (long a : {2L, 325L, 9375L, 28178L, 450775L, 9780504L, 1795265022L}) {
        if (!check(a)) return false;
    }
    return true;
}

namespace detail {
inline long pollard_rho(long n) {
    if (n % 2 == 0) return 2;
    unsigned long state = static_cast<unsigned long>(n) ^ 0x9e3779b97f4a7c15UL;
    auto next_rand = [&]() -> unsigned long {
        state += 0x9e3779b97f4a7c15UL;
        unsigned long z = state;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9UL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebUL;
        return z ^ (z >> 31);
    };

    constexpr long m = 128;
    while (true) {
        long c = static_cast<long>(next_rand() % static_cast<unsigned long>(n - 1)) + 1;
        long y = static_cast<long>(next_rand() % static_cast<unsigned long>(n));
        auto f = [&](long x) -> long { return static_cast<long>((static_cast<u128>(x) * x + c) % n); };
        long g = 1, r = 1, q = 1;
        long x = y, ys = y;
        while (g == 1) {
            x = y;
            long k = 0;
            while (k < r && g == 1) {
                ys = y;
                long lim = std::min(m, r - k);
                for (long i = 0; i < lim; ++i) {
                    y = f(y);
                    long diff = x > y ? x - y : y - x;
                    q = mul_mod(q, diff, n);
                }
                g = std::gcd(q, n);
                k += m;
            }
            r <<= 1;
        }
        if (g == n) {
            do {
                ys = f(ys);
                long diff = x > ys ? x - ys : ys - x;
                g = std::gcd(diff, n);
            } while (g == 1);
        }
        if (g != n) return g;
    }
}
}  // namespace detail

// n >= 1 の素因数分解 {素因数, 指数}（昇順）: 期待 O(n^(1/4) log n)
inline std::vector<std::pair<long, int>> factorize(long n) {
    assert(n >= 1);
    std::vector<std::pair<long, int>> res;
    for (long p : {2L, 3L, 5L, 7L, 11L, 13L, 17L, 19L, 23L, 29L, 31L, 37L}) {
        if (n % p == 0) {
            int cnt = 0;
            while (n % p == 0) {
                n /= p;
                ++cnt;
            }
            res.emplace_back(p, cnt);
        }
    }
    if (n == 1) return res;

    std::vector<long> stk = {n};
    std::vector<long> factors;
    while (!stk.empty()) {
        long x = stk.back();
        stk.pop_back();
        if (x == 1) continue;
        if (is_prime(x)) {
            factors.push_back(x);
            continue;
        }
        long d = detail::pollard_rho(x);
        stk.push_back(d);
        stk.push_back(x / d);
    }

    std::ranges::sort(factors);
    for (long p : factors) {
        if (res.empty() || res.back().first != p) {
            res.emplace_back(p, 1);
        } else {
            ++res.back().second;
        }
    }
    return res;
}

// n >= 1 の全約数（昇順）
inline std::vector<long> divisors(long n) {
    assert(n >= 1);
    std::vector<long> res = {1};
    for (auto [p, e] : factorize(n)) {
        int sz = static_cast<int>(res.size());
        res.reserve(sz * (e + 1));
        long mul = 1;
        for (int i = 0; i < e; ++i) {
            mul *= p;
            for (int j = 0; j < sz; ++j) {
                res.push_back(res[j] * mul);
            }
        }
    }
    std::ranges::sort(res);
    return res;
}

// n >= 1 のオイラーのトーシェント関数 phi(n)
inline long totient(long n) {
    assert(n >= 1);
    long res = n;
    for (auto [p, e] : factorize(n)) {
        res -= res / p;
    }
    return res;
}

// n >= 1 の最小の原始根（存在しない場合は -1）
inline long primitive_root(long n) {
    assert(n >= 1);
    if (n == 1) return -1;
    if (n == 2) return 1;
    if (n == 4) return 3;
    if (n % 4 == 0) return -1;

    long m = (n % 2 == 0) ? n / 2 : n;
    auto pf_m = factorize(m);
    if (pf_m.size() != 1 || pf_m[0].first == 2) {
        return -1;
    }

    long p = pf_m[0].first;
    int k = pf_m[0].second;
    long phi = m - m / p;

    std::vector<long> phi_primes;
    for (auto [q, e] : factorize(p - 1)) {
        phi_primes.push_back(q);
    }
    if (k >= 2) {
        phi_primes.push_back(p);
    }

    for (long g = 2; g < n; ++g) {
        if (std::gcd(g, n) != 1) continue;
        bool ok = true;
        for (long q : phi_primes) {
            if (pow_mod(g, phi / q, n) == 1) {
                ok = false;
                break;
            }
        }
        if (ok) return g;
    }
    return -1;
}

// 半開区間 [l, r) の全素数（昇順, 区間篩）
inline std::vector<long> primes_in_range(long l, long r) {
    assert(0 <= l && l <= r);
    if (r <= 2) return {};
    l = std::max(l, 2L);
    if (l >= r) return {};

    long limit = detail::isqrt(r - 1);
    constexpr long MAX_SIEVE_LIMIT = 10'000'000L;
    if (limit > MAX_SIEVE_LIMIT) {
        std::vector<long> res;
        for (long x = l; x < r; ++x) {
            if (is_prime(x)) res.push_back(x);
        }
        return res;
    }

    int small_lim = static_cast<int>(limit);
    std::vector<char> is_small_prime(small_lim + 1, 1);
    is_small_prime[0] = is_small_prime[1] = 0;
    std::vector<int> small_primes;
    for (int i = 2; i <= small_lim; ++i) {
        if (is_small_prime[i]) {
            small_primes.push_back(i);
            for (long j = static_cast<long>(i) * i; j <= small_lim; j += i) {
                is_small_prime[static_cast<int>(j)] = 0;
            }
        }
    }

    long len = r - l;
    std::vector<char> is_seg_prime(len, 1);
    for (int p_int : small_primes) {
        long p = p_int;
        long start = std::max(p * p, ((l + p - 1) / p) * p);
        for (long x = start; x < r; x += p) {
            is_seg_prime[x - l] = 0;
        }
    }

    std::vector<long> res;
    for (long i = 0; i < len; ++i) {
        if (is_seg_prime[i]) res.push_back(l + i);
    }
    return res;
}

// n 以上の最小の素数（long の上限を超える場合は -1）
inline long next_prime(long n) {
    if (n <= 2) return 2;
    if (n % 2 == 0) ++n;
    while (true) {
        if (is_prime(n)) return n;
        if (n > std::numeric_limits<long>::max() - 2) return -1;
        n += 2;
    }
}

// n 以下の最大の素数（存在しない場合は -1）
inline long prev_prime(long n) {
    if (n < 2) return -1;
    if (n == 2) return 2;
    if (n % 2 == 0) --n;
    while (n >= 3) {
        if (is_prime(n)) return n;
        n -= 2;
    }
    return 2;
}

}  // namespace prime
