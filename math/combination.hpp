#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

#include "modint.hpp"

template <typename Mint>
class Combination {
   private:
    std::vector<Mint> _fact;
    std::vector<Mint> _inv_fact;
    std::vector<Mint> _montmort_cache;

    // 内部テーブルを指定サイズまで拡張
    void extend(std::size_t n) {
        std::size_t old_sz = _fact.size();
        if (n < old_sz) return;
        std::size_t new_sz = std::max(n + 1, old_sz * 2);

        _fact.resize(new_sz);
        _inv_fact.resize(new_sz);

        for (std::size_t i = old_sz; i < new_sz; ++i) {
            _fact[i] = _fact[i - 1] * i;
        }
        _inv_fact[new_sz - 1] = _fact[new_sz - 1].inv();
        for (std::size_t i = new_sz - 1; i >= old_sz; --i) {
            _inv_fact[i - 1] = _inv_fact[i] * i;
        }
    }

   public:
    explicit Combination(std::size_t initial_size = 2) {
        std::size_t sz = std::max<std::size_t>(2, initial_size);
        _fact.assign(sz, 1);
        _inv_fact.assign(sz, 1);
        for (std::size_t i = 1; i < sz; ++i)
            _fact[i] = _fact[i - 1] * i;
        _inv_fact[sz - 1] = _fact[sz - 1].inv();
        for (std::size_t i = sz - 1; i >= 1; --i)
            _inv_fact[i - 1] = _inv_fact[i] * i;
    }

    // 階乗: n!
    Mint fact(std::int64_t n) {
        if (n < 0) return 0;
        extend(n);
        return _fact[n];
    }

    // 階乗の逆元: 1 / n!
    Mint inv_fact(std::int64_t n) {
        if (n < 0) return 0;
        extend(n);
        return _inv_fact[n];
    }

    // 単体の逆元: 1 / n (O(1))
    Mint inv(std::int64_t n) {
        assert(n > 0);
        extend(n);
        return _inv_fact[n] * _fact[n - 1];
    }

    // 順列: nPr
    Mint P(std::int64_t n, std::int64_t r) {
        if (r < 0 || n < r) return 0;
        extend(n);
        return _fact[n] * _inv_fact[n - r];
    }

    // 組合せ: nCr
    Mint C(std::int64_t n, std::int64_t r) {
        if (r < 0 || n < r) return 0;
        extend(n);
        return _fact[n] * _inv_fact[r] * _inv_fact[n - r];
    }

    // 重複組合せ: nHr = (n + r - 1) C r
    Mint H(std::int64_t n, std::int64_t r) {
        if (n < 0 || r < 0) return 0;
        if (n == 0 && r == 0) return 1;
        if (n == 0) return 0;
        return C(n + r - 1, r);
    }

    // カタラン数: C_n = (2n)! / ( (n+1)! * n! )
    Mint catalan(std::int64_t n) {
        if (n < 0) return 0;
        return C(2 * n, n) - C(2 * n, n - 1);
    }

    // 攪乱順列 (モンモール数): 完全順列の総数
    // D(n) = (n - 1) * (D(n - 1) + D(n - 2))
    Mint montmort(std::int64_t n) {
        if (n < 0) return 0;
        if (n == 0) return 1;
        if (n == 1) return 0;

        std::size_t old_sz = _montmort_cache.size();
        if (old_sz == 0) {
            _montmort_cache = {1, 0, 1};  // n = 0, 1, 2
            old_sz = 3;
        }

        if (static_cast<std::size_t>(n) < old_sz) {
            return _montmort_cache[n];
        }

        _montmort_cache.resize(n + 1);
        for (std::size_t i = old_sz; i <= static_cast<std::size_t>(n); ++i) {
            _montmort_cache[i] = Mint(i - 1) * (_montmort_cache[i - 1] + _montmort_cache[i - 2]);
        }
        return _montmort_cache[n];
    }

    // 第2種スターリング数: S(n, k)
    // n 個の区別できる要素を k 個の空でないグループに分割する方法の数 (O(k log n))
    Mint stirling2(std::int64_t n, std::int64_t k) {
        if (n < 0 || k < 0 || n < k) return 0;
        if (n == 0 && k == 0) return 1;

        Mint res = 0;
        for (std::int64_t j = 0; j <= k; ++j) {
            Mint term = C(k, j) * Mint(j).pow(n);
            if ((k - j) & 1) {
                res -= term;
            } else {
                res += term;
            }
        }
        return res * inv_fact(k);
    }

    // N が巨大（最大 10^18）かつ K が小さい（最大 10^6 程度）場合の nCr (O(k))
    Mint C_large(std::int64_t n, std::int64_t k) {
        if (k < 0 || n < k) return 0;
        extend(k);
        Mint num = 1;
        for (std::int64_t i = 0; i < k; ++i) {
            num *= (n - i);
        }
        return num * _inv_fact[k];
    }
};