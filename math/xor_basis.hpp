#include <algorithm>
#include <array>
#include <concepts>
#include <optional>
#include <print>
#include <vector>

/**
 * @brief XOR基底 (Linear Basis) ライブラリ
 * @tparam T 扱う整数型 (int または long)
 * @tparam BITS 扱う最大ビット数 (デフォルト: 60)
 */
template <std::integral T = long, int BITS = 60>
struct XorBasis {
    std::array<T, BITS> basis{};
    int sz = 0;  // 基底のサイズ (rank)

    XorBasis() = default;

    /**
     * @brief 要素 x を基底に追加する
     * @return 新たに基底が拡張された場合 true、既存の基底で表現可能な場合 false
     */
    bool insert(T x) {
        for (int i = BITS - 1; i >= 0; --i) {
            if ((x >> i) & 1) {
                if (!basis[i]) {
                    basis[i] = x;
                    ++sz;
                    return true;
                }
                x ^= basis[i];
            }
        }
        return false;
    }

    /**
     * @brief 要素 x が現在の基底の部分集合のXORで表現可能かを判定する
     */
    bool contains(T x) const {
        for (int i = BITS - 1; i >= 0; --i) {
            if ((x >> i) & 1) {
                if (!basis[i]) return false;
                x ^= basis[i];
            }
        }
        return true;
    }

    /**
     * @brief initial とのXORの最大値を求める
     */
    T max_xor(T initial = 0) const {
        T res = initial;
        for (int i = BITS - 1; i >= 0; --i) {
            if ((res ^ basis[i]) > res) {
                res ^= basis[i];
            }
        }
        return res;
    }

    /**
     * @brief initial とのXORの最小値を求める
     */
    T min_xor(T initial = 0) const {
        T res = initial;
        for (int i = BITS - 1; i >= 0; --i) {
            if ((res ^ basis[i]) < res) {
                res ^= basis[i];
            }
        }
        return res;
    }

    /**
     * @brief 行基本変形を行い、簡約化された基底 (各基底のMSB位置が他では0) を返す
     */
    std::vector<T> get_reduced_basis() const {
        auto b = basis;
        for (int i = 0; i < BITS; ++i) {
            if (!b[i]) continue;
            for (int j = i + 1; j < BITS; ++j) {
                if ((b[j] >> i) & 1) {
                    b[j] ^= b[i];
                }
            }
        }
        std::vector<T> res;
        res.reserve(sz);
        for (int i = 0; i < BITS; ++i) {
            if (b[i]) res.push_back(b[i]);
        }
        return res;
    }

    std::vector<T> all_values() const {
        if (sz > 22) {
            // 2^23 ≒ 8.3 * 10^6 要素。これを超える場合は競技プログラミングでは危険
            return {};
        }

        auto rref = get_reduced_basis();
        std::vector<T> res;
        res.reserve(1L << sz);
        res.push_back(0);  // 空集合のXOR値

        for (T b : rref) {
            int cur_sz = res.size();
            for (int i = 0; i < cur_sz; ++i) {
                res.push_back(res[i] ^ b);
            }
        }
        return res;
    }

    /**
     * @brief 作成可能なXOR値のうち、相異なる値の中で k 番目に小さい値を求める (1-indexed)
     * @note 空集合のXOR値 0 は 1 番目として扱われます。
     * @return 存在すればその値、範囲外なら std::nullopt
     */
    std::optional<T> kth_smallest(long k) const {
        if (k < 1) return std::nullopt;
        if (sz < 62 && k > (1L << sz)) return std::nullopt;

        auto rref = get_reduced_basis();
        long idx = k - 1;
        T res = 0;
        for (int i = 0; i < sz; ++i) {
            if ((idx >> i) & 1) {
                res ^= rref[i];
            }
        }
        return res;
    }

    /**
     * @brief 作成可能なXOR値のうち、相異なる値の中で k 番目に大きい値を求める (1-indexed)
     */
    std::optional<T> kth_largest(long k) const {
        if (k < 1) return std::nullopt;
        if (sz < 62 && k > (1L << sz)) return std::nullopt;

        long total = 1L << sz;
        return kth_smallest(total - k + 1);
    }

    /**
     * @brief 別の基底 other をマージする
     */
    void merge(const XorBasis& other) {
        for (int i = BITS - 1; i >= 0; --i) {
            if (other.basis[i]) {
                insert(other.basis[i]);
            }
        }
    }

    XorBasis operator+(const XorBasis& other) const {
        XorBasis res = *this;
        res.merge(other);
        return res;
    }

    XorBasis& operator+=(const XorBasis& other) {
        this->merge(other);
        return *this;
    }

    int size() const { return sz; }
    int rank() const { return sz; }
    bool empty() const { return sz == 0; }
};
// 型エイリアス (32bit用と64bit用)
using XorBasis32 = XorBasis<int, 30>;   // 10^9 程度まで
using XorBasis64 = XorBasis<long, 60>;  // 10^18 程度まで
