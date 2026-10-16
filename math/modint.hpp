#include <concepts>
#include <cstdint>
#include <format>
#include <iostream>
#include <print>
#include <utility>

namespace internal {

// 内部計算用の拡張GCD (ACL準拠: a*x + b*y = gcd(a, b))
constexpr std::pair<std::int64_t, std::int64_t> inv_gcd(std::int64_t a, std::int64_t b) {
    a %= b;
    if (a < 0) a += b;
    if (a == 0) return {b, 0};
    std::int64_t s = b, t = a;
    std::int64_t m0 = 0, m1 = 1;
    while (t) {
        std::int64_t u = s / t;
        s -= t * u;
        m0 -= m1 * u;
        std::swap(s, t);
        std::swap(m0, m1);
    }
    if (m0 < 0) m0 += b / s;
    return {s, m0};
}

// 高速剰余計算のための Barrett Reduction
struct Barrett {
    std::uint32_t _m;
    std::uint64_t im;

    explicit constexpr Barrett(std::uint32_t m = 1) : _m(m), im(static_cast<std::uint64_t>(-1) / m + 1) {}

    [[nodiscard]] constexpr std::uint32_t umod() const noexcept { return _m; }

    [[nodiscard]] constexpr std::uint32_t mul(std::uint32_t a, std::uint32_t b) const noexcept {
        std::uint64_t z = static_cast<std::uint64_t>(a) * b;
        std::uint64_t x = static_cast<std::uint64_t>((static_cast<unsigned __int128>(z) * im) >> 64);
        std::uint64_t v = z - x * _m;
        if (v >= _m) v += _m;
        return static_cast<std::uint32_t>(v >= _m ? v - _m : v);
    }
};

}  // namespace internal

// ==========================================
// 1. 静的 Modint (StaticModInt)
// ==========================================
template <std::uint32_t Mod>
struct StaticModInt {
    static_assert(Mod >= 1, "Modulus must be at least 1");
    static_assert(Mod < (1u << 31), "Modulus must be less than 2^31 to prevent addition overflow");

   private:
    std::uint32_t _v = 0;

   public:
    static constexpr std::uint32_t mod() noexcept { return Mod; }

    constexpr StaticModInt() noexcept : _v(0) {}

    template <std::integral T>
    constexpr StaticModInt(T v) noexcept {
        std::int64_t x = static_cast<std::int64_t>(v % static_cast<std::int64_t>(Mod));
        if (x < 0) x += Mod;
        _v = static_cast<std::uint32_t>(x);
    }

    [[nodiscard]] constexpr std::uint32_t val() const noexcept { return _v; }

    constexpr StaticModInt& operator++() noexcept {
        _v++;
        if (_v == Mod) _v = 0;
        return *this;
    }
    constexpr StaticModInt& operator--() noexcept {
        if (_v == 0) _v = Mod;
        _v--;
        return *this;
    }
    constexpr StaticModInt operator++(int) noexcept {
        StaticModInt res = *this;
        ++*this;
        return res;
    }
    constexpr StaticModInt operator--(int) noexcept {
        StaticModInt res = *this;
        --*this;
        return res;
    }

    constexpr StaticModInt& operator+=(const StaticModInt& rhs) noexcept {
        _v += rhs._v;
        if (_v >= Mod) _v -= Mod;
        return *this;
    }
    constexpr StaticModInt& operator-=(const StaticModInt& rhs) noexcept {
        if (_v < rhs._v) _v += Mod;
        _v -= rhs._v;
        return *this;
    }
    constexpr StaticModInt& operator*=(const StaticModInt& rhs) noexcept {
        _v = static_cast<std::uint32_t>(static_cast<std::uint64_t>(_v) * rhs._v % Mod);
        return *this;
    }
    constexpr StaticModInt& operator/=(const StaticModInt& rhs) noexcept { return *this = *this * rhs.inv(); }

    [[nodiscard]] constexpr StaticModInt operator+() const noexcept { return *this; }
    [[nodiscard]] constexpr StaticModInt operator-() const noexcept { return StaticModInt() - *this; }

    [[nodiscard]] friend constexpr StaticModInt operator+(StaticModInt lhs, const StaticModInt& rhs) noexcept { return lhs += rhs; }
    [[nodiscard]] friend constexpr StaticModInt operator-(StaticModInt lhs, const StaticModInt& rhs) noexcept { return lhs -= rhs; }
    [[nodiscard]] friend constexpr StaticModInt operator*(StaticModInt lhs, const StaticModInt& rhs) noexcept { return lhs *= rhs; }
    [[nodiscard]] friend constexpr StaticModInt operator/(StaticModInt lhs, const StaticModInt& rhs) noexcept { return lhs /= rhs; }

    [[nodiscard]] friend constexpr bool operator==(const StaticModInt& lhs, const StaticModInt& rhs) noexcept { return lhs._v == rhs._v; }

    [[nodiscard]] constexpr StaticModInt pow(std::int64_t n) const noexcept {
        if (n == 0) return 1;
        StaticModInt res = 1, base = *this;
        std::uint64_t exp = (n < 0) ? -static_cast<std::uint64_t>(n) : static_cast<std::uint64_t>(n);
        if (n < 0) base = base.inv();

        while (exp > 0) {
            if (exp & 1) res *= base;
            base *= base;
            exp >>= 1;
        }
        return res;
    }

    [[nodiscard]] constexpr StaticModInt inv() const {
        auto [g, x] = internal::inv_gcd(_v, Mod);
        return StaticModInt(x);
    }

    friend std::istream& operator>>(std::istream& is, StaticModInt& m) {
        std::int64_t val;
        is >> val;
        m = StaticModInt(val);
        return is;
    }
};

// ==========================================
// 2. 動的 Modint (DynamicModInt)
// ==========================================
template <int Id = -1>
struct DynamicModInt {
   private:
    std::uint32_t _v = 0;
    static inline internal::Barrett bt{998244353};  // デフォルト値

   public:
    static void set_mod(std::uint32_t m) {
        // 実行時に 2^31 未満であることを保証
        if (m < 1 || m >= (1u << 31)) {
            std::println(stderr, "Error: Modulo must be between 1 and 2^31 - 1");
            std::abort();
        }
        bt = internal::Barrett(m);
    }

    [[nodiscard]] static std::uint32_t mod() noexcept { return bt.umod(); }

    constexpr DynamicModInt() noexcept : _v(0) {}

    template <std::integral T>
    DynamicModInt(T v) noexcept {
        std::int64_t x = static_cast<std::int64_t>(v % static_cast<std::int64_t>(mod()));
        if (x < 0) x += mod();
        _v = static_cast<std::uint32_t>(x);
    }

    [[nodiscard]] constexpr std::uint32_t val() const noexcept { return _v; }

    DynamicModInt& operator++() noexcept {
        _v++;
        if (_v == mod()) _v = 0;
        return *this;
    }
    DynamicModInt& operator--() noexcept {
        if (_v == 0) _v = mod();
        _v--;
        return *this;
    }
    DynamicModInt operator++(int) noexcept {
        DynamicModInt res = *this;
        ++*this;
        return res;
    }
    DynamicModInt operator--(int) noexcept {
        DynamicModInt res = *this;
        --*this;
        return res;
    }

    DynamicModInt& operator+=(const DynamicModInt& rhs) noexcept {
        _v += rhs._v;
        if (_v >= mod()) _v -= mod();
        return *this;
    }
    DynamicModInt& operator-=(const DynamicModInt& rhs) noexcept {
        if (_v < rhs._v) _v += mod();
        _v -= rhs._v;
        return *this;
    }
    DynamicModInt& operator*=(const DynamicModInt& rhs) noexcept {
        _v = bt.mul(_v, rhs._v);
        return *this;
    }
    DynamicModInt& operator/=(const DynamicModInt& rhs) noexcept { return *this = *this * rhs.inv(); }

    [[nodiscard]] constexpr DynamicModInt operator+() const noexcept { return *this; }
    [[nodiscard]] DynamicModInt operator-() const noexcept { return DynamicModInt() - *this; }

    [[nodiscard]] friend DynamicModInt operator+(DynamicModInt lhs, const DynamicModInt& rhs) noexcept { return lhs += rhs; }
    [[nodiscard]] friend DynamicModInt operator-(DynamicModInt lhs, const DynamicModInt& rhs) noexcept { return lhs -= rhs; }
    [[nodiscard]] friend DynamicModInt operator*(DynamicModInt lhs, const DynamicModInt& rhs) noexcept { return lhs *= rhs; }
    [[nodiscard]] friend DynamicModInt operator/(DynamicModInt lhs, const DynamicModInt& rhs) noexcept { return lhs /= rhs; }

    [[nodiscard]] friend constexpr bool operator==(const DynamicModInt& lhs, const DynamicModInt& rhs) noexcept { return lhs._v == rhs._v; }

    [[nodiscard]] DynamicModInt pow(std::int64_t n) const noexcept {
        if (n < 0) return inv().pow(-n);
        DynamicModInt res = 1, base = *this;
        while (n > 0) {
            if (n & 1) res *= base;
            base *= base;
            n >>= 1;
        }
        return res;
    }

    [[nodiscard]] DynamicModInt inv() const {
        auto [g, x] = internal::inv_gcd(_v, mod());
        return DynamicModInt(x);
    }

    friend std::istream& operator>>(std::istream& is, DynamicModInt& m) {
        std::int64_t val;
        is >> val;
        m = DynamicModInt(val);
        return is;
    }
};

// ==========================================
// 3. std::print / std::format 対応
// ==========================================
template <std::uint32_t Mod>
struct std::formatter<StaticModInt<Mod>> : std::formatter<std::uint32_t> {
    auto format(const StaticModInt<Mod>& m, std::format_context& ctx) const { return std::formatter<std::uint32_t>::format(m.val(), ctx); }
};

template <int Id>
struct std::formatter<DynamicModInt<Id>> : std::formatter<std::uint32_t> {
    auto format(const DynamicModInt<Id>& m, std::format_context& ctx) const { return std::formatter<std::uint32_t>::format(m.val(), ctx); }
};

// 型エイリアス
using mint998244353 = StaticModInt<998244353>;
using mint1000000007 = StaticModInt<1000000007>;
using mint = mint998244353;

using dmint = DynamicModInt<-1>;  // デフォルトの動的modint

/*
    実行時に Mod を設定
    dmint::set_mod(13);

    複数の異なる動的 mod を併用する例
    using dmintA = DynamicModInt<0>;
    using dmintB = DynamicModInt<1>;

    dmintA::set_mod(7);
    dmintB::set_mod(11);

*/