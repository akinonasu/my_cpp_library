#include <istream>
#include <numeric>
#include <print>
#include <ranges>
#include <vector>

#define ALL(a) (a).begin(), (a).end()
using i128 = __int128;

template <typename T, typename U>
inline bool chmin(T& a, const U& b) {
    if (a > b) {
        a = b;
        return true;
    }
    return false;
}

template <typename T, typename U>
inline bool chmax(T& a, const U& b) {
    if (a < b) {
        a = b;
        return true;
    }
    return false;
}

template <typename T, typename U>
using common_int_t =
    std::conditional_t<std::is_signed_v<T> || std::is_signed_v<U>, std::make_signed_t<std::common_type_t<T, U>>, std::common_type_t<T, U>>;

template <std::integral T, std::integral U>
[[nodiscard]] constexpr auto div_floor(T a, U b) {
    using C = common_int_t<T, U>;
    C ca = static_cast<C>(a);
    C cb = static_cast<C>(b);
    C q = ca / cb;
    C r = ca % cb;
    if (r != 0 && ((ca < 0) != (cb < 0))) {
        --q;
    }
    return q;
}

template <std::integral T, std::integral U>
[[nodiscard]] constexpr auto div_ceil(T a, U b) {
    using C = common_int_t<T, U>;
    C ca = static_cast<C>(a);
    C cb = static_cast<C>(b);
    C q = ca / cb;
    C r = ca % cb;
    if (r != 0 && ((ca < 0) == (cb < 0))) {
        ++q;
    }
    return q;
}

template <std::integral T, std::integral U>
[[nodiscard]] constexpr auto mod(T a, U m) {
    using C = common_int_t<T, U>;
    C ca = static_cast<C>(a);
    C cm = static_cast<C>(m);
    C r = ca % cm;
    if (r < 0) {
        r += (cm < 0 ? -cm : cm);
    }
    return r;
}

template <typename T>
inline constexpr T INF = std::numeric_limits<T>::max() / 2;

template <>
inline constexpr float INF<float> = std::numeric_limits<float>::infinity();

template <>
inline constexpr double INF<double> = std::numeric_limits<double>::infinity();

template <>
inline constexpr long double INF<long double> = std::numeric_limits<long double>::infinity();

template <typename T>
auto make_vector(size_t size, T&& initial_value) {
    return std::vector<std::decay_t<T>>(size, std::forward<T>(initial_value));
}

template <typename... Args>
auto make_vector(size_t size, Args&&... args) {
    auto inner = make_vector(std::forward<Args>(args)...);
    return std::vector<decltype(inner)>(size, inner);
}

template <typename T = int>
inline std::vector<T> iota_vec(int n, T start = 0) {
    std::vector<T> v(n);
    std::iota(v.begin(), v.end(), start);
    return v;
}

template <typename T>
inline std::vector<T> doubled_vec(const std::vector<T>& v) {
    std::vector<T> res;
    res.reserve(v.size() * 2);
    res.insert(res.end(), v.begin(), v.end());
    res.insert(res.end(), v.begin(), v.end());
    return res;
}

inline bool is_palindrome(std::string_view s) {
    return std::ranges::equal(s, s | std::views::reverse);
}

constexpr int digit_sum(long n, int base = 10) {
    // 負の数の場合やオーバーフロー(LONG_MIN)を考慮して符号なしに変換
    unsigned long u = (n < 0) ? -static_cast<unsigned long>(n) : static_cast<unsigned long>(n);
    int sum = 0;
    while (u > 0) {
        sum += static_cast<int>(u % base);
        u /= base;
    }
    return sum;
}

// 文字列版: 64bit整数に収まらない巨大な数値の桁和を求める (10進法)
constexpr int digit_sum(std::string_view s) {
    int sum = 0;
    for (char c : s) {
        if ('0' <= c && c <= '9') {
            sum += c - '0';
        }
    }
    return sum;
}
