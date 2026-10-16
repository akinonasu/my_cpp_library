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

template <std::integral T>
inline T div_ceil(T a, T b) {
    if (a > 0) return a / b + (a % b != 0);
    return a / b;
}

template <std::integral T>
inline T div_floor(T a, T b) {
    if (a < 0) return a / b - (a % b != 0);
    return a / b;
}

template <std::integral T>
inline T mod(T a, T m) {
    a %= m;
    if (a < 0) a += m;
    return a;
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
