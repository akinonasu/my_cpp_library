#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

template <typename T>
std::istream& operator>>(std::istream& is, std::vector<T>& v);

template <typename T1, typename T2>
std::istream& operator>>(std::istream& is, std::pair<T1, T2>& p);

template <typename... Args>
std::istream& operator>>(std::istream& is, std::tuple<Args...>& t);

namespace input_detail {

#ifdef LOCAL

inline std::size_t& read_count() {
    static std::size_t count = 0;
    return count;
}

inline bool& failed_read() {
    static bool failed = false;
    return failed;
}

inline bool stdin_is_terminal() {
#if defined(__unix__) || defined(__APPLE__)
    return isatty(fileno(stdin));
#else
    return false;
#endif
}

struct EndOfInputChecker {
    ~EndOfInputChecker() {
        if (stdin_is_terminal() || failed_read() || std::cin.bad()) return;

        std::cin >> std::ws;
        std::string token;
        if (std::cin >> token) {
            std::cerr << "[input error] 入力が余っています。\n"
                      << "  最初に余った値: " << token << '\n'
                      << "  N と M、H と W、辺数 m、配列長 n などを間違えていないか確認してください。\n";
            std::abort();
        }
    }
};

inline void ensure_end_checker() {
    static EndOfInputChecker checker;
    (void)checker;
}

inline void begin_input() {
    ensure_end_checker();
}

template <typename T>
void read_one(std::istream& is, T& value) {
    if (&is != &std::cin) {
        is >> value;
        return;
    }
    ensure_end_checker();
    ++read_count();
    if (!(is >> value)) {
        failed_read() = true;
        std::cerr << "[input error] 入力が足りないか、型が合いません。\n"
                  << "  " << read_count() << " 個目の値を読み込めませんでした。\n"
                  << "  N と M、H と W、辺数 m、配列長 n などを間違えていないか確認してください。\n";
        std::abort();
    }
}

#else

inline void begin_input() {
}

template <typename T>
void read_one(std::istream& is, T& value) {
    is >> value;
}

#endif

template <typename... Args>
void read_values(Args&... args) {
    (read_one(std::cin, args), ...);
}

}  // namespace input_detail

struct FastIO {
    FastIO() {
        std::ios_base::sync_with_stdio(false);
        std::cin.tie(nullptr);
    }
};
inline FastIO fast_io_init;

template <typename T1, typename T2>
std::istream& operator>>(std::istream& is, std::pair<T1, T2>& p) {
    input_detail::read_one(is, p.first);
    input_detail::read_one(is, p.second);
    return is;
}

template <typename Tuple, std::size_t... I>
void read_tuple_impl(std::istream& is, Tuple& t, std::index_sequence<I...>) {
    (..., input_detail::read_one(is, std::get<I>(t)));
}

template <typename... Args>
std::istream& operator>>(std::istream& is, std::tuple<Args...>& t) {
    read_tuple_impl(is, t, std::index_sequence_for<Args...>{});
    return is;
}

template <typename T>
std::istream& operator>>(std::istream& is, std::vector<T>& v) {
    for (auto& elem : v) {
        input_detail::read_one(is, elem);
    }
    return is;
}

template <typename T>
struct is_pair : std::false_type {};
template <typename T1, typename T2>
struct is_pair<std::pair<T1, T2>> : std::true_type {};

template <typename T>
struct is_tuple : std::false_type {};
template <typename... Args>
struct is_tuple<std::tuple<Args...>> : std::true_type {};

template <typename T>
struct is_vector : std::false_type {};
template <typename T>
struct is_vector<std::vector<T>> : std::true_type {};

template <typename T>
void adjust_zero_indexed(T& val) {
    using DecayedT = std::decay_t<T>;
    if constexpr (is_pair<DecayedT>::value) {
        adjust_zero_indexed(val.first);
        adjust_zero_indexed(val.second);
    } else if constexpr (is_tuple<DecayedT>::value) {
        std::apply([](auto&... args) { (adjust_zero_indexed(args), ...); }, val);
    } else if constexpr (is_vector<DecayedT>::value) {
        for (auto& elem : val) {
            adjust_zero_indexed(elem);
        }
    } else if constexpr (std::is_arithmetic_v<DecayedT> && !std::is_same_v<DecayedT, char> && !std::is_same_v<DecayedT, signed char> &&
                         !std::is_same_v<DecayedT, unsigned char> && !std::is_same_v<DecayedT, wchar_t> &&
#if defined(__cpp_char8_t)
                         !std::is_same_v<DecayedT, char8_t> &&
#endif
                         !std::is_same_v<DecayedT, char16_t> && !std::is_same_v<DecayedT, char32_t> && !std::is_same_v<DecayedT, bool>) {
        --val;
    } else {
    }
}

using default_type = long;

template <typename... Args>
void read(Args&... args) {
    input_detail::read_values(args...);
}

template <typename T = default_type>
T read_val() {
    T val;
    input_detail::read_one(std::cin, val);
    return val;
}

template <typename T1 = default_type, typename T2 = default_type>
std::pair<T1, T2> read_pair() {
    std::pair<T1, T2> p;
    input_detail::begin_input();
    std::cin >> p;
    return p;
}

template <typename... Args>
std::tuple<Args...> read_tuple() {
    std::tuple<Args...> t;
    input_detail::begin_input();
    std::cin >> t;
    return t;
}

template <typename T = default_type, typename Size, std::enable_if_t<std::is_integral_v<Size> && !std::is_same_v<Size, bool>, int> = 0>
std::vector<T> read_vec(Size n, bool zero_indexed = false) {
    std::vector<T> v(n);
    input_detail::begin_input();
    std::cin >> v;
    if (zero_indexed) {
        adjust_zero_indexed(v);
    }
    return v;
}

template <typename T = default_type>
std::vector<T> read_vec(bool zero_indexed = false) {
    int n;
    input_detail::read_one(std::cin, n);
    return read_vec<T>(n, zero_indexed);
}

template <typename T1 = default_type, typename T2 = default_type, typename Size,
          std::enable_if_t<std::is_integral_v<Size> && !std::is_same_v<Size, bool>, int> = 0>
std::vector<std::pair<T1, T2>> read_vec_pair(Size n, bool zero_indexed = false) {
    return read_vec<std::pair<T1, T2>>(n, zero_indexed);
}

template <typename T = default_type, typename Size, std::enable_if_t<std::is_integral_v<Size> && !std::is_same_v<Size, bool>, int> = 0>
std::vector<T> read_sorted_vec(Size n, bool zero_indexed = false) {
    auto v = read_vec<T>(n, zero_indexed);
    std::sort(v.begin(), v.end());
    return v;
}

template <typename T = default_type>
std::vector<T> read_sorted_vec(bool zero_indexed = false) {
    auto v = read_vec<T>(zero_indexed);
    std::sort(v.begin(), v.end());
    return v;
}

template <typename T = default_type, typename Size, typename Compare,
          std::enable_if_t<std::is_integral_v<Size> && !std::is_same_v<Size, bool>, int> = 0>
std::vector<T> read_sorted_vec(Size n, Compare comp, bool zero_indexed = false) {
    auto v = read_vec<T>(n, zero_indexed);
    std::sort(v.begin(), v.end(), comp);
    return v;
}

template <typename T = default_type, typename Compare>
std::vector<T> read_sorted_vec(Compare comp, bool zero_indexed = false) {
    auto v = read_vec<T>(zero_indexed);
    std::sort(v.begin(), v.end(), comp);
    return v;
}

template <typename T = default_type, typename Size, std::enable_if_t<std::is_integral_v<Size> && !std::is_same_v<Size, bool>, int> = 0>
std::vector<std::pair<T, int>> read_indexed_vec(Size n, bool zero_indexed = false) {
    std::vector<std::pair<T, int>> v(n);
    input_detail::begin_input();
    for (int i = 0; i < static_cast<int>(n); ++i) {
        input_detail::read_one(std::cin, v[i].first);
        v[i].second = i;
    }
    if (zero_indexed) {
        for (auto& elem : v) {
            adjust_zero_indexed(elem.first);
        }
    }
    return v;
}

template <typename T = default_type>
std::vector<std::pair<T, int>> read_indexed_vec(bool zero_indexed = false) {
    int n;
    input_detail::read_one(std::cin, n);
    return read_indexed_vec<T>(n, zero_indexed);
}

template <typename T = default_type, typename Size, std::enable_if_t<std::is_integral_v<Size> && !std::is_same_v<Size, bool>, int> = 0>
std::vector<std::pair<T, int>> read_sorted_indexed_vec(Size n, bool zero_indexed = false) {
    auto v = read_indexed_vec<T>(n, zero_indexed);
    std::stable_sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    return v;
}

template <typename T = default_type>
std::vector<std::pair<T, int>> read_sorted_indexed_vec(bool zero_indexed = false) {
    auto v = read_indexed_vec<T>(zero_indexed);
    std::stable_sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    return v;
}

template <
    typename T = default_type, typename Size, typename Compare,
    std::enable_if_t<std::is_integral_v<Size> && !std::is_same_v<Size, bool> && !std::is_same_v<std::decay_t<Compare>, bool>, int> = 0>
std::vector<std::pair<T, int>> read_sorted_indexed_vec(Size n, Compare comp, bool zero_indexed = false) {
    auto v = read_indexed_vec<T>(n, zero_indexed);
    std::stable_sort(v.begin(), v.end(), [&](const auto& a, const auto& b) { return comp(a.first, b.first); });
    return v;
}

template <typename T = default_type, typename Compare, std::enable_if_t<!std::is_same_v<std::decay_t<Compare>, bool>, int> = 0>
std::vector<std::pair<T, int>> read_sorted_indexed_vec(Compare comp, bool zero_indexed = false) {
    auto v = read_indexed_vec<T>(zero_indexed);
    std::stable_sort(v.begin(), v.end(), [&](const auto& a, const auto& b) { return comp(a.first, b.first); });
    return v;
}

template <typename T1 = default_type, typename T2 = default_type>
std::vector<std::pair<T1, T2>> read_vec_pair(bool zero_indexed = false) {
    int n;
    input_detail::read_one(std::cin, n);
    return read_vec_pair<T1, T2>(n, zero_indexed);
}

template <typename... Args, typename Size, std::enable_if_t<std::is_integral_v<Size> && !std::is_same_v<Size, bool>, int> = 0>
std::vector<std::tuple<Args...>> read_vec_tuple(Size n, bool zero_indexed = false) {
    return read_vec<std::tuple<Args...>>(n, zero_indexed);
}

template <typename... Args>
std::vector<std::tuple<Args...>> read_vec_tuple(bool zero_indexed = false) {
    int n;
    input_detail::read_one(std::cin, n);
    return read_vec_tuple<Args...>(n, zero_indexed);
}

template <typename T = default_type>
std::vector<std::vector<T>> read_vec_grid(int h, int w, bool zero_indexed = false) {
    std::vector<std::vector<T>> grid(h, std::vector<T>(w));
    input_detail::begin_input();
    std::cin >> grid;
    if (zero_indexed) {
        adjust_zero_indexed(grid);
    }
    return grid;
}

template <typename T = default_type>
std::vector<std::vector<T>> read_vec_grid(bool zero_indexed = false) {
    int h, w;
    input_detail::read_values(h, w);
    return read_vec_grid<T>(h, w, zero_indexed);
}

template <typename T = default_type, typename Size, std::enable_if_t<std::is_integral_v<Size> && !std::is_same_v<Size, bool>, int> = 0>
std::vector<std::vector<T>> read_vec_var(Size n, bool zero_indexed = false) {
    std::vector<std::vector<T>> res(n);
    input_detail::begin_input();
    for (int i = 0; i < static_cast<int>(n); ++i) {
        int m;
        input_detail::read_one(std::cin, m);
        res[i] = read_vec<T>(m, zero_indexed);
    }
    return res;
}

template <typename T = default_type>
std::vector<std::vector<T>> read_vec_var(bool zero_indexed = false) {
    int n;
    input_detail::read_one(std::cin, n);
    return read_vec_var<T>(n, zero_indexed);
}

template <typename T = default_type>
T read_zero_idx() {
    T val;
    input_detail::read_one(std::cin, val);
    adjust_zero_indexed(val);
    return val;
}

inline std::vector<std::vector<int>> read_graph(int n, int m, bool directed = false) {
    std::vector<std::vector<int>> g(n);
    input_detail::begin_input();
    for (int i = 0; i < m; ++i) {
        int u = read_zero_idx<int>();
        int v = read_zero_idx<int>();
        g[u].push_back(v);
        if (!directed) {
            g[v].push_back(u);
        }
    }
    return g;
}

inline std::vector<std::vector<int>> read_graph(bool directed = false) {
    int n, m;
    input_detail::read_values(n, m);
    return read_graph(n, m, directed);
}

template <typename Cost = default_type>
struct Edge {
    int to;
    Cost cost;

    Edge() = default;
    Edge(int to, Cost cost) : to(to), cost(cost) {}
};

template <typename Cost = default_type>
inline std::vector<std::vector<Edge<Cost>>> read_weighted_graph(int n, int m, bool directed = false) {
    std::vector<std::vector<Edge<Cost>>> g(n);
    input_detail::begin_input();
    for (int i = 0; i < m; ++i) {
        int u = read_zero_idx<int>();
        int v = read_zero_idx<int>();
        Cost w;
        input_detail::read_one(std::cin, w);
        g[u].push_back(Edge<Cost>{v, w});
        if (!directed) {
            g[v].push_back(Edge<Cost>{u, w});
        }
    }
    return g;
}

template <typename Cost = default_type>
inline std::vector<std::vector<Edge<Cost>>> read_weighted_graph(bool directed = false) {
    int n, m;
    input_detail::read_values(n, m);
    return read_weighted_graph<Cost>(n, m, directed);
}
