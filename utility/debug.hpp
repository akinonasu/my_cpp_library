#include <algorithm>
#include <bitset>
#include <concepts>
#include <iostream>
#include <optional>
#include <print>
#include <queue>
#include <ranges>
#include <stack>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace debug_internal {

// =========================================================================
// 1. 各種 Concept & 前方宣言
// =========================================================================
template <typename T>
concept is_string_like = std::is_convertible_v<T, std::string_view>;

template <typename T>
concept is_2d_range =
    std::ranges::range<T> && !is_string_like<T> && requires { requires std::ranges::range<std::ranges::range_reference_t<T>>; };

// 各種オーバーロードの前方宣言
template <typename T>
void print_val(const T& val);
template <typename T, typename U>
void print_val(const std::pair<T, U>& p);
template <typename... Ts>
void print_val(const std::tuple<Ts...>& t);
template <typename T>
void print_val(const std::optional<T>& opt);
template <typename T, typename C>
void print_val(const std::queue<T, C>& q_orig);
template <typename T, typename C>
void print_val(const std::stack<T, C>& s_orig);
template <typename T, typename C, typename Comp>
void print_val(const std::priority_queue<T, C, Comp>& pq_orig);
template <size_t N>
void print_val(const std::bitset<N>& b);

// =========================================================================
// 2. 個別型の出力実装
// =========================================================================
template <typename T, typename U>
void print_val(const std::pair<T, U>& p) {
    std::print(stderr, "(");
    print_val(p.first);
    std::print(stderr, ", ");
    print_val(p.second);
    std::print(stderr, ")");
}

template <typename Tuple, size_t... Is>
void print_tuple_helper(const Tuple& t, std::index_sequence<Is...>) {
    std::print(stderr, "(");
    size_t i = 0;
    auto print_elem = [&](const auto& elem) {
        if (i++ > 0) std::print(stderr, ", ");
        print_val(elem);
    };
    (print_elem(std::get<Is>(t)), ...);
    std::print(stderr, ")");
}

template <typename... Ts>
void print_val(const std::tuple<Ts...>& t) {
    print_tuple_helper(t, std::index_sequence_for<Ts...>{});
}

template <typename T>
void print_val(const std::optional<T>& opt) {
    if (opt) {
        std::print(stderr, "opt(");
        print_val(*opt);
        std::print(stderr, ")");
    } else {
        std::print(stderr, "nullopt");
    }
}

template <typename T, typename C>
void print_val(const std::queue<T, C>& q_orig) {
    auto q = q_orig;
    std::print(stderr, "queue[");
    bool first = true;
    while (!q.empty()) {
        if (!first) std::print(stderr, ", ");
        first = false;
        print_val(q.front());
        q.pop();
    }
    std::print(stderr, "]");
}

template <typename T, typename C>
void print_val(const std::stack<T, C>& s_orig) {
    auto s = s_orig;
    std::vector<T> v;
    v.reserve(s.size());
    while (!s.empty()) {
        v.push_back(s.top());
        s.pop();
    }
    std::ranges::reverse(v);
    std::print(stderr, "stack[");
    bool first = true;
    for (const auto& elem : v) {
        if (!first) std::print(stderr, ", ");
        first = false;
        print_val(elem);
    }
    std::print(stderr, "]");
}

template <typename T, typename C, typename Comp>
void print_val(const std::priority_queue<T, C, Comp>& pq_orig) {
    auto pq = pq_orig;
    std::print(stderr, "priority_queue[");
    bool first = true;
    while (!pq.empty()) {
        if (!first) std::print(stderr, ", ");
        first = false;
        print_val(pq.top());
        pq.pop();
    }
    std::print(stderr, "]");
}

template <size_t N>
void print_val(const std::bitset<N>& b) {
    std::print(stderr, "{}", b.to_string());
}

// 汎用（レンジ・プリミティブ型）の出力
template <typename T>
void print_val(const T& val) {
    using U = std::remove_cvref_t<T>;
    if constexpr (is_string_like<U>) {
        std::print(stderr, "{}", std::string_view(val));
    } else if constexpr (std::ranges::range<U>) {
        std::print(stderr, "[");
        bool first = true;
        for (const auto& elem : val) {
            if (!first) std::print(stderr, ", ");
            first = false;
            print_val(elem);
        }
        std::print(stderr, "]");
    } else if constexpr (std::is_same_v<U, bool>) {
        std::print(stderr, "{}", val ? "true" : "false");
    } else if constexpr (std::is_same_v<U, char>) {
        std::print(stderr, "'{}'", val);
    } else {
        std::print(stderr, "{}", val);
    }
}

// =========================================================================
// 3. 変数名パーサー & デバッグ出力本体
// =========================================================================
inline std::vector<std::string_view> split_names(std::string_view s) {
    std::vector<std::string_view> names;
    int paren = 0;  // (), [], {} 用
    int angle = 0;  // <> (テンプレート引数) 用
    bool in_str = false, in_char = false;
    size_t start = 0;

    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if (c == '\\' && (in_str || in_char)) {
            i++;
            continue;
        }
        if (c == '"' && !in_char)
            in_str = !in_str;
        else if (c == '\'' && !in_str)
            in_char = !in_char;
        else if (!in_str && !in_char) {
            if (c == '(' || c == '[' || c == '{')
                paren++;
            else if (c == ')' || c == ']' || c == '}')
                paren--;
            else if (c == '<')
                angle++;
            else if (c == '>')
                angle--;
            else if (c == ',' && paren == 0 && angle == 0) {
                auto sub = s.substr(start, i - start);
                while (!sub.empty() && (sub.front() == ' ' || sub.front() == '\t'))
                    sub.remove_prefix(1);
                while (!sub.empty() && (sub.back() == ' ' || sub.back() == '\t'))
                    sub.remove_suffix(1);
                if (!sub.empty()) names.push_back(sub);
                start = i + 1;
            }
        }
    }
    auto sub = s.substr(start);
    while (!sub.empty() && (sub.front() == ' ' || sub.front() == '\t'))
        sub.remove_prefix(1);
    while (!sub.empty() && (sub.back() == ' ' || sub.back() == '\t'))
        sub.remove_suffix(1);
    if (!sub.empty()) names.push_back(sub);
    return names;
}

template <typename... Args>
void print_debug(int line, std::string_view names_str, const Args&... args) {
    auto names = split_names(names_str);
    std::print(stderr, "[Line {:>3}] ", line);

    size_t idx = 0;
    bool prev_was_2d = false;

    auto print_one = [&](const auto& val) {
        using T = std::remove_cvref_t<decltype(val)>;
        std::string_view name = (idx < names.size()) ? names[idx] : "";
        constexpr bool is_2d = is_2d_range<T>;

        if (idx == 0) {
            std::print(stderr, "{} = ", name);
        } else {
            if (is_2d || prev_was_2d) {
                std::print(stderr, "\n  {} = ", name);
            } else {
                std::print(stderr, ", {} = ", name);
            }
        }

        if constexpr (is_2d) {
            if (std::ranges::empty(val)) {
                std::print(stderr, "[]");
            } else {
                std::print(stderr, "[\n");
                int width = 2;
                if constexpr (std::ranges::sized_range<T>) {
                    size_t sz = std::ranges::size(val);
                    if (sz >= 100) width = static_cast<int>(std::to_string(sz - 1).size());
                }
                size_t r = 0;
                for (const auto& row : val) {
                    std::print(stderr, "    [{:>{}}] ", r++, width);
                    print_val(row);
                    std::print(stderr, "\n");
                }
                std::print(stderr, "  ]");
            }
        } else {
            print_val(val);
        }

        prev_was_2d = is_2d;
        idx++;
    };

    (print_one(args), ...);
    std::print(stderr, "\n");
}

inline void print_debug(int line, std::string_view) {
    std::print(stderr, "[Line {:>3}]\n", line);
}
}  // namespace debug_internal

// =========================================================================
// 4. debug マクロ定義
// =========================================================================
#if defined(LOCAL) || !defined(ONLINE_JUDGE)
#define debug(...) debug_internal::print_debug(__LINE__, #__VA_ARGS__ __VA_OPT__(, ) __VA_ARGS__)
#else
#define debug(...) ((void)0)
#endif