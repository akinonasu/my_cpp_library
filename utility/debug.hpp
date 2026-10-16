#include <algorithm>
#include <bitset>
#include <concepts>
#include <format>
#include <optional>
#include <print>
#include <queue>
#include <ranges>
#include <stack>
#include <string_view>
#include <utility>
#include <vector>

// =========================================================================
// 1. 標準で未対応なコンテナの formatter 拡張
// =========================================================================
template <typename T, typename C>
struct std::formatter<std::queue<T, C>> : std::formatter<std::string_view> {
    auto format(const std::queue<T, C>& q_orig, std::format_context& ctx) const {
        auto q = q_orig;
        std::vector<T> v;
        v.reserve(q.size());
        while (!q.empty()) {
            v.push_back(q.front());
            q.pop();
        }
        return std::formatter<std::string_view>::format(std::format("queue{}", v), ctx);
    }
};

template <typename T, typename C>
struct std::formatter<std::stack<T, C>> : std::formatter<std::string_view> {
    auto format(const std::stack<T, C>& s_orig, std::format_context& ctx) const {
        auto s = s_orig;
        std::vector<T> v;
        v.reserve(s.size());
        while (!s.empty()) {
            v.push_back(s.top());
            s.pop();
        }
        std::ranges::reverse(v);
        return std::formatter<std::string_view>::format(std::format("stack{}", v), ctx);
    }
};

template <typename T, typename C, typename Comp>
struct std::formatter<std::priority_queue<T, C, Comp>> : std::formatter<std::string_view> {
    auto format(const std::priority_queue<T, C, Comp>& pq_orig, std::format_context& ctx) const {
        auto pq = pq_orig;
        std::vector<T> v;
        v.reserve(pq.size());
        while (!pq.empty()) {
            v.push_back(pq.top());
            pq.pop();
        }
        return std::formatter<std::string_view>::format(std::format("priority_queue{}", v), ctx);
    }
};

template <typename T>
struct std::formatter<std::optional<T>> : std::formatter<std::string_view> {
    auto format(const std::optional<T>& opt, std::format_context& ctx) const {
        if (opt) {
            return std::formatter<std::string_view>::format(std::format("opt({})", *opt), ctx);
        } else {
            return std::formatter<std::string_view>::format("nullopt", ctx);
        }
    }
};

template <size_t N>
struct std::formatter<std::bitset<N>> : std::formatter<std::string_view> {
    auto format(const std::bitset<N>& b, std::format_context& ctx) const {
        return std::formatter<std::string_view>::format(b.to_string(), ctx);
    }
};

// =========================================================================
// 2. 2次元レンジ（グリッド）判定 Concept & 出力ロジック
// =========================================================================
namespace debug_internal {
template <typename T>
concept is_string_like = std::is_convertible_v<T, std::string_view>;

// 「自身がレンジ」かつ「文字列ではない」かつ「要素もレンジ」である型を2次元グリッドと判定
template <typename T>
concept is_2d_range =
    std::ranges::range<T> && !is_string_like<T> && requires { requires std::ranges::range<std::ranges::range_reference_t<T>>; };

// カンマ区切りの変数名パース
inline std::vector<std::string_view> split_names(std::string_view s) {
    std::vector<std::string_view> names;
    int paren = 0;
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
            else if (c == ',' && paren == 0) {
                auto sub = s.substr(start, i - start);
                while (!sub.empty() && (sub.front() == ' ' || sub.front() == '\t'))
                    sub.remove_prefix(1);
                while (!sub.empty() && (sub.back() == ' ' || sub.back() == '\t'))
                    sub.remove_suffix(1);
                names.push_back(sub);
                start = i + 1;
            }
        }
    }
    auto sub = s.substr(start);
    while (!sub.empty() && (sub.front() == ' ' || sub.front() == '\t'))
        sub.remove_prefix(1);
    while (!sub.empty() && (sub.back() == ' ' || sub.back() == '\t'))
        sub.remove_suffix(1);
    names.push_back(sub);
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

        // 改行位置の自動調整
        if (idx == 0) {
            std::print(stderr, "{} = ", name);
        } else {
            if (is_2d || prev_was_2d) {
                std::print(stderr, "\n  {} = ", name);
            } else {
                std::print(stderr, ", {} = ", name);
            }
        }

        // 2次元グリッドの場合のフォーマット出力
        if constexpr (is_2d) {
            if (std::ranges::empty(val)) {
                std::print(stderr, "[]");
            } else {
                std::print(stderr, "[\n");
                int width = 2;
                if constexpr (std::ranges::sized_range<T>) {
                    size_t sz = std::ranges::size(val);
                    if (sz >= 100) width = std::to_string(sz - 1).size();
                }
                size_t r = 0;
                for (const auto& row : val) {
                    // 行番号を揃えて出力: "    [ 0] [1, 2, 3]"
                    std::print(stderr, "    [{:>{}}] {}\n", r++, width, row);
                }
                std::print(stderr, "  ]");
            }
        } else {
            std::print(stderr, "{}", val);
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
// 3. debug マクロ定義
// =========================================================================
#if defined(LOCAL) || !defined(ONLINE_JUDGE)
#define debug(...) debug_internal::print_debug(__LINE__, #__VA_ARGS__ __VA_OPT__(, ) __VA_ARGS__)
#else
#define debug(...) ((void)0)
#endif