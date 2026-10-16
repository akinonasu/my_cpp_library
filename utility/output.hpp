#include <iterator>
#include <optional>
#include <print>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace output_detail {

template <typename T>
struct is_string_like {
   private:
    using U = std::remove_cvref_t<T>;

    static constexpr bool is_char_array = std::is_array_v<U> && std::is_same_v<std::remove_cv_t<std::remove_extent_t<U>>, char>;

    static constexpr bool is_char_pointer = std::is_pointer_v<U> && std::is_same_v<std::remove_cv_t<std::remove_pointer_t<U>>, char>;

   public:
    static constexpr bool value = std::is_same_v<U, std::string> || std::is_same_v<U, std::string_view> || is_char_array || is_char_pointer;
};

template <typename T>
inline constexpr bool is_string_like_v = is_string_like<T>::value;

template <typename T, typename = void>
struct is_iterable : std::false_type {};

template <typename T>
struct is_iterable<T, std::void_t<decltype(std::begin(std::declval<const T&>())), decltype(std::end(std::declval<const T&>()))>>
    : std::bool_constant<!is_string_like_v<T>> {};

template <typename T>
inline constexpr bool is_iterable_v = is_iterable<T>::value;

}  // namespace output_detail

template <typename T>
struct Plus1Wrapper {
    const T& val;
};

namespace output_detail {

template <typename T, typename = void>
struct is_plus1_supported : std::false_type {};

template <typename T>
struct is_plus1_supported<
    T, std::enable_if_t<std::is_arithmetic_v<std::remove_cvref_t<T>> && !std::is_same_v<std::remove_cvref_t<T>, bool> &&
                        !std::is_same_v<std::remove_cvref_t<T>, char> && !std::is_same_v<std::remove_cvref_t<T>, signed char> &&
                        !std::is_same_v<std::remove_cvref_t<T>, unsigned char> && !std::is_same_v<std::remove_cvref_t<T>, wchar_t> &&
#if defined(__cpp_char8_t)
                        !std::is_same_v<std::remove_cvref_t<T>, char8_t> &&
#endif
                        !std::is_same_v<std::remove_cvref_t<T>, char16_t> && !std::is_same_v<std::remove_cvref_t<T>, char32_t>>>
    : std::true_type {
};

template <typename T1, typename T2>
struct is_plus1_supported<std::pair<T1, T2>> : std::bool_constant<is_plus1_supported<T1>::value && is_plus1_supported<T2>::value> {};

template <typename... Args>
struct is_plus1_supported<std::tuple<Args...>> : std::bool_constant<(is_plus1_supported<Args>::value && ...)> {};

template <typename T>
struct is_plus1_supported<T, std::void_t<decltype(std::begin(std::declval<const T&>())), decltype(std::end(std::declval<const T&>())),
                                         decltype(*std::begin(std::declval<const T&>()))>>
    : std::bool_constant<!is_string_like_v<T> &&
                         is_plus1_supported<std::remove_cvref_t<decltype(*std::begin(std::declval<const T&>()))>>::value> {};

template <typename T>
inline constexpr bool is_plus1_supported_v = is_plus1_supported<std::remove_cvref_t<T>>::value;

}  // namespace output_detail

template <typename T>
    requires output_detail::is_plus1_supported_v<T>
constexpr Plus1Wrapper<T> plus1(const T& x) {
    return Plus1Wrapper<T>{x};
}

template <typename T>
void print_val(const T& val);

template <typename T, typename U>
void print_val(const std::pair<T, U>& p);

template <typename... Args>
void print_val(const std::tuple<Args...>& t);

template <typename T>
void print_val(const std::optional<T>& opt);

template <typename T, typename U>
void print_val(const Plus1Wrapper<std::pair<T, U>>& w);

template <typename... Args>
void print_val(const Plus1Wrapper<std::tuple<Args...>>& w);

template <typename T>
void print_val(const Plus1Wrapper<T>& w);

namespace output_detail {

template <typename Func, typename... Args>
void print_separated(Func&& print_one, const Args&... args) {
    bool first = true;

    auto emit = [&](const auto& x) {
        if (!first) {
            std::print(" ");
        }

        first = false;
        print_one(x);
    };

    (emit(args), ...);
}

}  // namespace output_detail

template <typename T>
void print_val(const T& val) {
    if constexpr (output_detail::is_iterable_v<T>) {
        bool first = true;

        for (const auto& item : val) {
            if (!first) {
                std::print(" ");
            }

            first = false;
            print_val(item);
        }
    } else {
        std::print("{}", val);
    }
}

template <typename T, typename U>
void print_val(const std::pair<T, U>& p) {
    print_val(p.first);
    std::print(" ");
    print_val(p.second);
}

namespace output_detail {

template <typename Tuple, std::size_t... I>
void print_tuple_impl(const Tuple& t, std::index_sequence<I...>) {
    print_separated([](const auto& x) { print_val(x); }, std::get<I>(t)...);
}

}  // namespace output_detail

template <typename... Args>
void print_val(const std::tuple<Args...>& t) {
    output_detail::print_tuple_impl(t, std::index_sequence_for<Args...>{});
}

template <typename T>
void print_val(const std::optional<T>& opt) {
    if (opt.has_value()) {
        print_val(*opt);
    } else {
        std::print("-1");
    }
}

template <typename T, typename U>
void print_val(const Plus1Wrapper<std::pair<T, U>>& w) {
    print_val(plus1(w.val.first));
    std::print(" ");
    print_val(plus1(w.val.second));
}

namespace output_detail {

template <typename Tuple, std::size_t... I>
void print_plus1_tuple_impl(const Tuple& t, std::index_sequence<I...>) {
    print_separated([](const auto& x) { print_val(plus1(x)); }, std::get<I>(t)...);
}

}  // namespace output_detail

template <typename... Args>
void print_val(const Plus1Wrapper<std::tuple<Args...>>& w) {
    output_detail::print_plus1_tuple_impl(w.val, std::index_sequence_for<Args...>{});
}

template <typename T>
void print_val(const Plus1Wrapper<T>& w) {
    if constexpr (output_detail::is_iterable_v<T>) {
        bool first = true;

        for (const auto& item : w.val) {
            if (!first) {
                std::print(" ");
            }

            first = false;
            print_val(plus1(item));
        }
    } else {
        print_val(w.val + 1);
    }
}

template <typename Container>
void out_1indexed(const Container& c) {
    auto it = std::begin(c);

    if (it != std::end(c)) {
        ++it;
    }

    bool first = true;

    for (; it != std::end(c); ++it) {
        if (!first) {
            std::print(" ");
        }

        first = false;
        print_val(*it);
    }

    std::println();
}

template <typename Container>
void out_from1(const Container& c) {
    out_1indexed(c);
}

inline void out() {
    std::println();
}

template <typename... Args>
void out(const Args&... args) {
    output_detail::print_separated([](const auto& x) { print_val(x); }, args...);

    std::println();
}

template <typename... Args>
void pr(const Args&... args) {
    out(args...);
}

template <typename Container>
void out_lines(const Container& c) {
    for (const auto& x : c) {
        print_val(x);
        std::println();
    }
}

template <typename Grid>
void out_grid(const Grid& grid) {
    for (const auto& row : grid) {
        out(row);
    }
}

inline void Yes(bool condition = true, std::string_view true_str = "Yes", std::string_view false_str = "No") {
    std::println("{}", condition ? true_str : false_str);
}

inline void No(bool condition = true) {
    Yes(!condition);
}
