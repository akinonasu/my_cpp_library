#include <algorithm>
#include <concepts>
#include <print>
#include <vector>

// 半開区間 [left, right) のインデックスを表す構造体
struct IndexRange {
    int left;
    int right;

    // 条件を満たす要素の個数
    [[nodiscard]] int count() const noexcept { return right - left; }

    // 条件を満たす要素が存在するか
    [[nodiscard]] bool empty() const noexcept { return left >= right; }
};

// 内部実装用ヘルパー
namespace internal {
template <typename T, typename Val>
IndexRange get_range(const std::vector<T>& v, Val a, Val b, bool left_closed, bool right_closed) {
    auto it_l = left_closed ? std::ranges::lower_bound(v, a) : std::ranges::upper_bound(v, a);
    auto it_r = right_closed ? std::ranges::upper_bound(v, b) : std::ranges::lower_bound(v, b);

    int l = static_cast<int>(std::distance(v.begin(), it_l));
    int r = static_cast<int>(std::distance(v.begin(), it_r));

    if (l > r) r = l;  // a > b などの逆転防止
    return IndexRange{l, r};
}
}  // namespace internal

// 1. (a, b) : a < x < b
template <typename T, typename Val>
IndexRange range_open(const std::vector<T>& v, Val a, Val b) {
    return internal::get_range(v, a, b, false, false);
}

// 2. [a, b] : a <= x <= b
template <typename T, typename Val>
IndexRange range_closed(const std::vector<T>& v, Val a, Val b) {
    return internal::get_range(v, a, b, true, true);
}

// 3. [a, b) : a <= x < b
template <typename T, typename Val>
IndexRange range_left_closed(const std::vector<T>& v, Val a, Val b) {
    return internal::get_range(v, a, b, true, false);
}

// 4. (a, b] : a < x <= b
template <typename T, typename Val>
IndexRange range_right_closed(const std::vector<T>& v, Val a, Val b) {
    return internal::get_range(v, a, b, false, true);
}

// 短縮記法 (数学の区間表記 [a, b], (a, b), [a, b), (a, b] に対応)
template <typename T, typename Val>
IndexRange range_oo(const std::vector<T>& v, Val a, Val b) {
    return range_open(v, a, b);
}
template <typename T, typename Val>
IndexRange range_cc(const std::vector<T>& v, Val a, Val b) {
    return range_closed(v, a, b);
}
template <typename T, typename Val>
IndexRange range_co(const std::vector<T>& v, Val a, Val b) {
    return range_left_closed(v, a, b);
}
template <typename T, typename Val>
IndexRange range_oc(const std::vector<T>& v, Val a, Val b) {
    return range_right_closed(v, a, b);
}
