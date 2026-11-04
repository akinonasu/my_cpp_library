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

namespace detail::index_comb_impl {

template <std::size_t K, std::integral Int>
    requires(K > 0)
class view : public std::ranges::view_interface<view<K, Int>> {
   private:
    Int first_{0};
    Int last_{0};

    class iterator {
       private:
        Int last_{0};
        std::array<Int, K> idx_{};
        bool is_end_ = false;

       public:
        using iterator_concept = std::forward_iterator_tag;
        using iterator_category = std::forward_iterator_tag;
        using value_type = decltype([]<std::size_t... Is>(std::index_sequence<Is...>) { return std::tuple<decltype(((void)Is, Int{}))...>{}; }(std::make_index_sequence<K>{}));
        using difference_type = std::ptrdiff_t;

        iterator() = default;

        constexpr iterator(Int first, Int last) : last_(last) {
            // 要素数が K 未満なら何もしない（即終了）
            if (last - first < static_cast<Int>(K)) {
                is_end_ = true;
                return;
            }
            for (std::size_t i = 0; i < K; ++i) {
                idx_[i] = first + static_cast<Int>(i);
            }
        }

        constexpr auto operator*() const noexcept {
            return [&]<std::size_t... Is>(std::index_sequence<Is...>) { return std::make_tuple(idx_[Is]...); }(std::make_index_sequence<K>{});
        }

        constexpr iterator& operator++() noexcept {
            for (std::size_t i = K; i-- > 0;) {
                // 上限チェック: i 番目が到達可能な最大値未満か？
                if (idx_[i] + 1 <= last_ - static_cast<Int>(K - i)) {
                    ++idx_[i];
                    for (std::size_t j = i + 1; j < K; ++j) {
                        idx_[j] = idx_[j - 1] + 1;
                    }
                    return *this;
                }
            }
            is_end_ = true;
            return *this;
        }

        constexpr iterator operator++(int) noexcept {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        constexpr bool operator==(std::default_sentinel_t) const noexcept { return is_end_; }
        constexpr bool operator==(const iterator& other) const noexcept {
            if (is_end_ && other.is_end_) return true;
            if (is_end_ != other.is_end_) return false;
            return idx_ == other.idx_;
        }
    };

   public:
    constexpr view() = default;
    constexpr explicit view(Int n) : first_(0), last_(n) {}
    constexpr view(Int first, Int last) : first_(first), last_(last) {}

    constexpr auto begin() const { return iterator(first_, last_); }
    constexpr auto end() const { return std::default_sentinel; }
};

}  // namespace detail::index_comb_impl

// =========================================================
// グローバル関数: どこからでも直接呼べる
// =========================================================

// 0 <= i < j < k < n
template <std::size_t K, std::integral Int>
[[nodiscard]] constexpr auto combinations(Int n) {
    return detail::index_comb_impl::view<K, Int>(n);
}

// first <= i < j < k < last
template <std::size_t K, std::integral Int>
[[nodiscard]] constexpr auto combinations(Int first, Int last) {
    return detail::index_comb_impl::view<K, Int>(first, last);
}

// 別名としても使えるようにエイリアス関数も提供
template <std::size_t K, std::integral Int>
[[nodiscard]] constexpr auto index_combinations(Int n) {
    return combinations<K, Int>(n);
}

template <std::size_t K, std::integral Int>
[[nodiscard]] constexpr auto index_combinations(Int first, Int last) {
    return combinations<K, Int>(first, last);
}
