#include <algorithm>
#include <limits>
#include <optional>
#include <print>
#include <set>
#include <utility>
#include <vector>

/**
 * @brief 互いに素な半開区間 [l, r) を管理するデータ構造
 * @tparam T 座標の型 (int または long)
 */
template <typename T = long>
class IntervalSet {
   private:
    std::set<std::pair<T, T>> intervals;
    T total_length = 0;

   public:
    IntervalSet() = default;

    /**
     * @brief 半開区間 [l, r) を追加する (重複・隣接区間はマージされる)
     * @return 新たにカバーされた長さ (増分)
     * @complexity ならし O(log N)
     */
    T insert(T l, T r) {
        if (l >= r) return 0;
        auto it = intervals.lower_bound({l, l});
        if (it != intervals.begin()) {
            auto pit = std::prev(it);
            if (pit->second >= l) {  // 隣接または重複
                it = pit;
            }
        }

        T nl = l, nr = r;
        T removed_len = 0;
        while (it != intervals.end() && it->first <= r) {
            nl = std::min(nl, it->first);
            nr = std::max(nr, it->second);
            removed_len += it->second - it->first;
            it = intervals.erase(it);
        }

        intervals.emplace(nl, nr);
        T added_len = (nr - nl) - removed_len;
        total_length += added_len;
        return added_len;
    }

    /**
     * @brief 半開区間 [l, r) を削除する (重なっている部分を切り取る)
     * @return 実際に削除された長さ
     * @complexity ならし O(log N)
     */
    T erase(T l, T r) {
        if (l >= r) return 0;
        auto it = intervals.lower_bound({l, l});
        if (it != intervals.begin()) {
            auto pit = std::prev(it);
            if (pit->second > l) {  // 厳密に交差している場合のみ
                it = pit;
            }
        }

        std::vector<std::pair<T, T>> to_add;
        T erased_len = 0;
        while (it != intervals.end() && it->first < r) {
            auto [il, ir] = *it;
            it = intervals.erase(it);
            erased_len += ir - il;
            if (il < l) {
                to_add.emplace_back(il, l);
            }
            if (ir > r) {
                to_add.emplace_back(r, ir);
            }
        }

        for (const auto& iv : to_add) {
            intervals.insert(iv);
            erased_len -= (iv.second - iv.first);
        }

        total_length -= erased_len;
        return erased_len;
    }

    /**
     * @brief 点 x を含む区間 [l, r) を取得する
     */
    std::optional<std::pair<T, T>> find(T x) const {
        auto it = intervals.upper_bound({x, std::numeric_limits<T>::max()});
        if (it != intervals.begin()) {
            --it;
            if (it->first <= x && x < it->second) {
                return *it;
            }
        }
        return std::nullopt;
    }

    /**
     * @brief 点 x がカバーされているかを判定する
     */
    bool covered(T x) const { return find(x).has_value(); }

    bool contains(T x) const { return covered(x); }

    /**
     * @brief 区間 [l, r) 全体が単一の区間にカバーされているかを判定する
     */
    bool covered(T l, T r) const {
        if (l >= r) return true;
        auto it = intervals.upper_bound({l, std::numeric_limits<T>::max()});
        if (it != intervals.begin()) {
            --it;
            if (it->first <= l && r <= it->second) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief x 以上の最小のカバーされていない点 (mex) を取得する
     */
    T mex(T x = 0) const {
        auto it = intervals.upper_bound({x, std::numeric_limits<T>::max()});
        if (it != intervals.begin()) {
            --it;
            if (it->first <= x && x < it->second) {
                return it->second;
            }
        }
        return x;
    }

    /**
     * @brief カバーされている全体の長さを取得する
     */
    T length() const noexcept { return total_length; }

    /**
     * @brief 存在する区間の個数を取得する
     */
    int size() const noexcept { return static_cast<int>(intervals.size()); }

    bool empty() const noexcept { return intervals.empty(); }

    void clear() noexcept {
        intervals.clear();
        total_length = 0;
    }

    // イテレータのサポート (範囲for文で [l, r) を走査可能)
    auto begin() const noexcept { return intervals.begin(); }
    auto end() const noexcept { return intervals.end(); }
    auto rbegin() const noexcept { return intervals.rbegin(); }
    auto rend() const noexcept { return intervals.rend(); }
};
