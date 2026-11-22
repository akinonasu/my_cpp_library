#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <utility>

// GCC拡張 PBDS のヘッダ
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

// ==========================================
// 1. OrderedSet (重複なし)
// ==========================================
template <typename T, typename Compare = std::less<T>>
class OrderedSet {
   private:
    using pbds_tree =
        __gnu_pbds::tree<T, __gnu_pbds::null_type, Compare, __gnu_pbds::rb_tree_tag, __gnu_pbds::tree_order_statistics_node_update>;
    pbds_tree tree;

   public:
    OrderedSet() = default;

    bool insert(const T& val) { return tree.insert(val).second; }
    bool erase(const T& val) { return tree.erase(val); }
    bool contains(const T& val) const { return tree.find(val) != tree.end(); }
    size_t size() const noexcept { return tree.size(); }
    bool empty() const noexcept { return tree.empty(); }
    void clear() noexcept { tree.clear(); }

    // 0-indexed で k 番目に小さい要素の値を取得
    const T& get_kth(size_t k) const {
        auto it = tree.find_by_order(k);
        if (it == tree.end()) throw std::out_of_range("OrderedSet: out of range");
        return *it;
    }

    // val 未満の要素の個数を取得
    size_t order_of_key(const T& val) const { return tree.order_of_key(val); }

    auto begin() const noexcept { return tree.begin(); }
    auto end() const noexcept { return tree.end(); }
};

// ==========================================
// 2. OrderedMultiset (重複あり・安全版)
// ==========================================
template <typename T>
class OrderedMultiset {
   private:
    // 要素を一意にするため (値, ID) で管理
    using pbds_tree = __gnu_pbds::tree<std::pair<T, uint64_t>, __gnu_pbds::null_type, std::less<std::pair<T, uint64_t>>,
                                       __gnu_pbds::rb_tree_tag, __gnu_pbds::tree_order_statistics_node_update>;

    pbds_tree tree;
    uint64_t uid = 0;  // 追加順のユニークID

   public:
    OrderedMultiset() = default;

    // 要素を追加 (O(log N))
    void insert(const T& val) { tree.insert({val, uid++}); }

    // 要素を1つだけ削除 (存在すれば true, しなければ false) (O(log N))
    bool erase_one(const T& val) {
        auto it = tree.lower_bound({val, 0});
        if (it != tree.end() && it->first == val) {
            tree.erase(it);
            return true;
        }
        return false;
    }

    // val 未満の要素の個数を取得 (O(log N))
    size_t order_of_key(const T& val) const {
        // IDが 0 以上の要素のみ入っているため、(val, 0) 未満を数えれば良い
        return tree.order_of_key({val, 0});
    }

    // 0-indexed で k 番目に小さい要素の値を取得 (O(log N))
    const T& get_kth(size_t k) const {
        auto it = tree.find_by_order(k);
        if (it == tree.end()) throw std::out_of_range("OrderedMultiset: out of range");
        return it->first;
    }

    // 指定した値の個数を取得 (O(log N))
    size_t count(const T& val) const {
        // (val 以下の個数) - (val 未満の個数)
        return tree.order_of_key({val, UINT64_MAX}) - tree.order_of_key({val, 0});
    }

    bool contains(const T& val) const { return count(val) > 0; }
    size_t size() const noexcept { return tree.size(); }
    bool empty() const noexcept { return tree.empty(); }
    void clear() noexcept {
        tree.clear();
        uid = 0;
    }
};
