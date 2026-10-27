#include <cstddef>
#include <iostream>
#include <stdexcept>

// GCC拡張 PBDS のヘッダ
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

template <typename T, typename Compare = std::less<T>>
class OrderedSet {
   private:
    // PBDSの内部型定義
    using pbds_tree =
        __gnu_pbds::tree<T, __gnu_pbds::null_type, Compare, __gnu_pbds::rb_tree_tag, __gnu_pbds::tree_order_statistics_node_update>;

    pbds_tree tree;

   public:
    OrderedSet() = default;

    // --- 基本的な操作 ---
    bool insert(const T& val) { return tree.insert(val).second; }

    bool erase(const T& val) { return tree.erase(val); }

    bool contains(const T& val) const { return tree.find(val) != tree.end(); }

    size_t size() const noexcept { return tree.size(); }

    bool empty() const noexcept { return tree.empty(); }

    void clear() noexcept { tree.clear(); }

    // --- PBDS特有の拡張機能 ---

    // 0-indexed で k 番目に小さい要素のイテレータを返す
    auto find_by_order(size_t k) const { return tree.find_by_order(k); }

    // 0-indexed で k 番目に小さい要素の値を直接取得する（範囲外は例外送出）
    const T& get_kth(size_t k) const {
        auto it = tree.find_by_order(k);
        if (it == tree.end()) {
            throw std::out_of_range("OrderedSet::get_kth index out of range");
        }
        return *it;
    }

    // val 未満の要素の個数を返す（val が集合内に存在しなくても動く）
    size_t order_of_key(const T& val) const { return tree.order_of_key(val); }

    // --- イテレータ ---
    auto begin() const noexcept { return tree.begin(); }
    auto end() const noexcept { return tree.end(); }
    auto rbegin() const noexcept { return tree.rbegin(); }
    auto rend() const noexcept { return tree.rend(); }
};
