#include <algorithm>
#include <chrono>
#include <print>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

template <typename Key, typename Val>
struct OrderedNode {
    Key first{};
    Val second{};
    int forward[20]{};  // レベルごとの次ノードのインデックス
};

// 構造化束縛用設定
namespace std {
template <typename K, typename V>
struct tuple_size<OrderedNode<K, V>> : std::integral_constant<size_t, 2> {};
template <typename K, typename V>
struct tuple_size<const OrderedNode<K, V>> : std::integral_constant<size_t, 2> {};
template <size_t I, typename K, typename V>
struct tuple_element<I, OrderedNode<K, V>> {
    using type = std::conditional_t<I == 0, const K, V>;
};
template <size_t I, typename K, typename V>
struct tuple_element<I, const OrderedNode<K, V>> {
    using type = std::conditional_t<I == 0, const K, const V>;
};
}  // namespace std
template <size_t I, typename K, typename V>
decltype(auto) get(OrderedNode<K, V>& b) {
    if constexpr (I == 0)
        return (b.first);
    else
        return (b.second);
}
template <size_t I, typename K, typename V>
decltype(auto) get(const OrderedNode<K, V>& b) {
    if constexpr (I == 0)
        return (b.first);
    else
        return (b.second);
}

// ------------------- FastOrderedMap (SkipList) -------------------
template <typename Key, typename Val>
class FastOrderedMap {
   public:
    using Node = OrderedNode<Key, Val>;
    static constexpr int MAX_LEVEL = 20;

   private:
    std::vector<Node> pool;
    int num_entries = 0;

    static unsigned long xor_shift() {
        static unsigned long x = 88172645463325252ULL;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        return x;
    }

    int random_level() {
        int lvl = 1;
        while ((xor_shift() & 1) && lvl < MAX_LEVEL)
            ++lvl;
        return lvl;
    }

   public:
    FastOrderedMap() {
        pool.reserve(1024);
        pool.emplace_back();  // index 0 をダミーヘッドノードとして使用
    }

    int size() const { return num_entries; }
    bool empty() const { return num_entries == 0; }
    void clear() {
        pool.clear();
        pool.emplace_back();
        num_entries = 0;
    }

    Val& operator[](const Key& key) {
        int cur = 0;
        int update[MAX_LEVEL];
        for (int i = MAX_LEVEL - 1; i >= 0; --i) {
            while (pool[cur].forward[i] != 0 && pool[pool[cur].forward[i]].first < key) {
                cur = pool[cur].forward[i];
            }
            update[i] = cur;
        }

        int next = pool[cur].forward[0];
        if (next != 0 && pool[next].first == key) {
            return pool[next].second;  // 既存キー
        }

        // 新規ノード挿入
        int lvl = random_level();
        int new_idx = static_cast<int>(pool.size());
        pool.emplace_back();
        pool[new_idx].first = key;
        pool[new_idx].second = Val{};

        for (int i = 0; i < lvl; ++i) {
            pool[new_idx].forward[i] = pool[update[i]].forward[i];
            pool[update[i]].forward[i] = new_idx;
        }
        num_entries++;
        return pool[new_idx].second;
    }

    bool contains(const Key& key) const {
        int cur = 0;
        for (int i = MAX_LEVEL - 1; i >= 0; --i) {
            while (pool[cur].forward[i] != 0 && pool[pool[cur].forward[i]].first < key) {
                cur = pool[cur].forward[i];
            }
        }
        int next = pool[cur].forward[0];
        return (next != 0 && pool[next].first == key);
    }

    // イテレータ (レベル0を辿ることで自動的にソート順になる)
    struct Iterator {
        int idx;
        std::vector<Node>* p_pool;

        Node& operator*() const { return (*p_pool)[idx]; }
        Node* operator->() const { return &(*p_pool)[idx]; }
        Iterator& operator++() {
            idx = (*p_pool)[idx].forward[0];
            return *this;
        }
        bool operator==(const Iterator& o) const { return idx == o.idx; }
        bool operator!=(const Iterator& o) const { return idx != o.idx; }
    };

    Iterator begin() { return Iterator(pool[0].forward[0], &pool); }
    Iterator end() { return Iterator(0, &pool); }

    // 二分探索: key 以上の最小要素を指すイテレータを返す
    Iterator lower_bound(const Key& key) {
        int cur = 0;
        for (int i = MAX_LEVEL - 1; i >= 0; --i) {
            while (pool[cur].forward[i] != 0 && pool[pool[cur].forward[i]].first < key) {
                cur = pool[cur].forward[i];
            }
        }
        return Iterator(pool[cur].forward[0], &pool);
    }
};
