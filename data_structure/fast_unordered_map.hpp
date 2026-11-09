#include <algorithm>
#include <chrono>
#include <print>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// ------------------- Anti-Hack Custom Hash -------------------
struct FastHash {
    static long splitmix64(long x) {
        unsigned long z = static_cast<unsigned long>(x) + 0x9e3779b97f4a7c15ULL;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return static_cast<long>(z ^ (z >> 31));
    }

    template <typename T>
    long operator()(const T& x) const {
        static const long FIXED_RANDOM = std::chrono::steady_clock::now().time_since_epoch().count();
        if constexpr (std::is_integral_v<T>) {
            return splitmix64(static_cast<long>(x) + FIXED_RANDOM);
        } else {
            return splitmix64(std::hash<T>{}(x) + FIXED_RANDOM);
        }
    }
};

template <typename Key, typename Val>
struct UnorderedBucket {
    Key first{};
    Val second{};
    bool occupied = false;
};

// 構造化束縛用設定
namespace std {
template <typename K, typename V>
struct tuple_size<UnorderedBucket<K, V>> : std::integral_constant<size_t, 2> {};
template <typename K, typename V>
struct tuple_size<const UnorderedBucket<K, V>> : std::integral_constant<size_t, 2> {};
template <size_t I, typename K, typename V>
struct tuple_element<I, UnorderedBucket<K, V>> {
    using type = std::conditional_t<I == 0, const K, V>;
};
template <size_t I, typename K, typename V>
struct tuple_element<I, const UnorderedBucket<K, V>> {
    using type = std::conditional_t<I == 0, const K, const V>;
};
}  // namespace std
template <size_t I, typename K, typename V>
decltype(auto) get(UnorderedBucket<K, V>& b) {
    if constexpr (I == 0)
        return (b.first);
    else
        return (b.second);
}
template <size_t I, typename K, typename V>
decltype(auto) get(const UnorderedBucket<K, V>& b) {
    if constexpr (I == 0)
        return (b.first);
    else
        return (b.second);
}

// ------------------- FastUnorderedMap -------------------
template <typename Key, typename Val, typename Hash = FastHash>
class FastUnorderedMap {
   public:
    using Bucket = UnorderedBucket<Key, Val>;

   private:
    int cap;
    int mask;
    int num_entries;
    std::vector<Bucket> table;
    Hash hasher;

    void rehash(int new_cap) {
        std::vector<Bucket> old_table = std::move(table);
        cap = new_cap;
        mask = cap - 1;
        num_entries = 0;
        table.assign(cap, Bucket{});
        for (auto& b : old_table) {
            if (b.occupied) {
                int idx = hasher(b.first) & mask;
                while (table[idx].occupied)
                    idx = (idx + 1) & mask;
                table[idx].occupied = true;
                table[idx].first = std::move(b.first);
                table[idx].second = std::move(b.second);
                num_entries++;
            }
        }
    }

   public:
    FastUnorderedMap(int initial_cap = 8) : num_entries(0) {
        int c = 1;
        while (c < initial_cap)
            c <<= 1;
        cap = std::max(8, c);
        mask = cap - 1;
        table.resize(cap);
    }

    void reserve(int expected_size) {
        int target_cap = 8;
        while (target_cap < expected_size * 2)
            target_cap <<= 1;
        if (target_cap > cap) rehash(target_cap);
    }

    int size() const { return num_entries; }
    bool empty() const { return num_entries == 0; }
    void clear() {
        for (auto& b : table)
            b.occupied = false;
        num_entries = 0;
    }

    Val& operator[](const Key& key) {
        if (num_entries * 2 >= cap) rehash(cap * 2);
        int idx = hasher(key) & mask;
        while (table[idx].occupied) {
            if (table[idx].first == key) return table[idx].second;
            idx = (idx + 1) & mask;
        }
        table[idx].occupied = true;
        table[idx].first = key;
        table[idx].second = Val{};
        num_entries++;
        return table[idx].second;
    }

    bool contains(const Key& key) const {
        if (num_entries == 0) return false;
        int idx = hasher(key) & mask;
        while (table[idx].occupied) {
            if (table[idx].first == key) return true;
            idx = (idx + 1) & mask;
        }
        return false;
    }

    Val* find_ptr(const Key& key) {
        if (num_entries == 0) return nullptr;
        int idx = hasher(key) & mask;
        while (table[idx].occupied) {
            if (table[idx].first == key) return &table[idx].second;
            idx = (idx + 1) & mask;
        }
        return nullptr;
    }

    struct Iterator {
        Bucket* ptr;
        Bucket* end_ptr;
        void advance() {
            while (ptr != end_ptr && !ptr->occupied)
                ++ptr;
        }
        Iterator(Bucket* p, Bucket* end) : ptr(p), end_ptr(end) { advance(); }
        Bucket& operator*() const { return *ptr; }
        Bucket* operator->() const { return ptr; }
        Iterator& operator++() {
            ++ptr;
            advance();
            return *this;
        }
        bool operator==(const Iterator& o) const { return ptr == o.ptr; }
        bool operator!=(const Iterator& o) const { return ptr != o.ptr; }
    };

    Iterator begin() { return Iterator(table.data(), table.data() + cap); }
    Iterator end() { return Iterator(table.data() + cap, table.data() + cap); }
};
