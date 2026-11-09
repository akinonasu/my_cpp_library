#include <algorithm>
#include <chrono>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <print>
#include <type_traits>
#include <utility>
#include <vector>

// ============================================================================
// 1. Anti-Hack Custom Hash (SplitMix64)
// ============================================================================
struct FastHash {
    static long splitmix64(long x) {
        // 符号付き整数のオーバーフロー未定義動作(UB)を防ぐため unsigned long で演算
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

    // std::pair をキーにする場合のサポート（未使用だった FIXED_RANDOM を削除）
    template <typename T1, typename T2>
    long operator()(const std::pair<T1, T2>& p) const {
        long h1 = (*this)(p.first);
        long h2 = (*this)(p.second);
        return splitmix64(h1 ^ (h2 + static_cast<long>(0x9e3779b97f4a7c15ULL) + (h1 << 6) + (h1 >> 2)));
    }
};
// ============================================================================
// 2. Bucket Definition & Structured Binding Support
// ============================================================================
template <typename Key, typename Val>
struct FastMapBucket {
    Key first{};
    Val second{};
    bool occupied = false;
};

// 構造化束縛 (for (auto& [k, v] : map)) のための特殊化
namespace std {
template <typename K, typename V>
struct tuple_size<FastMapBucket<K, V>> : std::integral_constant<size_t, 2> {};

template <typename K, typename V>
struct tuple_size<const FastMapBucket<K, V>> : std::integral_constant<size_t, 2> {};

template <size_t I, typename K, typename V>
struct tuple_element<I, FastMapBucket<K, V>> {
    using type = std::conditional_t<I == 0, const K, V>;
};

template <size_t I, typename K, typename V>
struct tuple_element<I, const FastMapBucket<K, V>> {
    using type = std::conditional_t<I == 0, const K, const V>;
};
}  // namespace std

template <size_t I, typename K, typename V>
decltype(auto) get(FastMapBucket<K, V>& b) {
    if constexpr (I == 0)
        return (b.first);
    else
        return (b.second);
}

template <size_t I, typename K, typename V>
decltype(auto) get(const FastMapBucket<K, V>& b) {
    if constexpr (I == 0)
        return (b.first);
    else
        return (b.second);
}

// ============================================================================
// 3. Fast Flat Hash Map
// ============================================================================
template <typename Key, typename Val, typename Hash = FastHash>
class FastMap {
   public:
    using Bucket = FastMapBucket<Key, Val>;

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
                while (table[idx].occupied) {
                    idx = (idx + 1) & mask;
                }
                table[idx].occupied = true;
                table[idx].first = std::move(b.first);
                table[idx].second = std::move(b.second);
                num_entries++;
            }
        }
    }

   public:
    FastMap(int initial_cap = 8) : cap(initial_cap), num_entries(0) {
        int c = 1;
        while (c < initial_cap)
            c <<= 1;
        cap = std::max(8, c);
        mask = cap - 1;
        table.resize(cap);
    }

    // 要素数をあらかじめ予測して確保（負荷率0.5以下を維持する2の累乗サイズ）
    void reserve(int expected_size) {
        int target_cap = 8;
        while (target_cap < expected_size * 2) {
            target_cap <<= 1;
        }
        if (target_cap > cap) {
            rehash(target_cap);
        }
    }

    int size() const { return num_entries; }
    int capacity() const { return cap; }
    bool empty() const { return num_entries == 0; }

    void clear() {
        for (auto& b : table) {
            b.occupied = false;
        }
        num_entries = 0;
    }

    Val& operator[](const Key& key) {
        if (num_entries * 2 >= cap) {
            rehash(cap * 2);
        }
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

    Val& operator[](Key&& key) {
        if (num_entries * 2 >= cap) {
            rehash(cap * 2);
        }
        int idx = hasher(key) & mask;
        while (table[idx].occupied) {
            if (table[idx].first == key) return table[idx].second;
            idx = (idx + 1) & mask;
        }
        table[idx].occupied = true;
        table[idx].first = std::move(key);
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

    int count(const Key& key) const { return contains(key) ? 1 : 0; }

    Val* find_ptr(const Key& key) {
        if (num_entries == 0) return nullptr;
        int idx = hasher(key) & mask;
        while (table[idx].occupied) {
            if (table[idx].first == key) return &table[idx].second;
            idx = (idx + 1) & mask;
        }
        return nullptr;
    }

    const Val* find_ptr(const Key& key) const {
        if (num_entries == 0) return nullptr;
        int idx = hasher(key) & mask;
        while (table[idx].occupied) {
            if (table[idx].first == key) return &table[idx].second;
            idx = (idx + 1) & mask;
        }
        return nullptr;
    }

    // 後方シフト削除（Tombstoneを残さず連鎖クラスタを維持）
    bool erase(const Key& key) {
        if (num_entries == 0) return false;
        int idx = hasher(key) & mask;
        while (table[idx].occupied) {
            if (table[idx].first == key) {
                int curr = idx;
                int next = (curr + 1) & mask;
                while (table[next].occupied) {
                    int natural = hasher(table[next].first) & mask;
                    if (((next - natural) & mask) >= ((next - curr) & mask)) {
                        table[curr] = std::move(table[next]);
                        curr = next;
                    }
                    next = (next + 1) & mask;
                }
                table[curr].occupied = false;
                table[curr].first = Key{};
                table[curr].second = Val{};
                num_entries--;
                return true;
            }
            idx = (idx + 1) & mask;
        }
        return false;
    }

    // イテレータ実装
    template <typename B>
    struct IterBase {
        B* ptr;
        B* end_ptr;

        void advance() {
            while (ptr != end_ptr && !ptr->occupied) {
                ++ptr;
            }
        }

        IterBase(B* p, B* end) : ptr(p), end_ptr(end) { advance(); }

        B& operator*() const { return *ptr; }
        B* operator->() const { return ptr; }

        IterBase& operator++() {
            ++ptr;
            advance();
            return *this;
        }

        bool operator==(const IterBase& other) const { return ptr == other.ptr; }
        bool operator!=(const IterBase& other) const { return ptr != other.ptr; }
    };

    using Iterator = IterBase<Bucket>;
    using ConstIterator = IterBase<const Bucket>;

    Iterator begin() { return Iterator(table.data(), table.data() + cap); }
    Iterator end() { return Iterator(table.data() + cap, table.data() + cap); }
    ConstIterator begin() const { return ConstIterator(table.data(), table.data() + cap); }
    ConstIterator end() const { return ConstIterator(table.data() + cap, table.data() + cap); }
};
