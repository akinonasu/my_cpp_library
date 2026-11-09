#include <algorithm>
#include <print>
#include <utility>
#include <vector>

// ============================================================================
// Fast B-Tree Map (std::map の高速代替)
// ============================================================================
template <typename Key, typename Val, int B = 8>
class FastBTreeMap {
   public:
    static constexpr int MAX_K = 2 * B - 1;  // 1ノードあたりの最大キー数 (15)
    static constexpr int MAX_C = 2 * B;      // 最大子ノード数 (16)

    struct Entry {
        Key key{};
        Val val{};
    };

    struct Node {
        int num = 0;
        bool leaf = true;
        Entry entries[MAX_K];
        int children[MAX_C];
    };

   private:
    int root = -1;
    int num_entries = 0;
    std::vector<Node> pool;

    int new_node(bool leaf) {
        int idx = static_cast<int>(pool.size());
        pool.emplace_back();
        pool.back().leaf = leaf;
        return idx;
    }

    void split_child(int p_idx, int i) {
        int c_idx = pool[p_idx].children[i];
        int z_idx = new_node(pool[c_idx].leaf);

        auto& c = pool[c_idx];
        auto& z = pool[z_idx];

        z.num = B - 1;
        for (int j = 0; j < B - 1; ++j) {
            z.entries[j] = std::move(c.entries[j + B]);
        }
        if (!c.leaf) {
            for (int j = 0; j < B; ++j) {
                z.children[j] = c.children[j + B];
            }
        }
        c.num = B - 1;

        auto& p = pool[p_idx];
        for (int j = p.num; j > i; --j) {
            p.children[j + 1] = p.children[j];
            p.entries[j] = std::move(p.entries[j - 1]);
        }
        p.children[i + 1] = z_idx;
        p.entries[i] = std::move(c.entries[B - 1]);
        p.num++;
    }

    Val* insert_non_full(int u, const Key& key) {
        auto& node = pool[u];
        // ノード内を連続配列上で二分探索（L1キャッシュヒット）
        int i = std::lower_bound(node.entries, node.entries + node.num, key, [](const Entry& e, const Key& k) { return e.key < k; }) - node.entries;

        if (i < node.num && node.entries[i].key == key) {
            return &node.entries[i].val;
        }

        if (node.leaf) {
            for (int j = node.num; j > i; --j) {
                node.entries[j] = std::move(node.entries[j - 1]);
            }
            node.entries[i] = Entry{key, Val{}};
            node.num++;
            num_entries++;
            return &node.entries[i].val;
        } else {
            if (pool[node.children[i]].num == MAX_K) {
                split_child(u, i);
                if (key > pool[u].entries[i].key) {
                    i++;
                } else if (key == pool[u].entries[i].key) {
                    return &pool[u].entries[i].val;
                }
            }
            return insert_non_full(pool[u].children[i], key);
        }
    }

   public:
    FastBTreeMap() { pool.reserve(1024); }

    int size() const { return num_entries; }
    bool empty() const { return num_entries == 0; }
    void clear() {
        pool.clear();
        root = -1;
        num_entries = 0;
    }

    Val& operator[](const Key& key) {
        if (root == -1) {
            root = new_node(true);
        }
        if (pool[root].num == MAX_K) {
            int s = new_node(false);
            pool[s].children[0] = root;
            split_child(s, 0);
            root = s;
        }
        return *insert_non_full(root, key);
    }

    // 高速検索 (ポインタ取得)
    Val* find_ptr(const Key& key) {
        int u = root;
        while (u != -1) {
            auto& node = pool[u];
            int i = std::lower_bound(node.entries, node.entries + node.num, key, [](const Entry& e, const Key& k) { return e.key < k; }) - node.entries;
            if (i < node.num && node.entries[i].key == key) {
                return &node.entries[i].val;
            }
            if (node.leaf) break;
            u = node.children[i];
        }
        return nullptr;
    }

    bool contains(const Key& key) { return find_ptr(key) != nullptr; }

    // 二分探索: key 以上の最小の要素へのポインタ (Entry: {key, val}) を返す
    Entry* lower_bound(const Key& key) {
        Entry* best = nullptr;
        int u = root;
        while (u != -1) {
            auto& node = pool[u];
            int i = std::lower_bound(node.entries, node.entries + node.num, key, [](const Entry& e, const Key& k) { return e.key < k; }) - node.entries;
            if (i < node.num) {
                best = &node.entries[i];
                if (node.entries[i].key == key) return best;  // 完全一致
            }
            if (node.leaf) break;
            u = node.children[i];
        }
        return best;
    }

    // 昇順走査
    template <typename F>
    void for_each(F&& func) const {
        auto dfs = [&](auto self, int u) -> void {
            if (u == -1) return;
            const auto& node = pool[u];
            for (int i = 0; i < node.num; ++i) {
                if (!node.leaf) self(self, node.children[i]);
                func(node.entries[i].key, node.entries[i].val);
            }
            if (!node.leaf) self(self, node.children[node.num]);
        };
        dfs(dfs, root);
    }
};
