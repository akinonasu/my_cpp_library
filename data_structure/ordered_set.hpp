#pragma once

#include <array>
#include <cassert>
#include <functional>
#include <utility>
#include <vector>

// High-speed order-statistic ordered set / multiset based on a split-node pooled B-Tree.
// Time complexity: O(log N) per operation (insert, erase, kth, order_of_key, bounds).
// Space complexity: O(N).
// Precondition: Compare defines a strict weak ordering on T; T is default-constructible and movable.
template <class T, class Compare = std::less<T>, bool Multi = false>
struct OrderedSet {
    explicit OrderedSet(Compare comp = Compare()) : comp_(comp) {}

    long size() const {
        if (root_ == -1) return 0;
        return root_is_leaf_ ? leaves_[root_].n : internals_[root_].sz;
    }

    bool empty() const { return size() == 0; }

    void clear() {
        leaves_.clear();
        internals_.clear();
        free_leaves_.clear();
        free_internals_.clear();
        root_ = -1;
        root_is_leaf_ = true;
    }

    // Reserves capacity for roughly n elements to reduce reallocations.
    void reserve(long n) {
        if (n > 0) {
            leaves_.reserve(static_cast<std::size_t>(n / T_MIN + 4));
            internals_.reserve(static_cast<std::size_t>(n / (T_MIN * T_MIN) + 4));
        }
    }

    bool insert(const T& x) {
        if (root_ == -1) {
            root_ = alloc_leaf();
            root_is_leaf_ = true;
            leaves_[root_].n = 1;
            leaves_[root_].keys[0] = x;
            return true;
        }
        int root_n = root_is_leaf_ ? leaves_[root_].n : internals_[root_].n;
        if (root_n == MAX_KEYS) {
            int s = alloc_internal();
            long old_sz = root_is_leaf_ ? leaves_[root_].n : internals_[root_].sz;
            internals_[s].n = 0;
            internals_[s].sz = old_sz;
            internals_[s].child_leaf = root_is_leaf_;
            internals_[s].ch[0] = root_;
            internals_[s].ch_sz[0] = old_sz;
            root_ = s;
            root_is_leaf_ = false;
            split_child(s, 0);
        }

        PathStep path[32];
        int depth = 0;
        int u = root_;
        bool is_leaf = root_is_leaf_;

        while (true) {
            if (is_leaf) {
                int i = lb_keys(leaves_[u].keys.data(), leaves_[u].n, x);
                if constexpr (!Multi) {
                    if (i < leaves_[u].n && eq(leaves_[u].keys[i], x)) {
                        return false;
                    }
                }
                for (int k = leaves_[u].n; k > i; --k) {
                    leaves_[u].keys[k] = std::move(leaves_[u].keys[k - 1]);
                }
                leaves_[u].keys[i] = x;
                leaves_[u].n++;
                break;
            }

            int i = lb_keys(internals_[u].keys.data(), internals_[u].n, x);
            if constexpr (!Multi) {
                if (i < internals_[u].n && eq(internals_[u].keys[i], x)) {
                    return false;
                }
            }
            int c = internals_[u].ch[i];
            bool c_leaf = internals_[u].child_leaf;
            int c_n = c_leaf ? leaves_[c].n : internals_[c].n;
            if (c_n == MAX_KEYS) {
                split_child(u, i);
                if constexpr (!Multi) {
                    if (eq(internals_[u].keys[i], x)) {
                        return false;
                    }
                }
                if (comp_(internals_[u].keys[i], x)) {
                    ++i;
                }
            }
            path[depth++] = {u, i};
            is_leaf = internals_[u].child_leaf;
            u = internals_[u].ch[i];
        }

        for (int d = 0; d < depth; ++d) {
            internals_[path[d].u].sz++;
            internals_[path[d].u].ch_sz[path[d].slot]++;
        }
        return true;
    }

    bool erase(const T& x) {
        if (empty()) return false;

        PathStep path[32];
        int depth = 0;
        int u = root_;
        bool is_leaf = root_is_leaf_;

        while (true) {
            if (is_leaf) {
                int i = lb_keys(leaves_[u].keys.data(), leaves_[u].n, x);
                if (i < leaves_[u].n && eq(leaves_[u].keys[i], x)) {
                    remove_leaf_key(u, i);
                    if (u == root_ && leaves_[u].n == 0) {
                        free_leaf(u);
                        root_ = -1;
                    }
                    break;
                }
                return false;
            }

            int i = lb_keys(internals_[u].keys.data(), internals_[u].n, x);
            if (i < internals_[u].n && eq(internals_[u].keys[i], x)) {
                int left_c = internals_[u].ch[i];
                int right_c = internals_[u].ch[i + 1];
                bool c_leaf = internals_[u].child_leaf;
                int left_n = c_leaf ? leaves_[left_c].n : internals_[left_c].n;
                int right_n = c_leaf ? leaves_[right_c].n : internals_[right_c].n;

                if (left_n >= T_MIN) {
                    path[depth++] = {u, i};
                    T pred = erase_max_down(left_c, c_leaf, path, depth);
                    internals_[u].keys[i] = std::move(pred);
                    break;
                }
                if (right_n >= T_MIN) {
                    path[depth++] = {u, i + 1};
                    T succ = erase_min_down(right_c, c_leaf, path, depth);
                    internals_[u].keys[i] = std::move(succ);
                    break;
                }
                merge_children(u, i);
                if (u == root_ && internals_[u].n == 0) {
                    root_ = left_c;
                    root_is_leaf_ = c_leaf;
                    free_internal(u);
                    u = left_c;
                    is_leaf = c_leaf;
                } else {
                    path[depth++] = {u, i};
                    u = left_c;
                    is_leaf = c_leaf;
                }
                continue;
            }

            int c = internals_[u].ch[i];
            bool c_leaf = internals_[u].child_leaf;
            int c_n = c_leaf ? leaves_[c].n : internals_[c].n;
            int slot = i;
            if (c_n == T_MIN - 1) {
                auto [next_c, next_slot] = fix_short_child(u, i);
                if (u == root_ && internals_[u].n == 0) {
                    root_ = next_c;
                    root_is_leaf_ = c_leaf;
                    free_internal(u);
                    u = next_c;
                    is_leaf = c_leaf;
                    continue;
                }
                c = next_c;
                slot = next_slot;
            }
            path[depth++] = {u, slot};
            u = c;
            is_leaf = c_leaf;
        }

        for (int d = 0; d < depth; ++d) {
            internals_[path[d].u].sz--;
            internals_[path[d].u].ch_sz[path[d].slot]--;
        }
        return true;
    }

    bool contains(const T& x) const {
        if (root_ == -1) return false;
        int u = root_;
        bool is_leaf = root_is_leaf_;
        while (!is_leaf) {
            const Internal& node = internals_[u];
            int i = lb_keys(node.keys.data(), node.n, x);
            if (i < node.n && eq(node.keys[i], x)) return true;
            is_leaf = node.child_leaf;
            u = node.ch[i];
        }
        const Leaf& leaf = leaves_[u];
        int i = lb_keys(leaf.keys.data(), leaf.n, x);
        return i < leaf.n && eq(leaf.keys[i], x);
    }

    long count(const T& x) const {
        if constexpr (!Multi) {
            return contains(x) ? 1 : 0;
        } else {
            return count_le(x) - order_of_key(x);
        }
    }

    // Returns the k-th smallest element (0-indexed).
    // Precondition: 0 <= k && k < size().
    T kth(long k) const {
        assert(0 <= k && k < size());
        int u = root_;
        bool is_leaf = root_is_leaf_;
        while (!is_leaf) {
            const Internal& node = internals_[u];
            for (int i = 0; i <= node.n; ++i) {
                long c_sz = node.ch_sz[i];
                if (k < c_sz) {
                    is_leaf = node.child_leaf;
                    u = node.ch[i];
                    break;
                }
                k -= c_sz;
                if (i < node.n) {
                    if (k == 0) return node.keys[i];
                    --k;
                }
            }
        }
        return leaves_[u].keys[k];
    }

    T operator[](long k) const { return kth(k); }

    // Returns the number of elements strictly less than x.
    long order_of_key(const T& x) const {
        if (root_ == -1) return 0;
        long res = 0;
        int u = root_;
        bool is_leaf = root_is_leaf_;
        while (!is_leaf) {
            const Internal& node = internals_[u];
            int i = lb_keys(node.keys.data(), node.n, x);
            res += i;
            for (int k = 0; k < i; ++k) {
                res += node.ch_sz[k];
            }
            is_leaf = node.child_leaf;
            u = node.ch[i];
        }
        res += lb_keys(leaves_[u].keys.data(), leaves_[u].n, x);
        return res;
    }

    // Returns the number of elements less than or equal to x.
    long count_le(const T& x) const {
        if (root_ == -1) return 0;
        long res = 0;
        int u = root_;
        bool is_leaf = root_is_leaf_;
        while (!is_leaf) {
            const Internal& node = internals_[u];
            int i = ub_keys(node.keys.data(), node.n, x);
            res += i;
            for (int k = 0; k < i; ++k) {
                res += node.ch_sz[k];
            }
            is_leaf = node.child_leaf;
            u = node.ch[i];
        }
        res += ub_keys(leaves_[u].keys.data(), leaves_[u].n, x);
        return res;
    }

    // Returns {found, min element >= x}.
    std::pair<bool, T> lower_bound(const T& x) const {
        if (root_ == -1) return {false, T{}};
        bool found = false;
        T res{};
        int u = root_;
        bool is_leaf = root_is_leaf_;
        while (!is_leaf) {
            const Internal& node = internals_[u];
            int i = lb_keys(node.keys.data(), node.n, x);
            if (i < node.n) {
                res = node.keys[i];
                found = true;
            }
            is_leaf = node.child_leaf;
            u = node.ch[i];
        }
        const Leaf& leaf = leaves_[u];
        int i = lb_keys(leaf.keys.data(), leaf.n, x);
        if (i < leaf.n) {
            res = leaf.keys[i];
            found = true;
        }
        return {found, res};
    }

    // Returns {found, min element > x}.
    std::pair<bool, T> upper_bound(const T& x) const {
        if (root_ == -1) return {false, T{}};
        bool found = false;
        T res{};
        int u = root_;
        bool is_leaf = root_is_leaf_;
        while (!is_leaf) {
            const Internal& node = internals_[u];
            int i = ub_keys(node.keys.data(), node.n, x);
            if (i < node.n) {
                res = node.keys[i];
                found = true;
            }
            is_leaf = node.child_leaf;
            u = node.ch[i];
        }
        const Leaf& leaf = leaves_[u];
        int i = ub_keys(leaf.keys.data(), leaf.n, x);
        if (i < leaf.n) {
            res = leaf.keys[i];
            found = true;
        }
        return {found, res};
    }

    // Returns {found, max element <= x}.
    std::pair<bool, T> prev_le(const T& x) const {
        if (root_ == -1) return {false, T{}};
        bool found = false;
        T res{};
        int u = root_;
        bool is_leaf = root_is_leaf_;
        while (!is_leaf) {
            const Internal& node = internals_[u];
            int i = ub_keys(node.keys.data(), node.n, x);
            if (i > 0) {
                res = node.keys[i - 1];
                found = true;
            }
            is_leaf = node.child_leaf;
            u = node.ch[i];
        }
        const Leaf& leaf = leaves_[u];
        int i = ub_keys(leaf.keys.data(), leaf.n, x);
        if (i > 0) {
            res = leaf.keys[i - 1];
            found = true;
        }
        return {found, res};
    }

    // Returns {found, max element < x}.
    std::pair<bool, T> prev_lt(const T& x) const {
        if (root_ == -1) return {false, T{}};
        bool found = false;
        T res{};
        int u = root_;
        bool is_leaf = root_is_leaf_;
        while (!is_leaf) {
            const Internal& node = internals_[u];
            int i = lb_keys(node.keys.data(), node.n, x);
            if (i > 0) {
                res = node.keys[i - 1];
                found = true;
            }
            is_leaf = node.child_leaf;
            u = node.ch[i];
        }
        const Leaf& leaf = leaves_[u];
        int i = lb_keys(leaf.keys.data(), leaf.n, x);
        if (i > 0) {
            res = leaf.keys[i - 1];
            found = true;
        }
        return {found, res};
    }

    std::vector<T> to_vector() const {
        std::vector<T> res;
        res.reserve(static_cast<std::size_t>(size()));
        if (root_ != -1) dfs(root_, root_is_leaf_, res);
        return res;
    }

   private:
    static constexpr int T_MIN = 16;
    static constexpr int MAX_KEYS = 2 * T_MIN - 1;
    static constexpr int MAX_CH = 2 * T_MIN;

    struct Leaf {
        int n = 0;
        std::array<T, MAX_KEYS> keys{};
    };

    struct Internal {
        long sz = 0;
        int n = 0;
        bool child_leaf = false;
        std::array<T, MAX_KEYS> keys{};
        std::array<int, MAX_CH> ch{};
        std::array<long, MAX_CH> ch_sz{};
    };

    struct PathStep {
        int u;
        int slot;
    };

    std::vector<Leaf> leaves_;
    std::vector<Internal> internals_;
    std::vector<int> free_leaves_;
    std::vector<int> free_internals_;
    int root_ = -1;
    bool root_is_leaf_ = true;
    Compare comp_;

    bool eq(const T& a, const T& b) const { return !comp_(a, b) && !comp_(b, a); }

    int alloc_leaf() {
        if (!free_leaves_.empty()) {
            int id = free_leaves_.back();
            free_leaves_.pop_back();
            leaves_[id].n = 0;
            return id;
        }
        leaves_.emplace_back();
        return static_cast<int>(leaves_.size()) - 1;
    }

    void free_leaf(int id) { free_leaves_.push_back(id); }

    int alloc_internal() {
        if (!free_internals_.empty()) {
            int id = free_internals_.back();
            free_internals_.pop_back();
            internals_[id].n = 0;
            internals_[id].sz = 0;
            return id;
        }
        internals_.emplace_back();
        return static_cast<int>(internals_.size()) - 1;
    }

    void free_internal(int id) { free_internals_.push_back(id); }

    int lb_keys(const T* keys, int n, const T& x) const {
        int len = n;
        int base = 0;
        while (len > 1) {
            int half = len >> 1;
            base += comp_(keys[base + half - 1], x) ? half : 0;
            len -= half;
        }
        if (len == 1 && comp_(keys[base], x)) ++base;
        return base;
    }

    int ub_keys(const T* keys, int n, const T& x) const {
        int len = n;
        int base = 0;
        while (len > 1) {
            int half = len >> 1;
            base += !comp_(x, keys[base + half - 1]) ? half : 0;
            len -= half;
        }
        if (len == 1 && !comp_(x, keys[base])) ++base;
        return base;
    }

    void remove_leaf_key(int u, int idx) {
        for (int k = idx; k + 1 < leaves_[u].n; ++k) {
            leaves_[u].keys[k] = std::move(leaves_[u].keys[k + 1]);
        }
        leaves_[u].n--;
    }

    void split_child(int u, int i) {
        if (internals_[u].child_leaf) {
            int y = internals_[u].ch[i];
            int z = alloc_leaf();

            leaves_[z].n = T_MIN - 1;
            for (int k = 0; k < T_MIN - 1; ++k) {
                leaves_[z].keys[k] = std::move(leaves_[y].keys[T_MIN + k]);
            }
            leaves_[y].n = T_MIN - 1;

            for (int k = internals_[u].n; k >= i + 1; --k) {
                internals_[u].ch[k + 1] = internals_[u].ch[k];
                internals_[u].ch_sz[k + 1] = internals_[u].ch_sz[k];
            }
            internals_[u].ch[i + 1] = z;
            internals_[u].ch_sz[i] = T_MIN - 1;
            internals_[u].ch_sz[i + 1] = T_MIN - 1;

            for (int k = internals_[u].n - 1; k >= i; --k) {
                internals_[u].keys[k + 1] = std::move(internals_[u].keys[k]);
            }
            internals_[u].keys[i] = std::move(leaves_[y].keys[T_MIN - 1]);
            internals_[u].n++;
        } else {
            int y = internals_[u].ch[i];
            int z = alloc_internal();

            internals_[z].child_leaf = internals_[y].child_leaf;
            internals_[z].n = T_MIN - 1;
            for (int k = 0; k < T_MIN - 1; ++k) {
                internals_[z].keys[k] = std::move(internals_[y].keys[T_MIN + k]);
            }
            long z_sz = T_MIN - 1;
            for (int k = 0; k < T_MIN; ++k) {
                internals_[z].ch[k] = internals_[y].ch[T_MIN + k];
                long c_sz = internals_[y].ch_sz[T_MIN + k];
                internals_[z].ch_sz[k] = c_sz;
                z_sz += c_sz;
            }
            internals_[z].sz = z_sz;
            internals_[y].n = T_MIN - 1;
            internals_[y].sz -= (z_sz + 1);
            long y_sz = internals_[y].sz;

            for (int k = internals_[u].n; k >= i + 1; --k) {
                internals_[u].ch[k + 1] = internals_[u].ch[k];
                internals_[u].ch_sz[k + 1] = internals_[u].ch_sz[k];
            }
            internals_[u].ch[i + 1] = z;
            internals_[u].ch_sz[i] = y_sz;
            internals_[u].ch_sz[i + 1] = z_sz;

            for (int k = internals_[u].n - 1; k >= i; --k) {
                internals_[u].keys[k + 1] = std::move(internals_[u].keys[k]);
            }
            internals_[u].keys[i] = std::move(internals_[y].keys[T_MIN - 1]);
            internals_[u].n++;
        }
    }

    std::pair<int, int> fix_short_child(int u, int i) {
        bool c_leaf = internals_[u].child_leaf;
        if (i > 0) {
            int left_sib = internals_[u].ch[i - 1];
            int left_n = c_leaf ? leaves_[left_sib].n : internals_[left_sib].n;
            if (left_n >= T_MIN) {
                borrow_from_left(u, i);
                return {internals_[u].ch[i], i};
            }
        }
        if (i < internals_[u].n) {
            int right_sib = internals_[u].ch[i + 1];
            int right_n = c_leaf ? leaves_[right_sib].n : internals_[right_sib].n;
            if (right_n >= T_MIN) {
                borrow_from_right(u, i);
                return {internals_[u].ch[i], i};
            }
        }
        if (i < internals_[u].n) {
            merge_children(u, i);
            return {internals_[u].ch[i], i};
        }
        assert(i > 0);
        if (i > 0) {
            merge_children(u, i - 1);
            return {internals_[u].ch[i - 1], i - 1};
        }
        std::unreachable();
    }

    void borrow_from_left(int u, int i) {
        int c = internals_[u].ch[i];
        int sib = internals_[u].ch[i - 1];

        if (internals_[u].child_leaf) {
            for (int k = leaves_[c].n; k >= 1; --k) {
                leaves_[c].keys[k] = std::move(leaves_[c].keys[k - 1]);
            }
            leaves_[c].keys[0] = std::move(internals_[u].keys[i - 1]);
            internals_[u].keys[i - 1] = std::move(leaves_[sib].keys[leaves_[sib].n - 1]);
            leaves_[c].n++;
            leaves_[sib].n--;
            internals_[u].ch_sz[i]++;
            internals_[u].ch_sz[i - 1]--;
        } else {
            for (int k = internals_[c].n; k >= 1; --k) {
                internals_[c].keys[k] = std::move(internals_[c].keys[k - 1]);
            }
            for (int k = internals_[c].n + 1; k >= 1; --k) {
                internals_[c].ch[k] = internals_[c].ch[k - 1];
                internals_[c].ch_sz[k] = internals_[c].ch_sz[k - 1];
            }
            internals_[c].keys[0] = std::move(internals_[u].keys[i - 1]);
            internals_[u].keys[i - 1] = std::move(internals_[sib].keys[internals_[sib].n - 1]);

            int moved_ch = internals_[sib].ch[internals_[sib].n];
            long moved_ch_sz = internals_[sib].ch_sz[internals_[sib].n];
            internals_[c].ch[0] = moved_ch;
            internals_[c].ch_sz[0] = moved_ch_sz;

            long moved_sz = 1 + moved_ch_sz;
            internals_[c].n++;
            internals_[sib].n--;
            internals_[c].sz += moved_sz;
            internals_[sib].sz -= moved_sz;
            internals_[u].ch_sz[i] += moved_sz;
            internals_[u].ch_sz[i - 1] -= moved_sz;
        }
    }

    void borrow_from_right(int u, int i) {
        int c = internals_[u].ch[i];
        int sib = internals_[u].ch[i + 1];

        if (internals_[u].child_leaf) {
            leaves_[c].keys[leaves_[c].n] = std::move(internals_[u].keys[i]);
            internals_[u].keys[i] = std::move(leaves_[sib].keys[0]);
            for (int k = 0; k + 1 < leaves_[sib].n; ++k) {
                leaves_[sib].keys[k] = std::move(leaves_[sib].keys[k + 1]);
            }
            leaves_[c].n++;
            leaves_[sib].n--;
            internals_[u].ch_sz[i]++;
            internals_[u].ch_sz[i + 1]--;
        } else {
            internals_[c].keys[internals_[c].n] = std::move(internals_[u].keys[i]);
            internals_[u].keys[i] = std::move(internals_[sib].keys[0]);

            int moved_ch = internals_[sib].ch[0];
            long moved_ch_sz = internals_[sib].ch_sz[0];
            internals_[c].ch[internals_[c].n + 1] = moved_ch;
            internals_[c].ch_sz[internals_[c].n + 1] = moved_ch_sz;

            for (int k = 0; k + 1 < internals_[sib].n; ++k) {
                internals_[sib].keys[k] = std::move(internals_[sib].keys[k + 1]);
            }
            for (int k = 0; k < internals_[sib].n; ++k) {
                internals_[sib].ch[k] = internals_[sib].ch[k + 1];
                internals_[sib].ch_sz[k] = internals_[sib].ch_sz[k + 1];
            }

            long moved_sz = 1 + moved_ch_sz;
            internals_[c].n++;
            internals_[sib].n--;
            internals_[c].sz += moved_sz;
            internals_[sib].sz -= moved_sz;
            internals_[u].ch_sz[i] += moved_sz;
            internals_[u].ch_sz[i + 1] -= moved_sz;
        }
    }

    void merge_children(int u, int i) {
        int left = internals_[u].ch[i];
        int right = internals_[u].ch[i + 1];

        if (internals_[u].child_leaf) {
            leaves_[left].keys[T_MIN - 1] = std::move(internals_[u].keys[i]);
            for (int k = 0; k < T_MIN - 1; ++k) {
                leaves_[left].keys[T_MIN + k] = std::move(leaves_[right].keys[k]);
            }
            leaves_[left].n = MAX_KEYS;
            internals_[u].ch_sz[i] = MAX_KEYS;
            free_leaf(right);
        } else {
            internals_[left].keys[T_MIN - 1] = std::move(internals_[u].keys[i]);
            for (int k = 0; k < T_MIN - 1; ++k) {
                internals_[left].keys[T_MIN + k] = std::move(internals_[right].keys[k]);
            }
            for (int k = 0; k < T_MIN; ++k) {
                internals_[left].ch[T_MIN + k] = internals_[right].ch[k];
                internals_[left].ch_sz[T_MIN + k] = internals_[right].ch_sz[k];
            }
            internals_[left].n = MAX_KEYS;
            internals_[left].sz += 1 + internals_[right].sz;
            internals_[u].ch_sz[i] = internals_[left].sz;
            free_internal(right);
        }

        for (int k = i; k + 1 < internals_[u].n; ++k) {
            internals_[u].keys[k] = std::move(internals_[u].keys[k + 1]);
        }
        for (int k = i + 1; k < internals_[u].n; ++k) {
            internals_[u].ch[k] = internals_[u].ch[k + 1];
            internals_[u].ch_sz[k] = internals_[u].ch_sz[k + 1];
        }
        internals_[u].n--;
    }

    T erase_max_down(int cur, bool is_leaf, PathStep* path, int& depth) {
        while (!is_leaf) {
            int last_idx = internals_[cur].n;
            int next_c = internals_[cur].ch[last_idx];
            bool next_leaf = internals_[cur].child_leaf;
            int next_n = next_leaf ? leaves_[next_c].n : internals_[next_c].n;
            int slot = last_idx;
            if (next_n == T_MIN - 1) {
                auto [fixed_c, fixed_slot] = fix_short_child(cur, last_idx);
                next_c = fixed_c;
                slot = fixed_slot;
            }
            path[depth++] = {cur, slot};
            cur = next_c;
            is_leaf = next_leaf;
        }
        T pred = std::move(leaves_[cur].keys[leaves_[cur].n - 1]);
        leaves_[cur].n--;
        return pred;
    }

    T erase_min_down(int cur, bool is_leaf, PathStep* path, int& depth) {
        while (!is_leaf) {
            int next_c = internals_[cur].ch[0];
            bool next_leaf = internals_[cur].child_leaf;
            int next_n = next_leaf ? leaves_[next_c].n : internals_[next_c].n;
            int slot = 0;
            if (next_n == T_MIN - 1) {
                auto [fixed_c, fixed_slot] = fix_short_child(cur, 0);
                next_c = fixed_c;
                slot = fixed_slot;
            }
            path[depth++] = {cur, slot};
            cur = next_c;
            is_leaf = next_leaf;
        }
        T succ = std::move(leaves_[cur].keys[0]);
        remove_leaf_key(cur, 0);
        return succ;
    }

    void dfs(int u, bool is_leaf, std::vector<T>& res) const {
        if (is_leaf) {
            for (int i = 0; i < leaves_[u].n; ++i) {
                res.push_back(leaves_[u].keys[i]);
            }
            return;
        }
        const Internal& node = internals_[u];
        for (int i = 0; i < node.n; ++i) {
            dfs(node.ch[i], node.child_leaf, res);
            res.push_back(node.keys[i]);
        }
        dfs(node.ch[node.n], node.child_leaf, res);
    }
};

template <class T, class Compare = std::less<T>>
using OrderedMultiset = OrderedSet<T, Compare, true>;
