#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <set>
#include <utility>

// 半開区間 [l, r) の集合を管理するデータ構造（隣接・交差する区間は自動でマージされる）
// 空間計算量: O(N) (N は保持している極大区間の数)
// 時間計算量: insert / erase は償却 O(log N)、その他クエリは O(log N) (size / total_length は O(1))
template <class T = long>
struct IntervalSet {
   private:
    std::set<std::pair<T, T>> st_;
    T total_len_ = 0;

    // first > x を満たす最初のイテレータを返す
    auto upper_bound_left(T x) const {
        auto it = st_.lower_bound({x, x});
        if (it != st_.end() && it->first == x) ++it;
        return it;
    }

   public:
    long size() const { return std::ssize(st_); }
    bool empty() const { return st_.empty(); }
    T total_length() const { return total_len_; }

    void clear() {
        st_.clear();
        total_len_ = 0;
    }

    // 半開区間 [l, r) を追加し、新たに増加した長さを返す
    // 事前条件: l <= r
    T insert(T l, T r) {
        assert(l <= r);
        if (l == r) return 0;

        auto it = upper_bound_left(l);
        if (it != st_.begin()) {
            auto prev = std::prev(it);
            if (prev->second >= r) return 0;
            if (prev->second >= l) it = prev;
        }

        T removed = 0;
        while (it != st_.end() && it->first <= r) {
            l = std::min(l, it->first);
            r = std::max(r, it->second);
            removed += it->second - it->first;
            it = st_.erase(it);
        }
        st_.emplace_hint(it, l, r);

        T added = (r - l) - removed;
        total_len_ += added;
        return added;
    }

    // 1点 [x, x + 1) を追加し、増加した長さ (0 または 1) を返す
    T insert(T x) { return insert(x, x + 1); }

    // 半開区間 [l, r) を削除し、減少した長さを返す
    // 事前条件: l <= r
    T erase(T l, T r) {
        assert(l <= r);
        if (l == r) return 0;

        auto it = upper_bound_left(l);
        if (it != st_.begin() && std::prev(it)->second > l) --it;

        T removed = 0;
        while (it != st_.end() && it->first < r) {
            auto [il, ir] = *it;
            it = st_.erase(it);
            total_len_ -= ir - il;
            removed += std::min(r, ir) - std::max(l, il);
            if (il < l) {
                st_.emplace_hint(it, il, l);
                total_len_ += l - il;
            }
            if (r < ir) {
                st_.emplace_hint(it, r, ir);
                total_len_ += ir - r;
            }
        }
        return removed;
    }

    // 1点 [x, x + 1) を削除し、減少した長さ (0 または 1) を返す
    T erase(T x) { return erase(x, x + 1); }

    // 点 x が集合に含まれるか
    bool contains(T x) const { return find(x).first; }

    // 半開区間 [l, r) 全体が集合に含まれるか
    // 事前条件: l <= r
    bool contains(T l, T r) const {
        assert(l <= r);
        if (l == r) return true;
        auto it = upper_bound_left(l);
        if (it == st_.begin()) return false;
        --it;
        return it->first <= l && r <= it->second;
    }

    // 半開区間 [l, r) と交わる部分が存在するか
    // 事前条件: l <= r
    bool intersects(T l, T r) const {
        assert(l <= r);
        if (l == r) return false;
        auto it = upper_bound_left(l);
        if (it != st_.begin() && std::prev(it)->second > l) return true;
        return it != st_.end() && it->first < r;
    }

    // 点 x を含む極大区間 [l, r) を返す（存在しない場合は {false, {0, 0}}）
    std::pair<bool, std::pair<T, T>> find(T x) const {
        auto it = upper_bound_left(x);
        if (it == st_.begin()) return {false, {T{}, T{}}};
        --it;
        if (it->second > x) return {true, *it};
        return {false, {T{}, T{}}};
    }

    // x 以上で集合に含まれない最小の値を返す
    T mex(T x = 0) const {
        auto it = upper_bound_left(x);
        if (it != st_.begin() && std::prev(it)->second > x) {
            return std::prev(it)->second;
        }
        return x;
    }

    const std::set<std::pair<T, T>>& intervals() const { return st_; }
    auto begin() const { return st_.begin(); }
    auto end() const { return st_.end(); }

    void dump() const {
#ifdef LOCAL
        std::cerr << "IntervalSet(len=" << total_len_ << "):";
        for (const auto& [l, r] : st_) {
            std::cerr << " [" << l << ", " << r << ")";
        }
        std::cerr << '\n';
#endif
    }
};
