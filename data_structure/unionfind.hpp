#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

// 基本的な Union-Find (経路圧縮 + Union by Size)
// 空間計算量: O(N)
// 時間計算量: merge, leader, same, size は償却 O(α(N))、count は O(1)、groups は O(N α(N))
// 事前条件: 頂点番号は 0-indexed で [0, N)
struct UnionFind {
    UnionFind() = default;
    explicit UnionFind(long n) : data_(n, -1), comp_(n) { assert(n >= 0); }

    long leader(long u) {
        assert(0 <= u && u < std::ssize(data_));
        if (data_[u] < 0) return u;
        return data_[u] = leader(data_[u]);
    }

    bool same(long u, long v) { return leader(u) == leader(v); }

    // すでに同じ連結成分なら false、新しく併合したなら true を返す
    bool merge(long u, long v) {
        u = leader(u);
        v = leader(v);
        if (u == v) return false;
        if (-data_[u] < -data_[v]) std::swap(u, v);
        data_[u] += data_[v];
        data_[v] = u;
        --comp_;
        return true;
    }

    long size(long u) { return -data_[leader(u)]; }

    long count() const { return comp_; }

    std::vector<std::vector<long>> groups() {
        long n = std::ssize(data_);
        std::vector<long> leader_buf(n), group_size(n);
        for (long i = 0; i < n; ++i) {
            leader_buf[i] = leader(i);
            ++group_size[leader_buf[i]];
        }
        std::vector<std::vector<long>> res(n);
        for (long i = 0; i < n; ++i)
            res[i].reserve(group_size[i]);
        for (long i = 0; i < n; ++i)
            res[leader_buf[i]].push_back(i);
        std::erase_if(res, [](const auto& v) { return v.empty(); });
        return res;
    }

   private:
    std::vector<long> data_;
    long comp_ = 0;
};

// ポテンシャル付き Union-Find (重み付き Union-Find)
// 空間計算量: O(N)
// 時間計算量: merge, leader, same, weight, diff, size は償却 O(α(N))
// 事前条件: 型 T は +, -, +=, -=, ==, 単項 - および T{} (単位元) を持つアーベル群であること
template <class T = long>
struct WeightedUnionFind {
    WeightedUnionFind() = default;
    explicit WeightedUnionFind(long n) : data_(n, -1), diff_weight_(n, T{}), comp_(n) { assert(n >= 0); }

    long leader(long u) {
        assert(0 <= u && u < std::ssize(data_));
        if (data_[u] < 0) return u;
        long r = leader(data_[u]);
        diff_weight_[u] += diff_weight_[data_[u]];
        return data_[u] = r;
    }

    bool same(long u, long v) { return leader(u) == leader(v); }

    // 根を基準とした u のポテンシャルを返す
    T weight(long u) {
        leader(u);
        return diff_weight_[u];
    }

    // weight(v) - weight(u) を返す (事前条件: same(u, v) == true)
    T diff(long u, long v) {
        assert(same(u, v));
        return weight(v) - weight(u);
    }

    // weight(v) - weight(u) == w となるように併合する
    // すでに同じ連結成分にある場合、矛盾しなければ true、矛盾すれば false を返す
    bool merge(long u, long v, T w) {
        long ru = leader(u);
        long rv = leader(v);
        if (ru == rv) return diff_weight_[v] - diff_weight_[u] == w;
        w += diff_weight_[u];
        w -= diff_weight_[v];
        if (-data_[ru] < -data_[rv]) {
            std::swap(ru, rv);
            w = -w;
        }
        data_[ru] += data_[rv];
        data_[rv] = ru;
        diff_weight_[rv] = w;
        --comp_;
        return true;
    }

    long size(long u) { return -data_[leader(u)]; }

    long count() const { return comp_; }

   private:
    std::vector<long> data_;
    std::vector<T> diff_weight_;
    long comp_ = 0;
};

// Undo / Rollback 可能な Union-Find (経路圧縮なし + Union by Size)
// 空間計算量: O(N + Q) (Q は merge 回数)
// 時間計算量: merge, leader, same, size は最悪 O(log N)、undo, count, time は O(1)
struct UndoUnionFind {
    UndoUnionFind() = default;
    explicit UndoUnionFind(long n) : data_(n, -1), comp_(n) { assert(n >= 0); }

    long leader(long u) const {
        assert(0 <= u && u < std::ssize(data_));
        while (data_[u] >= 0)
            u = data_[u];
        return u;
    }

    bool same(long u, long v) const { return leader(u) == leader(v); }

    // merge 1回につき履歴を1つ積む (失敗時もダミーを積むため、1 merge に 1 undo が対応)
    bool merge(long u, long v) {
        u = leader(u);
        v = leader(v);
        if (u == v) {
            history_.push_back({-1, 0});
            return false;
        }
        if (-data_[u] < -data_[v]) std::swap(u, v);
        history_.push_back({v, data_[v]});
        data_[u] += data_[v];
        data_[v] = u;
        --comp_;
        return true;
    }

    void undo() {
        assert(!history_.empty());
        auto [v, dv] = history_.back();
        history_.pop_back();
        if (v != -1) {
            long u = data_[v];
            data_[u] -= dv;
            data_[v] = dv;
            ++comp_;
        }
    }

    // 現在の履歴位置 (スナップショット) を返す
    long time() const { return std::ssize(history_); }

    // 履歴位置 t の状態まで巻き戻す (事前条件: 0 <= t <= time())
    void rollback(long t) {
        assert(0 <= t && t <= std::ssize(history_));
        while (std::ssize(history_) > t)
            undo();
    }

    long size(long u) const { return -data_[leader(u)]; }

    long count() const { return comp_; }

   private:
    std::vector<long> data_;
    std::vector<std::pair<long, long>> history_;
    long comp_ = 0;
};

// 部分永続 Union-Find (経路圧縮なし + Union by Size)
// 空間計算量: O(N + Q) (Q は merge 回数)
// 時間計算量: merge, leader(u, t), same(u, v, t), size(u, t) は最悪 O(log N)
struct PartiallyPersistentUnionFind {
    PartiallyPersistentUnionFind() = default;
    explicit PartiallyPersistentUnionFind(long n) : data_(n, -1), time_(n, inf_), size_hist_(n, {{0, 1}}) { assert(n >= 0); }

    // 時刻 t (0 <= t <= now()) における u の代表元を返す
    long leader(long u, long t) const {
        assert(0 <= u && u < std::ssize(data_));
        assert(0 <= t && t <= now_);
        while (time_[u] <= t)
            u = data_[u];
        return u;
    }

    bool same(long u, long v, long t) const { return leader(u, t) == leader(v, t); }

    // 時刻 t における u の連結成分のサイズを返す
    long size(long u, long t) const {
        u = leader(u, t);
        auto it = std::upper_bound(size_hist_[u].begin(), size_hist_[u].end(), std::pair<long, long>{t, std::numeric_limits<long>::max()});
        return std::prev(it)->second;
    }

    // 時刻を 1 進めて u と v を併合し、新しい時刻を返す
    long merge(long u, long v) {
        ++now_;
        u = leader(u, now_);
        v = leader(v, now_);
        if (u == v) return now_;
        if (-data_[u] < -data_[v]) std::swap(u, v);
        data_[u] += data_[v];
        data_[v] = u;
        time_[v] = now_;
        size_hist_[u].push_back({now_, -data_[u]});
        return now_;
    }

    long now() const { return now_; }

   private:
    static constexpr long inf_ = std::numeric_limits<long>::max() / 4;
    std::vector<long> data_;
    std::vector<long> time_;
    std::vector<std::vector<std::pair<long, long>>> size_hist_;
    long now_ = 0;
};
