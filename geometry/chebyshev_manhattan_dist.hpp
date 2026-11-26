#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

// 2点間のマンハッタン距離 |x1 - x2| + |y1 - y2|
// 時間計算量: O(1), 空間計算量: O(1)
inline long manhattan_dist(long x1, long y1, long x2, long y2) {
    return std::abs(x1 - x2) + std::abs(y1 - y2);
}

inline long manhattan_dist(std::pair<long, long> p, std::pair<long, long> q) {
    return manhattan_dist(p.first, p.second, q.first, q.second);
}

inline long double manhattan_dist_real(std::pair<long double, long double> p, std::pair<long double, long double> q) {
    return std::abs(p.first - q.first) + std::abs(p.second - q.second);
}

// 2点間のチェビシェフ距離 max(|x1 - x2|, |y1 - y2|)
// 時間計算量: O(1), 空間計算量: O(1)
inline long chebyshev_dist(long x1, long y1, long x2, long y2) {
    return std::max(std::abs(x1 - x2), std::abs(y1 - y2));
}

inline long chebyshev_dist(std::pair<long, long> p, std::pair<long, long> q) {
    return chebyshev_dist(p.first, p.second, q.first, q.second);
}

inline long double chebyshev_dist_real(std::pair<long double, long double> p, std::pair<long double, long double> q) {
    return std::max(std::abs(p.first - q.first), std::abs(p.second - q.second));
}

// マンハッタン距離 -> チェビシェフ距離（等倍変換）
// (x, y) -> (x + y, x - y)
// d_manhattan(p, q) == d_chebyshev(f(p), f(q))
// 時間計算量: 点 O(1) / 配列 O(N), 空間計算量: 点 O(1) / 配列 O(N)
// 事前条件: x + y, x - y が long に収まること
inline std::pair<long, long> manhattan_to_chebyshev(long x, long y) {
    return {x + y, x - y};
}

inline std::pair<long, long> manhattan_to_chebyshev(std::pair<long, long> p) {
    return manhattan_to_chebyshev(p.first, p.second);
}

template <class T>
inline std::vector<std::pair<long, long>> manhattan_to_chebyshev(const std::vector<std::pair<T, T>>& pts) {
    std::vector<std::pair<long, long>> res;
    res.reserve(pts.size());
    for (const auto& p : pts) {
        res.push_back(manhattan_to_chebyshev(long(p.first), long(p.second)));
    }
    return res;
}

// チェビシェフ空間の座標 (u, v) が元の整数格子点へ復元可能か判定 (u ≡ v mod 2)
// 時間計算量: O(1), 空間計算量: O(1)
inline bool can_restore_from_chebyshev(long u, long v) {
    return ((u ^ v) & 1) == 0;
}

inline bool can_restore_from_chebyshev(std::pair<long, long> p) {
    return can_restore_from_chebyshev(p.first, p.second);
}

// チェビシェフ空間の座標 (u, v) から元の整数座標 ((u + v) / 2, (u - v) / 2) へ復元
// 復元できない場合は {false, {0, 0}} を返す
// 時間計算量: O(1), 空間計算量: O(1)
inline std::pair<bool, std::pair<long, long>> restore_from_chebyshev(long u, long v) {
    if (!can_restore_from_chebyshev(u, v)) return {false, {0, 0}};
    long x = long((__int128_t(u) + v) / 2);
    long y = long((__int128_t(u) - v) / 2);
    return {true, {x, y}};
}

inline std::pair<bool, std::pair<long, long>> restore_from_chebyshev(std::pair<long, long> p) {
    return restore_from_chebyshev(p.first, p.second);
}

// チェビシェフ距離 -> マンハッタン距離（整数 2倍スケール変換）
// (x, y) -> (x + y, x - y)
// 2 * d_chebyshev(p, q) == d_manhattan(f(p), f(q))
// 時間計算量: 点 O(1) / 配列 O(N), 空間計算量: 点 O(1) / 配列 O(N)
// 事前条件: x + y, x - y が long に収まること
inline std::pair<long, long> chebyshev_to_manhattan_2x(long x, long y) {
    return {x + y, x - y};
}

inline std::pair<long, long> chebyshev_to_manhattan_2x(std::pair<long, long> p) {
    return chebyshev_to_manhattan_2x(p.first, p.second);
}

template <class T>
inline std::vector<std::pair<long, long>> chebyshev_to_manhattan_2x(const std::vector<std::pair<T, T>>& pts) {
    std::vector<std::pair<long, long>> res;
    res.reserve(pts.size());
    for (const auto& p : pts) {
        res.push_back(chebyshev_to_manhattan_2x(long(p.first), long(p.second)));
    }
    return res;
}

// 2倍スケール変換された座標 (u, v) から元の整数座標へ復元
// 時間計算量: O(1), 空間計算量: O(1)
inline std::pair<bool, std::pair<long, long>> restore_from_manhattan_2x(long u, long v) {
    return restore_from_chebyshev(u, v);
}

inline std::pair<bool, std::pair<long, long>> restore_from_manhattan_2x(std::pair<long, long> p) {
    return restore_from_chebyshev(p);
}

// チェビシェフ距離 -> マンハッタン距離（実数等倍変換）
// (x, y) -> ((x + y) / 2, (x - y) / 2)
// d_chebyshev(p, q) == d_manhattan(f(p), f(q))
// 時間計算量: 点 O(1) / 配列 O(N), 空間計算量: 点 O(1) / 配列 O(N)
inline std::pair<long double, long double> chebyshev_to_manhattan(long double x, long double y) {
    return {(x + y) * 0.5L, (x - y) * 0.5L};
}

inline std::pair<long double, long double> chebyshev_to_manhattan(std::pair<long double, long double> p) {
    return chebyshev_to_manhattan(p.first, p.second);
}

template <class T>
inline std::vector<std::pair<long double, long double>> chebyshev_to_manhattan(const std::vector<std::pair<T, T>>& pts) {
    std::vector<std::pair<long double, long double>> res;
    res.reserve(pts.size());
    for (const auto& p : pts) {
        res.push_back(chebyshev_to_manhattan(p.first, p.second));
    }
    return res;
}

// 実数等倍変換された座標 (u, v) から元の座標 (u + v, u - v) へ復元
// 時間計算量: O(1), 空間計算量: O(1)
inline std::pair<long double, long double> restore_from_manhattan_real(long double u, long double v) {
    return {u + v, u - v};
}

inline std::pair<long double, long double> restore_from_manhattan_real(std::pair<long double, long double> p) {
    return restore_from_manhattan_real(p.first, p.second);
}

// 半開矩形領域 [xl, xr) x [yl, yr)
struct Rect {
    long xl, xr, yl, yr;
};

// 中心 (x, y)、マンハッタン距離 r 以下の領域をチェビシェフ空間の半開矩形へ変換
// 注意: 矩形内の整数点 (u, v) のうち、元の格子点に対応するのは u ≡ v (mod 2) の点のみ
// 時間計算量: O(1), 空間計算量: O(1)
// 事前条件: r >= 0
inline Rect manhattan_disk_to_chebyshev_rect(long x, long y, long r) {
    assert(r >= 0);
    auto [u, v] = manhattan_to_chebyshev(x, y);
    return {u - r, u + r + 1, v - r, v + r + 1};
}

// 点集合に対する最大マンハッタン距離および直径を O(1) で求める構造体
// 空間計算量: O(1)
struct ManhattanFarthest {
    void add(long x, long y) {
        auto [u, v] = manhattan_to_chebyshev(x, y);
        if (cnt_ == 0) {
            min_u_ = max_u_ = u;
            min_v_ = max_v_ = v;
        } else {
            min_u_ = std::min(min_u_, u);
            max_u_ = std::max(max_u_, u);
            min_v_ = std::min(min_v_, v);
            max_v_ = std::max(max_v_, v);
        }
        ++cnt_;
    }

    void add(std::pair<long, long> p) { add(p.first, p.second); }

    // 登録済み点群から (x, y) への最大マンハッタン距離
    // 時間計算量: O(1), 事前条件: !empty()
    long max_dist(long x, long y) const {
        assert(cnt_ > 0);
        auto [u, v] = manhattan_to_chebyshev(x, y);
        return std::max({u - min_u_, max_u_ - u, v - min_v_, max_v_ - v});
    }

    long max_dist(std::pair<long, long> p) const { return max_dist(p.first, p.second); }

    // 登録済み点群内の最大マンハッタン距離（直径）
    // 時間計算量: O(1), 事前条件: !empty()
    long diameter() const {
        assert(cnt_ > 0);
        return std::max(max_u_ - min_u_, max_v_ - min_v_);
    }

    long size() const { return cnt_; }
    bool empty() const { return cnt_ == 0; }
    void clear() { cnt_ = 0; }

    void dump() const {
#ifdef LOCAL
        if (cnt_ == 0) {
            std::cerr << "ManhattanFarthest(empty)\n";
            return;
        }
        std::cerr << "ManhattanFarthest(n=" << cnt_ << ", u=[" << min_u_ << ", " << max_u_ << "]"
                  << ", v=[" << min_v_ << ", " << max_v_ << "]"
                  << ", diameter=" << diameter() << ")\n";
#endif
    }

   private:
    long cnt_ = 0;
    long min_u_ = 0, max_u_ = 0;
    long min_v_ = 0, max_v_ = 0;
};
