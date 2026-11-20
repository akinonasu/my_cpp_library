#pragma once

#include <cassert>
#include <vector>

// ==========================================
// 2次元累積和 (2D Cumulative Sum)
// 時間計算量: 構築 O(HW), クエリ O(1)
// 空間計算量: O(HW)
// 事前条件:
//   - 0 <= h, 0 <= w
//   - grid から構築する場合、各行の長さはすべて等しいこと
//   - set / add は build() 前、sum / get は build() 後に呼び出すこと
// ==========================================
template <class T = long>
struct CumulativeSum2D {
    CumulativeSum2D() : CumulativeSum2D(0, 0) {}

    CumulativeSum2D(long h, long w) : h_(h), w_(w), data_(h + 1, std::vector<T>(w + 1, T{})), built_(false) { assert(0 <= h && 0 <= w); }

    CumulativeSum2D(const std::vector<std::vector<T>>& grid) : h_(std::ssize(grid)), w_(grid.empty() ? 0 : std::ssize(grid[0])), data_(h_ + 1, std::vector<T>(w_ + 1, T{})), built_(false) {
        for (long i = 0; i < h_; ++i) {
            assert(std::ssize(grid[i]) == w_);
            for (long j = 0; j < w_; ++j) {
                data_[i + 1][j + 1] = grid[i][j];
            }
        }
        build();
    }

    void set(long r, long c, T val) {
        assert(!built_);
        assert(0 <= r && r < h_ && 0 <= c && c < w_);
        data_[r + 1][c + 1] = val;
    }

    void add(long r, long c, T val) {
        assert(!built_);
        assert(0 <= r && r < h_ && 0 <= c && c < w_);
        data_[r + 1][c + 1] += val;
    }

    void build() {
        assert(!built_);
        for (long i = 0; i < h_; ++i) {
            for (long j = 0; j < w_; ++j) {
                data_[i + 1][j + 1] += data_[i + 1][j] + data_[i][j + 1] - data_[i][j];
            }
        }
        built_ = true;
    }

    // 半開区間 [r1, r2) x [c1, c2) の総和を取得
    T sum(long r1, long c1, long r2, long c2) const {
        assert(built_);
        assert(0 <= r1 && r1 <= r2 && r2 <= h_);
        assert(0 <= c1 && c1 <= c2 && c2 <= w_);
        return data_[r2][c2] - data_[r1][c2] - data_[r2][c1] + data_[r1][c1];
    }

    T get(long r, long c) const { return sum(r, c, r + 1, c + 1); }

   private:
    long h_;
    long w_;
    std::vector<std::vector<T>> data_;
    bool built_;
};

// ==========================================
// 2次元いもす法 (2D Imos)
// 時間計算量: add O(1), build O(HW), get / operator[] O(1), to_cumulative_sum O(HW)
// 空間計算量: O(HW)
// 事前条件:
//   - 0 <= h, 0 <= w
//   - add は build() 前、get / operator[] / result / to_cumulative_sum は build() 後に呼び出すこと
// ==========================================
template <class T = long>
struct Imos2D {
    Imos2D() : Imos2D(0, 0) {}

    Imos2D(long h, long w) : h_(h), w_(w), diff_(h + 1, std::vector<T>(w + 1, T{})), built_(false) { assert(0 <= h && 0 <= w); }

    // 半開区間 [r1, r2) x [c1, c2) に val を加算
    void add(long r1, long c1, long r2, long c2, T val) {
        assert(!built_);
        assert(0 <= r1 && r1 <= r2 && r2 <= h_);
        assert(0 <= c1 && c1 <= c2 && c2 <= w_);
        if (r1 == r2 || c1 == c2) return;

        diff_[r1][c1] += val;
        diff_[r1][c2] -= val;
        diff_[r2][c1] -= val;
        diff_[r2][c2] += val;
    }

    void build() {
        assert(!built_);
        for (long i = 0; i < h_; ++i) {
            for (long j = 0; j + 1 < w_; ++j) {
                diff_[i][j + 1] += diff_[i][j];
            }
        }
        for (long i = 0; i + 1 < h_; ++i) {
            for (long j = 0; j < w_; ++j) {
                diff_[i + 1][j] += diff_[i][j];
            }
        }
        diff_.pop_back();
        for (auto& row : diff_) {
            row.pop_back();
        }
        built_ = true;
    }

    T get(long r, long c) const {
        assert(built_);
        assert(0 <= r && r < h_ && 0 <= c && c < w_);
        return diff_[r][c];
    }

    const std::vector<T>& operator[](long r) const {
        assert(built_);
        assert(0 <= r && r < h_);
        return diff_[r];
    }

    const std::vector<std::vector<T>>& result() const {
        assert(built_);
        return diff_;
    }

    CumulativeSum2D<T> to_cumulative_sum() const {
        assert(built_);
        return CumulativeSum2D<T>(diff_);
    }

   private:
    long h_;
    long w_;
    std::vector<std::vector<T>> diff_;
    bool built_;
};