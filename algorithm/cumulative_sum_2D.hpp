#include <cassert>
#include <print>
#include <vector>

// ==========================================
// 2次元累積和 (2D Cumulative Sum)
// ==========================================
template <typename T = long>
class CumulativeSum2D {
   private:
    int H, W;
    std::vector<std::vector<T>> data;
    bool built = false;

   public:
    CumulativeSum2D() : H(0), W(0), built(false) {}
    CumulativeSum2D(int H, int W) : H(H), W(W), data(H + 1, std::vector<T>(W + 1, 0)), built(false) {}

    // 初期グリッドから一括構築
    CumulativeSum2D(const std::vector<std::vector<T>>& grid) : H(grid.size()), W(grid.empty() ? 0 : grid[0].size()), data(H + 1, std::vector<T>(W + 1, 0)), built(false) {
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                data[i + 1][j + 1] = grid[i][j];
            }
        }
        build();
    }

    // セル (r, c) に値を設定 (build前のみ呼出可能)
    void set(int r, int c, T val) {
        assert(!built && "Cannot call set() after build()");
        assert(0 <= r && r < H && 0 <= c && c < W);
        data[r + 1][c + 1] = val;
    }

    // セル (r, c) に値を加算 (build前のみ呼出可能)
    void add(int r, int c, T val) {
        assert(!built && "Cannot call add() after build()");
        assert(0 <= r && r < H && 0 <= c && c < W);
        data[r + 1][c + 1] += val;
    }

    // 累積和テーブルを構築
    void build() {
        assert(!built && "build() has already been called");
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                data[i + 1][j + 1] += data[i + 1][j] + data[i][j + 1] - data[i][j];
            }
        }
        built = true;
    }

    // 半開区間 [r1, r2) x [c1, c2) の総和を取得
    T sum(int r1, int c1, int r2, int c2) const {
        assert(built && "Must call build() before sum()");
        assert(0 <= r1 && r1 <= r2 && r2 <= H);
        assert(0 <= c1 && c1 <= c2 && c2 <= W);
        return data[r2][c2] - data[r1][c2] - data[r2][c1] + data[r1][c1];
    }

    // 特定セルの値を取得
    T get(int r, int c) const { return sum(r, c, r + 1, c + 1); }
};

// ==========================================
// 2次元いもす法 (2D Imos)
// ==========================================
template <typename T = long>
class Imos2D {
   private:
    int H, W;
    std::vector<std::vector<T>> diff;
    bool built = false;

   public:
    Imos2D(int H, int W) : H(H), W(W), diff(H + 1, std::vector<T>(W + 1, 0)), built(false) {}

    // 半開区間 [r1, r2) x [c1, c2) に val を加算 (build前のみ呼出可能)
    void add(int r1, int c1, int r2, int c2, T val) {
        assert(!built && "Cannot call add() after build()");
        assert(0 <= r1 && r1 <= r2 && r2 <= H);
        assert(0 <= c1 && c1 <= c2 && c2 <= W);
        if (r1 == r2 || c1 == c2) return;

        diff[r1][c1] += val;
        diff[r1][c2] -= val;
        diff[r2][c1] -= val;
        diff[r2][c2] += val;
    }

    // 二次元累積和を行って各セルの最終値を確定
    void build() {
        assert(!built && "build() has already been called");

        // 横方向の累積和
        for (int i = 0; i <= H; ++i) {
            for (int j = 0; j < W; ++j) {
                diff[i][j + 1] += diff[i][j];
            }
        }
        // 縦方向の累積和
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j <= W; ++j) {
                diff[i + 1][j] += diff[i][j];
            }
        }

        // 不要な境界部分を削除して H x W サイズにする
        diff.pop_back();
        for (auto& row : diff) {
            row.pop_back();
        }
        built = true;
    }

    // 確定後のセル (r, c) の値を取得
    T get(int r, int c) const {
        assert(built && "Must call build() before get()");
        assert(0 <= r && r < H && 0 <= c && c < W);
        return diff[r][c];
    }

    // [] 演算子によるアクセス (imos[r][c])
    const std::vector<T>& operator[](int r) const {
        assert(built && "Must call build() before access");
        assert(0 <= r && r < H);
        return diff[r];
    }

    // 全結果の2次元配列を取得
    const std::vector<std::vector<T>>& result() const {
        assert(built && "Must call build() before result()");
        return diff;
    }

    // 計算結果から CumulativeSum2D を構築して返す
    CumulativeSum2D<T> to_cumulative_sum() const {
        assert(built && "Must call build() before converting");
        return CumulativeSum2D<T>(diff);
    }
};
