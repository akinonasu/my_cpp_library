#include <algorithm>
#include <cassert>
#include <iostream>
#include <print>
#include <string>
#include <utility>
#include <vector>

struct GridRect {
    int r1 = -1, c1 = -1, r2 = -1, c2 = -1;
    bool is_empty() const { return r1 == -1; }
    int height() const { return is_empty() ? 0 : r2 - r1; }
    int width() const { return is_empty() ? 0 : c2 - c1; }
};

class StringGrid {
   public:
    int H = 0;
    int W = 0;
    std::vector<std::string> g;

    static constexpr int DR4[4] = {-1, 0, 1, 0};
    static constexpr int DC4[4] = {0, 1, 0, -1};
    static constexpr int DR8[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    static constexpr int DC8[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

    StringGrid() = default;
    StringGrid(int h, int w, char fill = '.') : H(h), W(w), g(h, std::string(w, fill)) { assert(h >= 0 && w >= 0); }
    StringGrid(const std::vector<std::string>& grid) : g(grid) {
        H = static_cast<int>(grid.size());
        W = H > 0 ? static_cast<int>(grid[0].size()) : 0;
        for (const auto& row : g) {
            assert(static_cast<int>(row.size()) == W && "All rows must have the same length");
        }
    }

    static StringGrid read(int h, int w) {
        assert(h >= 0 && w >= 0);
        std::vector<std::string> s(h);
        for (int i = 0; i < h; ++i) {
            std::cin >> s[i];
        }
        return StringGrid(s);
    }

    static StringGrid from_coordinates(const std::vector<std::pair<int, int>>& coords) {
        if (coords.empty()) return StringGrid(0, 0);
        int min_r = coords[0].first, max_r = coords[0].first;
        int min_c = coords[0].second, max_c = coords[0].second;
        for (const auto& [r, c] : coords) {
            min_r = std::min(min_r, r);
            max_r = std::max(max_r, r);
            min_c = std::min(min_c, c);
            max_c = std::max(max_c, c);
        }
        int h = max_r - min_r + 1;
        int w = max_c - min_c + 1;
        StringGrid res(h, w, '.');
        for (const auto& [r, c] : coords) {
            res(r - min_r, c - min_c) = '#';
        }
        return res;
    }

    bool in_bounds(int r, int c) const { return 0 <= r && r < H && 0 <= c && c < W; }
    std::string& operator[](int r) { return g[r]; }
    const std::string& operator[](int r) const { return g[r]; }
    char& operator()(int r, int c) { return g[r][c]; }
    char operator()(int r, int c) const { return g[r][c]; }
    auto operator<=>(const StringGrid&) const = default;

    // ==================== 幾何操作 (基本) ====================

    StringGrid rotate90() const {
        std::vector<std::string> res(W, std::string(H, ' '));
        for (int r = 0; r < H; ++r) {
            for (int c = 0; c < W; ++c)
                res[c][H - 1 - r] = g[r][c];
        }
        return StringGrid(res);
    }

    StringGrid rotate180() const {
        std::vector<std::string> res = g;
        std::ranges::reverse(res);
        for (auto& row : res)
            std::ranges::reverse(row);
        return StringGrid(res);
    }

    StringGrid rotate270() const {
        std::vector<std::string> res(W, std::string(H, ' '));
        for (int r = 0; r < H; ++r) {
            for (int c = 0; c < W; ++c)
                res[W - 1 - c][r] = g[r][c];
        }
        return StringGrid(res);
    }

    StringGrid flip_h() const {
        std::vector<std::string> res = g;
        for (auto& row : res)
            std::ranges::reverse(row);
        return StringGrid(res);
    }

    StringGrid flip_v() const {
        std::vector<std::string> res = g;
        std::ranges::reverse(res);
        return StringGrid(res);
    }

    std::vector<StringGrid> all_transformations(bool unique = false) const {
        std::vector<StringGrid> res;
        res.reserve(8);
        StringGrid cur = *this;
        for (int i = 0; i < 4; ++i) {
            res.push_back(cur);
            res.push_back(cur.flip_h());
            cur = cur.rotate90();
        }
        if (unique) {
            std::ranges::sort(res);
            auto [first, last] = std::ranges::unique(res);
            res.erase(first, last);
        }
        return res;
    }

    StringGrid subgrid(int r1, int c1, int r2, int c2) const {
        assert(0 <= r1 && r1 <= r2 && r2 <= H);
        assert(0 <= c1 && c1 <= c2 && c2 <= W);
        std::vector<std::string> sub(r2 - r1);
        for (int r = r1; r < r2; ++r)
            sub[r - r1] = g[r].substr(c1, c2 - c1);
        return StringGrid(sub);
    }

    GridRect bounding_box(char ignore_char = '.') const {
        int r_min = H, r_max = -1, c_min = W, c_max = -1;
        for (int r = 0; r < H; ++r) {
            for (int c = 0; c < W; ++c) {
                if (g[r][c] != ignore_char) {
                    r_min = std::min(r_min, r);
                    r_max = std::max(r_max, r);
                    c_min = std::min(c_min, c);
                    c_max = std::max(c_max, c);
                }
            }
        }
        if (r_max == -1) return {};
        return {r_min, c_min, r_max + 1, c_max + 1};
    }

    StringGrid trim(char ignore_char = '.') const {
        auto rect = bounding_box(ignore_char);
        if (rect.is_empty()) return StringGrid(0, 0);
        return subgrid(rect.r1, rect.c1, rect.r2, rect.c2);
    }

    void paste(int r, int c, const StringGrid& other, char transparent_char = '\0') {
        for (int dr = 0; dr < other.H; ++dr) {
            for (int dc = 0; dc < other.W; ++dc) {
                int nr = r + dr, nc = c + dc;
                if (in_bounds(nr, nc)) {
                    char ch = other(dr, dc);
                    if (ch != transparent_char) g[nr][nc] = ch;
                }
            }
        }
    }

    // ==================== 拡大・周期的充填 (Scale & Tiling) ====================

    // 各マスを縦 kr 倍、横 kc 倍に拡大
    StringGrid scale(int kr, int kc) const {
        assert(kr >= 1 && kc >= 1);
        StringGrid res(H * kr, W * kc);
        for (int r = 0; r < H; ++r) {
            for (int c = 0; c < W; ++c) {
                char ch = g[r][c];
                for (int dr = 0; dr < kr; ++dr) {
                    for (int dc = 0; dc < kc; ++dc) {
                        res[r * kr + dr][c * kc + dc] = ch;
                    }
                }
            }
        }
        return res;
    }
    StringGrid scale(int k) const { return scale(k, k); }

    // グリッド全体を縦 rep_r 回、横 rep_c 回繰り返す
    StringGrid repeat(int rep_r, int rep_c) const {
        assert(rep_r >= 1 && rep_c >= 1);
        StringGrid res(H * rep_r, W * rep_c);
        for (int rr = 0; rr < rep_r; ++rr) {
            for (int rc = 0; rc < rep_c; ++rc) {
                for (int r = 0; r < H; ++r) {
                    for (int c = 0; c < W; ++c) {
                        res[rr * H + r][rc * W + c] = g[r][c];
                    }
                }
            }
        }
        return res;
    }

    // target_h x target_w の盤面を、本パターンの周期充填で埋める
    StringGrid tile_to(int target_h, int target_w, int offset_r = 0, int offset_c = 0) const {
        assert(target_h >= 0 && target_w >= 0 && H > 0 && W > 0);
        StringGrid res(target_h, target_w);
        for (int r = 0; r < target_h; ++r) {
            int orig_r = ((r + offset_r) % H + H) % H;
            for (int c = 0; c < target_w; ++c) {
                int orig_c = ((c + offset_c) % W + W) % W;
                res[r][c] = g[orig_r][orig_c];
            }
        }
        return res;
    }

    // ==================== 小盤面での平面充填ソルバー (Exact Tiling) ====================

   private:
    // ピースの相対座標（辞書順最小マスを (0, 0) としたもの）
    using PieceShape = std::vector<std::pair<int, int>>;

    static PieceShape extract_normalized_shape(const StringGrid& piece, char ignore_char = '.') {
        std::vector<std::pair<int, int>> pts;
        for (int r = 0; r < piece.H; ++r) {
            for (int c = 0; c < piece.W; ++c) {
                if (piece(r, c) != ignore_char) pts.emplace_back(r, c);
            }
        }
        if (pts.empty()) return {};
        std::ranges::sort(pts);
        auto [base_r, base_c] = pts[0];
        for (auto& [r, c] : pts) {
            r -= base_r;
            c -= base_c;
        }
        return pts;
    }

   public:
    // 盤面上の empty_char マスを、pieces の全パーツで過不足なく敷き詰められるか判定
    // 返り値: {成否, 敷き詰め結果のグリッド}
    std::pair<bool, StringGrid> solve_tiling(const std::vector<StringGrid>& pieces, char empty_char = '.', bool allow_rotate = true, bool allow_flip = true) const {
        int empty_cells = 0;
        for (int r = 0; r < H; ++r) {
            for (int c = 0; c < W; ++c) {
                if (g[r][c] == empty_char) empty_cells++;
            }
        }

        int n_pieces = static_cast<int>(pieces.size());
        std::vector<std::vector<PieceShape>> piece_variants(n_pieces);
        int total_piece_cells = 0;

        for (int i = 0; i < n_pieces; ++i) {
            std::vector<StringGrid> forms;
            if (allow_rotate && allow_flip)
                forms = pieces[i].all_transformations(true);
            else if (allow_rotate) {
                StringGrid cur = pieces[i];
                for (int rot = 0; rot < 4; ++rot) {
                    forms.push_back(cur);
                    cur = cur.rotate90();
                }
            } else if (allow_flip) {
                forms = {pieces[i], pieces[i].flip_h(), pieces[i].flip_v(), pieces[i].flip_h().flip_v()};
            } else {
                forms = {pieces[i]};
            }

            for (const auto& form : forms) {
                PieceShape shape = extract_normalized_shape(form, '.');
                if (!shape.empty()) piece_variants[i].push_back(shape);
            }
            std::ranges::sort(piece_variants[i]);
            auto [first, last] = std::ranges::unique(piece_variants[i]);
            piece_variants[i].erase(first, last);

            if (!piece_variants[i].empty()) {
                total_piece_cells += static_cast<int>(piece_variants[i][0].size());
            }
        }

        // マス数の合計が合わなければ即座に不可
        if (empty_cells != total_piece_cells) {
            return {false, StringGrid(0, 0)};
        }

        StringGrid cur_grid = *this;
        std::vector<bool> used(n_pieces, false);

        // ラベル文字（A, B, ..., Z, a, b, ..., 0, 1, ...）
        auto get_label = [](int idx) -> char {
            if (idx < 26) return static_cast<char>('A' + idx);
            if (idx < 52) return static_cast<char>('a' + (idx - 26));
            return static_cast<char>('0' + (idx % 10));
        };

        // バックトラッキング探索（左上優先の枝刈り）
        auto dfs = [&](auto self, int start_r, int start_c, int placed_cnt) -> bool {
            if (placed_cnt == n_pieces) return true;

            // 辞書順で最初の空きマスを探す
            int target_r = -1, target_c = -1;
            for (int r = start_r; r < H; ++r) {
                int c_begin = (r == start_r) ? start_c : 0;
                for (int c = c_begin; c < W; ++c) {
                    if (cur_grid(r, c) == empty_char) {
                        target_r = r;
                        target_c = c;
                        break;
                    }
                }
                if (target_r != -1) break;
            }

            if (target_r == -1) return (placed_cnt == n_pieces);

            for (int i = 0; i < n_pieces; ++i) {
                if (used[i]) continue;
                // 同一形状ピースの順序固定（対称性の枝刈り）
                if (i > 0 && !used[i - 1] && piece_variants[i] == piece_variants[i - 1]) {
                    continue;
                }

                used[i] = true;
                char label = get_label(i);

                for (const auto& shape : piece_variants[i]) {
                    bool can_place = true;
                    for (const auto& [dr, dc] : shape) {
                        int nr = target_r + dr;
                        int nc = target_c + dc;
                        if (!cur_grid.in_bounds(nr, nc) || cur_grid(nr, nc) != empty_char) {
                            can_place = false;
                            break;
                        }
                    }

                    if (can_place) {
                        for (const auto& [dr, dc] : shape)
                            cur_grid(target_r + dr, target_c + dc) = label;
                        if (self(self, target_r, target_c, placed_cnt + 1)) return true;
                        for (const auto& [dr, dc] : shape)
                            cur_grid(target_r + dr, target_c + dc) = empty_char;
                    }
                }
                used[i] = false;
            }
            return false;
        };

        if (dfs(dfs, 0, 0, 0)) {
            return {true, cur_grid};
        }
        return {false, StringGrid(0, 0)};
    }

    // 指定の単一ピース複数個で盤面の空きマス (empty_char) を過不足なく埋められるか
    std::pair<bool, StringGrid> tile_with_single_piece(const StringGrid& piece, char empty_char = '.', bool allow_rotate = true, bool allow_flip = true) const {
        int empty_cnt = 0;
        for (int r = 0; r < H; ++r) {
            for (int c = 0; c < W; ++c) {
                if (g[r][c] == empty_char) empty_cnt++;
            }
        }
        int piece_cnt = 0;
        for (int r = 0; r < piece.H; ++r) {
            for (int c = 0; c < piece.W; ++c) {
                if (piece(r, c) != '.') piece_cnt++;
            }
        }
        if (piece_cnt == 0 || empty_cnt % piece_cnt != 0) {
            return {false, StringGrid(0, 0)};
        }
        int need = empty_cnt / piece_cnt;
        std::vector<StringGrid> pieces(need, piece);
        return solve_tiling(pieces, empty_char, allow_rotate, allow_flip);
    }

    // ==================== 定番ピースのプリセット生成 ====================
    static StringGrid domino() { return StringGrid({"##"}); }
    static StringGrid tromino_L() { return StringGrid({"#.", "##"}); }
    static StringGrid tromino_I() { return StringGrid({"###"}); }
    static StringGrid tetromino(char type) {
        switch (type) {
            case 'I':
                return StringGrid({"####"});
            case 'O':
                return StringGrid({"##", "##"});
            case 'T':
                return StringGrid({".#.", "###"});
            case 'L':
                return StringGrid({"#..", "###"});
            case 'J':
                return StringGrid({"..#", "###"});
            case 'S':
                return StringGrid({".##", "##."});
            case 'Z':
                return StringGrid({"##.", ".##"});
            default:
                assert(false && "Invalid tetromino type");
        }
        return {};
    }

    void print() const {
        for (const auto& row : g) {
            std::println("{}", row);
        }
    }
};
