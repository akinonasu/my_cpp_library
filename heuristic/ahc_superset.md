# `ahc_superset.hpp` 完全活用ガイド

AtCoder Heuristic Contest (AHC) 向けの 1 ファイル完結 C++23 スーパーセットライブラリ `ahc_superset.hpp` の使い方・設計意図・実戦テンプレート集です。

- **対応環境**: C++23 (`-std=c++23 -O2 -Wall -Wextra`)
- **ローカルデバッグ**: `-DLOCAL` を付けてコンパイルすると、各構造体の `dump()` が `std::cerr` に統計・内部状態を出力します。本番提出時（`-DLOCAL` なし）は空関数になり、`-O2` で完全に最適化消滅（ゼロオーバーヘッド）します。

---

## 目次

1. [収録コンポーネント一覧](#収録コンポーネント一覧)
2. [パターン別：実戦テンプレート](#パターン別実戦テンプレート)
   - [A. 焼きなまし法（SA）＋近傍選択の自動適応](#a-焼きなまし法sa近傍選択の自動適応)
   - [B. オイラーツアー型・差分更新ビームサーチ](#b-オイラーツアー型差分更新ビームサーチ)
   - [C. Chokudai サーチ（時間いっぱい自動拡張）](#c-chokudai-サーチ時間いっぱい自動拡張)
   - [D. 巡回経路最適化（2-opt / Or-opt）](#d-巡回経路最適化2-opt--or-opt)
   - [E. ゼロアロケーション反復 BFS（SA 内部用）](#e-ゼロアロケーション反復-bfssa-内部用)
3. [各コンポーネント詳細リファレンス](#各コンポーネント詳細リファレンス)
4. [AHC 本番で勝つためのチューニング Tips](#ahc-本番で勝つためのチューニング-tips)

---

## 収録コンポーネント一覧

| # | 構造体 | 役割 | 主な計算量 |
|---|---|---|---|
| 1 | `Timer` | 高精度経過時間・進捗率・時間切れ判定 | $O(1)$ |
| 2 | `XorShift` | `xoshiro256**` 超高速乱数・シャッフル | $O(1)$ / $O(N)$ |
| 3 | `ZobristTable` | 盤面・状態用の $O(1)$ 差分ハッシュテーブル | 参照 $O(1)$ |
| 4 | `StampArray` | $O(1)$ クリア可能な訪問済み配列 | 全操作 $O(1)$ |
| 5 | `FastQueue<T>` | ヒープ確保なしの固定長高速キュー | 全操作 $O(1)$ |
| 6 | `FastHashSet` | 世代管理付きオープンアドレス法ハッシュセット | 全操作期待 $O(1)$、`clear` $O(1)$ |
| 7 | `FastSet` | $[0, n)$ の集合（追加・削除・ランダム取得） | 全操作 $O(1)$ |
| 8 | `NeighborSelector` | 複数近傍の選択確率を改善率で自動調整 | `select` $O(K)$, `feedback` $O(1)$ |
| 9 | `Annealer` | 時間間引き・対数テーブル内蔵の焼きなまし | ならし $O(1)$ |
| 10 | `EulerTourBeam` | 状態コピーゼロの差分更新型ビームサーチ | 1ターン $O(W \cdot (\text{分岐} + \Delta))$ |
| 11 | `ChokudaiSearch` | 制限時間いっぱい自動で幅を広げる探索 | 時間制限まで反復 |
| 12 | `TourOptimizer` | 巡回経路向け 2-opt / Or-opt 局所最適化 | 1パス $O(N^2)$ |
| 13 | `Pos` / `Grid<T>` | 1次元フラット配列の高速 2D グリッド | アクセス $O(1)$ |

---

## パターン別：実戦テンプレート

### A. 焼きなまし法（SA）＋近傍選択の自動適応

`Annealer` と `NeighborSelector` を組み合わせることで、複数近傍のうち「今効いている近傍」を自動で優先しながら毎秒数千万イテレーション回せます。

```cpp
#include "ahc_superset.hpp"

void solve_sa() {
    XorShift rng(42);
    // 制限時間 1900ms, 開始温度 100.0, 終了温度 0.1, 最大化(true)
    Annealer sa(1900, 100.0L, 0.1L, true);
    NeighborSelector sel(3);  // 3種類の近傍を自動調整

    long double cur_score = 0.0L;
    long double best_score = cur_score;

    while (sa.step()) {
        int nb = sel.select(rng);
        long double delta = 0.0L;

        if (nb == 0) {
            // 近傍 0: 1点変更の差分スコア計算
            delta = /* ... */ 0.0L;
        } else if (nb == 1) {
            // 近傍 1: 2点スワップの差分スコア計算
            delta = /* ... */ 0.0L;
        } else {
            // 近傍 2: 区間反転などの差分スコア計算
            delta = /* ... */ 0.0L;
        }

        if (sa.accept(delta, rng)) {
            // 変更を確定
            cur_score += delta;
            if (cur_score > best_score) {
                best_score = cur_score;
                // ベスト解を保存
            }
            sel.feedback(nb, delta > 0.0L);
        } else {
            // 変更をロールバック（必要な場合）
            sel.feedback(nb, false);
        }
    }

    sa.dump();   // -DLOCAL 時のみイテレーション数・受理率等を stderr 出力
    sel.dump();  // -DLOCAL 時のみ各近傍の重み・改善回数を stderr 出力
}
```

> **ポイント**: 最小化問題の場合は `Annealer sa(1900, start_t, end_t, false);` と第 4 引数を `false` にします。`delta` は常に `new_score - cur_score` を渡してください。

---

### B. オイラーツアー型・差分更新ビームサーチ

盤面サイズが大きい問題で通常のビームサーチを行うと、状態コピーがボトルネックになります。  
`EulerTourBeam` は **単一の `State` に対して `apply` と `undo` だけを呼ぶ** ため、コピーコストゼロでビーム幅を大幅に広げられます。

```cpp
#include "ahc_superset.hpp"

struct Op {
    int pos = 0;
    int old_val = 0;
    int new_val = 0;
};

struct BoardState {
    Grid<int> grid;
    ZobristTable zt;
    long score = 0;
    std::uint64_t hash = 0;

    void apply(const Op& op) {
        hash ^= zt(op.pos, op.old_val) ^ zt(op.pos, op.new_val);
        Pos p = grid.pos(op.pos);
        grid[p] = op.new_val;
        score += (op.new_val - op.old_val);  // 差分スコア更新
    }

    void undo(const Op& op) {
        hash ^= zt(op.pos, op.new_val) ^ zt(op.pos, op.old_val);
        Pos p = grid.pos(op.pos);
        grid[p] = op.old_val;
        score -= (op.new_val - op.old_val);  // 差分スコア巻き戻し
    }

    // 現在状態からの遷移候補を push_cand(op, next_score, next_hash) で列挙
    template <class F>
    void expand(F push_cand) const {
        for (int idx = 0; idx < grid.size(); ++idx) {
            Pos p = grid.pos(idx);
            int cur_v = grid[p];
            int nxt_v = cur_v + 1;
            if (nxt_v >= 4) continue;
            Op op{idx, cur_v, nxt_v};
            long nxt_score = score + (nxt_v - cur_v);
            std::uint64_t nxt_hash = hash ^ zt(idx, cur_v) ^ zt(idx, nxt_v);
            push_cand(op, nxt_score, nxt_hash);
        }
    }
};

void solve_beam() {
    XorShift rng(123);
    BoardState state;
    state.grid = Grid<int>(5, 5, 0);
    state.zt = ZobristTable(25, 4, rng);
    for (int i = 0; i < 25; ++i) state.hash ^= state.zt(i, 0);

    // ビーム幅 200, ハッシュテーブル容量 2^18
    EulerTourBeam<BoardState, Op, long> beam(200, 1 << 18);
    // 最大 50 ターン, 制限時間 1800ms
    auto [best_score, ops] = beam.run(state, 50, 1800);
    beam.dump();

    // ops に最良解までの操作列が順番に入っている
}
```

---

### C. Chokudai サーチ（時間いっぱい自動拡張）

「ビーム幅をいくつにすれば制限時間ぴったりに終わるか分からない」というときに最強です。  
`State` に `eval()`, `hash()`, `expand(push_next)` の 3 つを実装するだけで使えます。

```cpp
#include "ahc_superset.hpp"

struct MyState {
    int pos = 0;
    long score = 0;
    std::uint64_t h = 0;

    long eval() const { return score; }
    std::uint64_t hash() const { return h; }

    template <class F>
    void expand(F push_next) const {
        for (int dir = 0; dir < 4; ++dir) {
            MyState nxt = *this;
            // nxt を更新 ...
            push_next(nxt, dir);
        }
    }
};

void solve_chokudai(const MyState& init) {
    // 1 反復あたり各深さから 2 ノード展開、ハッシュ容量 2^18
    ChokudaiSearch<MyState, int, long> cs(2, 1 << 18);
    auto [best_score, ops] = cs.run(init, /*max_turns=*/100, /*limit_ms=*/1800);
    cs.dump();
}
```

---

### D. 巡回経路最適化（2-opt / Or-opt）

配達・訪問順序の最適化（TSP 系）で初期解や最終仕上げに使えます。

```cpp
#include "ahc_superset.hpp"

void solve_tsp(const std::vector<Pos>& pts) {
    int n = static_cast<int>(std::ssize(pts));
    std::vector<int> tour(n);
    std::iota(tour.begin(), tour.end(), 0);

    auto dist = [&](int u, int v) -> long double {
        long double dr = pts[u].r - pts[v].r;
        long double dc = pts[u].c - pts[v].c;
        return std::sqrt(dr * dr + dc * dc);
    };

    // closed=true: 閉路（終点から始点に戻る）
    // closed=false: 開路（始点・終点固定）
    TourOptimizer::two_opt(tour, dist, true);
    TourOptimizer::or_opt(tour, dist, true, /*max_seg=*/3);

    long double cost = TourOptimizer::total_cost(tour, dist, true);
    (void)cost;
}
```

---

### E. ゼロアロケーション反復 BFS（SA 内部用）

焼きなましの中で毎回 `std::vector<bool> vis(N)` や `std::queue<int>` を作ると確保・初期化コストで激重になります。  
`StampArray` + `FastQueue` をループ外で 1 回だけ確保し、毎回 `clear()`（$O(1)$）するだけで爆速になります。

```cpp
#include "ahc_superset.hpp"

void fast_bfs_in_sa(const Grid<int>& wall) {
    int n = static_cast<int>(wall.size());
    StampArray vis(n);
    FastQueue<int> q(n);

    auto bfs = [&](Pos start) {
        vis.clear();  // O(1)
        q.clear();    // O(1)
        int sid = wall.id(start);
        vis.mark(sid);
        q.push(sid);

        while (!q.empty()) {
            int cur = q.pop();
            Pos p = wall.pos(cur);
            for (Pos d : kDir4) {
                Pos np = p + d;
                if (!wall.in_bounds(np) || wall[np] == 1) continue;
                int nid = wall.id(np);
                if (vis.try_mark(nid)) {
                    q.push(nid);
                }
            }
        }
    };

    bfs({0, 0});
}
```

---

## 各コンポーネント詳細リファレンス

### `FastSet`
- **用途**: 「現在空いているマス一覧からランダムに 1 つ選ぶ」「未配置のアイテムからランダムに選ぶ」など。
- **操作**:
  - `insert(x)` / `erase(x)` / `contains(x)`: すべて $O(1)$
  - `random_elem(rng)`: $O(1)$ で一様ランダムに 1 要素を取得
  - `data()`: 内部の密配列 `const std::vector<int>&` をそのままイテレート可能

### `FastHashSet`
- **注意**: `capacity` は **2 の冪乗** にし、格納予定数の **4 倍以上** を確保してください（負荷率 25% 以下で衝突をほぼゼロに抑えます）。
- `clear()` は世代カウンタを `+1` するだけなので、容量が $2^{20}$ でも $O(1)$ で一瞬で終わります。

### `Annealer`
- `check_mask`（デフォルト `1023` = $2^{10}-1$）: 何回に 1 回時刻取得と温度更新を行うかを制御します。
  - 1 反復が非常に重い問題（1 回あたり 1ms 近くかかる等）では `check_mask = 7` や `15` に下げると時間超過を防げます。
  - 1 反復が超軽量な問題ではデフォルトの `1023` や `4095` が最適です。

---

## AHC 本番で勝つためのチューニング Tips

1. **温度設定の目安 (`start_temp`, `end_temp`)**:
   - `start_temp`: 1 回の近傍遷移で起きるスコア変化の「大きめの値」程度（序盤の悪化受理率が 30〜60% になるように）。
   - `end_temp`: スコア変化の最小単位の $0.1$ 倍程度（終盤はほぼ山登りになるように）。
2. **`dump()` の活用**:
   - ローカルテスト時に `-DLOCAL` を付けて実行し、`Annealer::dump()` の `iter` と `accepted`、`NeighborSelector::dump()` の各近傍の `imp` を見るだけで、温度設定のミスや死んでいる近傍が一瞬で分かります。
3. **制限時間のマージン**:
   - 実行時間制限が 2000ms の問題なら、`limit_ms` は `1850`〜`1900` に設定し、出力処理や後処理の余裕を残すのが安全です。
