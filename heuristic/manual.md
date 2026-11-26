# ヒューリスティックコンテスト用 C++23 ライブラリ (`heuristic.hpp`) 実装マニュアル

定数倍高速化に特化した、1ファイル完結の C++23 ヒューリスティック用ライブラリ（`Timer` / `Random` / `FastHashSet` / `Annealer` / `BeamSearch` / `EulerBeamSearch`）の実践マニュアルです。

---

## 目次

1. [環境・ビルド方法](#1-環境ビルド方法)
2. [手法選択ガイドライン](#2-手法選択ガイドライン)
3. [API クイックリファレンス](#3-api-クイックリファレンス)
4. [共通部品 (`Timer`, `Random`, `FastHashSet`)](#4-共通部品-timer-random-fasthashset)
5. [焼きなまし法 (`Annealer`) 実装マニュアル](#5-焼きなまし法-annealer-実装マニュアル)
6. [標準ビームサーチ (`BeamSearch`) 実装マニュアル](#6-標準ビームサーチ-beamsearch-実装マニュアル)
7. [ロールバック型ビームサーチ (`EulerBeamSearch`) 実装マニュアル](#7-ロールバック型ビームサーチ-eulerbeamsearch-実装マニュアル)
8. [実戦サンプルコード集](#8-実戦サンプルコード集)
9. [定数倍高速化チェックリスト](#9-定数倍高速化チェックリスト)

---

## 1. 環境・ビルド方法

- **対応規格**: C++23 (`-std=c++23`)
- **規約**: すべて **0-indexed**・**半開区間 `[l, r)`**・浮動小数点数は `long double`・添字とサイズは `long` で統一しています。
- **本番ゼロオーバーヘッド**: すべての構造体の `dump()` メンバ関数と焼きなまし内部カウンタは `-DLOCAL` が定義されているときだけ有効になります。本番ジャッジ時は `-DLOCAL` を外すことで最適化により完全に消滅します。

```bash
# ローカルテスト時（アサート・dump デバッグ出力が有効）
g++ -std=c++23 -O2 -Wall -Wextra -DLOCAL main.cpp -o main

# 本番提出時（dump 等はゼロオーバーヘッドで消滅）
g++ -std=c++23 -O2 -Wall -Wextra main.cpp -o main
```

---

## 2. 手法選択ガイドライン

| 手法 | 適している問題 | 使用クラス |
| :--- | :--- | :--- |
| **焼きなまし法** | 経路・配置・割当などの組合せ最適化。局所的な変更（近傍）と差分スコア計算が作れる。 | `Annealer` |
| **標準ビームサーチ** | 手順を先頭から決める問題で、**状態サイズが小さい**（数バイト〜数百バイト程度）。`undo` を書かずに高速実装したい。 | `BeamSearch` |
| **ロールバック型ビームサーチ** | 手順を先頭から決める問題で、**状態サイズが大きい**（盤面・グラフ・複数配列など）。コピーをゼロにしてビーム幅を稼ぎたい。 | `EulerBeamSearch` |

---

## 3. API クイックリファレンス

### `Timer`
| メソッド | 説明 | 計算量 |
| :--- | :--- | :--- |
| `reset()` | 計測開始時刻を現在時刻にリセット | $O(1)$ |
| `elapsed_ns()` / `elapsed_us()` / `elapsed_ms()` | 経過時間を整数 (`long`) で取得 | $O(1)$ |
| `elapsed_sec()` | 経過時間を秒 (`long double`) で取得 | $O(1)$ |
| `progress(long limit_ms)` | `limit_ms` に対する進捗率 (`0.0L` 〜 `1.0L` 以上) | $O(1)$ |
| `is_over(long limit_ms)` | 制限時間を超過したか判定 | $O(1)$ |
| `dump(os = std::cerr)` | 経過時間を出力 (`-DLOCAL` 時のみ) | $O(1)$ |

### `Random` (`xoroshiro128++`)
| メソッド | 説明 | 計算量 |
| :--- | :--- | :--- |
| `Random(uint64_t seed = 42)` / `set_seed(seed)` | シード設定（再現性のため固定値推奨） | $O(1)$ |
| `next_u64()` / `next_u32()` | 64bit / 32bit 符号なし一様乱数 | $O(1)$ |
| `next_int(int n)` / `next_int(int l, int r)` | `[0, n)` / `[l, r)` の `int` 乱数（除算なし） | $O(1)$ |
| `next_long(long n)` / `next_long(long l, long r)` | `[0, n)` / `[l, r)` の `long` 乱数（`__uint128_t` 乗算） | $O(1)$ |
| `next_real()` / `next_real(l, r)` | `[0, 1)` / `[l, r)` の `long double` 乱数 | $O(1)$ |
| `next_bool()` / `next_bool(long double p)` | 50% / 確率 `p` で `true` | $O(1)$ |
| `shuffle(std::vector<T>& v)` | Fisher-Yates 法で高速シャッフル | $O(N)$ |

### `FastHashSet`
| メソッド | 説明 | 計算量 |
| :--- | :--- | :--- |
| `FastHashSet(int init_bits = 10)` | 初期容量 `2^init_bits` で構築（負荷率 70% 超で自動拡張） | $O(2^{\text{init\_bits}})$ |
| `clear()` | 世代番号更新による高速クリア | **$O(1)$** |
| `insert(uint64_t key)` | 未登録なら挿入し `true`、既出なら `false` | 平均 $O(1)$ |
| `contains(uint64_t key)` | 登録済みか判定 | 平均 $O(1)$ |

### `Annealer`
| メンバ / メソッド | 説明 |
| :--- | :--- |
| `rng` / `timer` | 公開メンバ。`sa.rng.next_int(...)` のように直接利用可能 |
| `Annealer(start_t, end_t, limit_ms, schedule = Exp, seed = 42)` | 内部新規タイマーで構築 |
| `Annealer(start_t, end_t, limit_ms, base_timer, schedule = Exp, seed = 42)` | 既存 `Timer` を引き継いで構築（前処理時間を自動考慮） |
| `Annealer::by_iters(start_t, end_t, max_iters, schedule = Exp, seed = 42)` | 反復回数固定モード（テスト・ベンチマーク向け） |
| `step()` | 継続判定＋必要時の温度更新（適応的間引きにより超軽量） |
| `accept(delta)` / `accept(new_s, cur_s)` | 最大化問題の遷移判定（`accept_max` と同じ） |
| `accept_min(delta)` / `accept_min(new_s, cur_s)` | 最小化問題の遷移判定 |
| `temp()` / `progress()` / `iters()` / `elapsed_ms()` | 現在温度・進捗 `[0, 1]`・正確な反復回数・経過時間 |
| `dump(os = std::cerr)` | 反復回数・改善数・悪化採択数・採択率を出力 (`-DLOCAL`) |

### `BeamSearch<State, Action, Score = long, Maximize = true>`
| メソッド | 説明 |
| :--- | :--- |
| `BeamSearch(init_state, init_score = Score{})` / `reset(...)` | 初期状態と初期スコアで初期化 |
| `reserve_cands(long n)` | 候補バッファの事前予約（`ビーム幅 * 分岐数` 推奨） |
| `size()` / `turn()` / `state(i)` / `score(i)` | 現在のビーム数・ターン・状態参照・スコア参照 |
| `push(score, parent, action)` | 重複統合なしで候補追加 |
| `push(score, parent, action, uint64_t hash)` | `hash` で重複を統合（良スコアのみ残す） |
| `advance(next_w, apply, sort_beam = true)` | `std::nth_element` で上位 `next_w` 個を選び、生存候補のみ `apply(State&, const Action&)` を実行 |
| `run(max_turns, w, expand, apply, sort_beam = true)` | `max_turns` 分を一括実行（`expand(long parent, const State& s)`） |
| `best_index()` / `best_state()` / `best_score()` | 現在のビーム内の最良解を取得 |
| `restore_actions(idx = -1)` | 探索木から最良解（または `idx`）の手順を復元 |

### `EulerBeamSearch<Action, Score = long, Maximize = true>`
| メソッド | 説明 |
| :--- | :--- |
| `EulerBeamSearch()` / `reset()` / `reserve_cands(n)` | 初期化・バッファ予約 |
| `turn()` / `size()` | 現在ターン数・現在の葉（生存状態）数 |
| `push(score, action)` / `push(score, action, uint64_t hash)` | `expand()` 内で候補を追加（親ノードは自動管理） |
| `record_best(score)` | `expand()` 内で現在状態の解を全ターン中ベスト候補として記録 |
| `step(w, apply, undo, expand, track_overall_best = false)` | 1ターン進める（`apply(const Action&)`, `undo(const Action&)`, `expand()`） |
| `run(max_turns, w, apply, undo, expand, track_overall_best = false)` | `max_turns` 分を一括実行 |
| `has_best()` / `best_score()` / `best_actions()` | 最終生き残り中の最良スコアと手順を取得 |
| `has_overall_best()` / `overall_best_score()` / `overall_best_actions()` | 途中ターン含む全体の最良スコアと手順を取得 |

---

## 4. 共通部品 (`Timer`, `Random`, `FastHashSet`)

### 4.1 `main()` 先頭での時間・乱数・Zobrist Hash 準備

```cpp
int main() {
    Timer global_timer; // プログラム開始直後に宣言する
    Random rng(42);     // シード固定（再現性確保）

    // 入力読み込み・前処理...

    // Zobrist Hash 用テーブルの例（盤面 HxW × 種類数 K）
    std::vector<std::uint64_t> zobrist(H * W * K);
    for (auto& h : zobrist) {
        h = rng.next_u64();
    }
}
```

### 4.2 任意の構造体・配列から 64bit ハッシュを作るコツ

`splitmix64(x)` を組み合わせることで、複数要素のハッシュ衝突を極めて低く抑えられます。

```cpp
std::uint64_t h = 0;
for (int i = 0; i < n; ++i) {
    h ^= splitmix64(static_cast<std::uint64_t>(i) * 1000003ULL + static_cast<std::uint64_t>(arr[i]));
}
```

---

## 5. 焼きなまし法 (`Annealer`) 実装マニュアル

### Step 1: 4要素の設計
1. **状態**: 現在の解を表す最小限のデータ構造
2. **初期解**: 貪欲法またはランダム構築
3. **近傍**: 1点変更・2点スワップ・区間反転など
4. **差分評価 (`delta`)**: 状態全体を再計算せず、変更に関わる要素だけを再計算する

### Step 2: 実装テンプレート

#### パターン A: 変更前に `delta` を計算できる場合（推奨・最速）
棄却時のロールバック処理が一切不要になります。

```cpp
// 前処理の経過時間を引き継いで合計 1950ms で終了させる
Annealer sa(start_temp, end_temp, 1950, global_timer);

long cur_score = calc_full_score();
long best_score = cur_score;
auto best_state = cur_state; // ※最高記録は必ず別変数に退避する

while (sa.step()) {
    int i = sa.rng.next_int(n);
    int new_val = sa.rng.next_int(m);

    long delta = calc_delta(i, new_val);

    // 最大化: accept(delta) / accept_max(delta)
    // 最小化: accept_min(delta)
    if (sa.accept(delta)) {
        apply_move(i, new_val);
        cur_score += delta;

        if (cur_score > best_score) {
            best_score = cur_score;
            best_state = cur_state;
        }
    }
}
sa.dump(); // ローカル時のみ出力
```

#### パターン B: 一旦変更を適用しないと評価できない場合
変更して新スコアを計算し、棄却されたら元の状態へ戻します。

```cpp
while (sa.step()) {
    int i = sa.rng.next_int(n);
    int old_val = cur_state[i]; // 復元用データを控えておく
    int new_val = sa.rng.next_int(m);

    apply_move(i, new_val);
    long new_score = eval_score();

    // 2引数版 accept(new_score, cur_score) はオーバーフローしない
    if (sa.accept(new_score, cur_score)) {
        cur_score = new_score;
        if (cur_score > best_score) {
            best_score = cur_score;
            best_state = cur_state;
        }
    } else {
        apply_move(i, old_val); // ロールバック
    }
}
```

### Step 3: 温度の決め方と `sa.dump()` の見方

ローカル実行すると以下のようなログが出ます。

```text
Annealer: iters=4521984, elapsed=1950 ms, temp=0.1, progress=1, improved=12054, bad_accept=215032, reject=4294898, accept_rate=5.02%
```

1. **開始温度 `start_temp`**: 序盤の典型的な悪化幅 $\Delta_{\text{typ}}$ に対して受理確率 $\exp(-\Delta_{\text{typ}} / T_{\text{start}}) \approx 0.3 \sim 0.5$ になるよう $T_{\text{start}} \approx \Delta_{\text{typ}}$ 前後から試します。
   - `bad_accept` が 0 付近なら温度が低すぎます（ただの山登り法になっています）。
2. **終了温度 `end_temp`**: 終盤に最小の悪化幅 $\Delta_{\min}$（整数スコアなら通常 1）をほぼ通さない値（$T_{\text{end}} \approx 0.1 \times \Delta_{\min}$ 程度）にします。
   - 終盤までスコアが収束せず上下に暴れるなら `end_temp` を下げます。
3. **高度なテクニック**:
   - **近傍の動的切り替え**: `if (sa.rng.next_real() < 0.5L * (1.0L - sa.progress()))` のように、序盤は大きな変更、終盤は局所微調整の比率を高めるとスコアが伸びやすくなります。
   - **実験用固定回数モード**: `Annealer::by_iters(start_temp, end_temp, 1000000)` を使うと、実行環境に左右されず毎回同一の反復回数で比較できます。

---

## 6. 標準ビームサーチ (`BeamSearch`) 実装マニュアル

### Step 1: 設計の鉄則
- **絶対に `State` の中に操作履歴 (`std::vector<Action>`) を持たせないでください。**  
  `BeamSearch` が内部の木で自動管理し、`restore_actions()` で最後に 1 回だけ復元します。
- 重複状態が出うる問題では、`push(score, parent, action, hash)` の第4引数に **Zobrist Hash** を渡してください。ハッシュ重複が自動的に最良スコアの1つへ統合されます。

### Step 2: 実装テンプレート

```cpp
BeamSearch<State, Action, long> bs(init_state, init_score);
bs.reserve_cands(beam_width * branch_factor);

for (long t = 0; t < max_turns; ++t) {
    for (long i = 0; i < bs.size(); ++i) {
        const State& s = bs.state(i);
        for (const Action& act : enum_actions(s)) {
            long next_score = calc_next_score(s, act);
            std::uint64_t next_hash = calc_next_hash(s, act);

            // 重複排除しない場合は hash を省略して bs.push(next_score, i, act)
            bs.push(next_score, i, act, next_hash);
        }
    }

    // 残り時間に応じた動的ビーム幅調整も可能
    long cur_w = global_timer.elapsed_ms() > 1700 ? 50 : beam_width;

    // 上位 cur_w 個に絞ってから、生き残り分だけ状態を複製・適用する
    if (!bs.advance(cur_w, [&](State& s, const Action& act) { s.apply(act); })) {
        break;
    }
}

// 結果の取り出し
const State& best = bs.best_state();
std::vector<Action> actions = bs.restore_actions();
```

- 手動制御が不要な場合は `bs.run(max_turns, beam_width, expand_fn, apply_fn)` で一括実行できます。
- 毎ターン上位ビームをソートする必要すらない場合、`advance(cur_w, apply, false)` とするとソートコストを省けます（`best_state()` や `restore_actions()` はソートなしでも正しく最大/最小値を検索します）。

---

## 7. ロールバック型ビームサーチ (`EulerBeamSearch`) 実装マニュアル

盤面が大きい問題（コピーコストが重い場合）では、`State` を **1 つだけ** 持ち、探索木上を `apply` と `undo` で移動する `EulerBeamSearch` を使います。

### Step 1: `Action` に「ロールバック情報」を含める
`undo(act)` で完全に元の状態に戻せるよう、変更前の値などを `Action` に持たせるのが定石です。

```cpp
struct Action {
    int r = 0;
    int c = 0;
    int new_val = 0;
    int old_val = 0;    // undo で元に戻すための値
    long old_score = 0; // スコア差分が複雑なら変更前スコアを保持
};
```
※ `Action` はデフォルト構築可能（初期値付き）である必要があります。

### Step 2: 実装テンプレート

```cpp
State state = init_state; // メモリ上の状態はこれ 1 つだけ

auto apply = [&](const Action& act) {
    state.board[act.r][act.c] = act.new_val;
    state.score = calc_new_score(...);
    state.hash ^= zobrist[act.r][act.c][act.old_val]
                ^ zobrist[act.r][act.c][act.new_val];
};

auto undo = [&](const Action& act) {
    state.hash ^= zobrist[act.r][act.c][act.new_val]
                ^ zobrist[act.r][act.c][act.old_val];
    state.board[act.r][act.c] = act.old_val;
    state.score = act.old_score;
};

EulerBeamSearch<Action, long> bs;
bs.reserve_cands(beam_width * branch_factor);

for (long t = 0; t < max_turns; ++t) {
    long cur_w = global_timer.elapsed_ms() > 1700 ? 100 : beam_width;

    bool ok = bs.step(cur_w, apply, undo, [&]() {
        // この中は各生き残り葉の状態に遷移した状態で呼ばれる
        for (...) {
            Action act = make_action(state, ...);
            long next_score = ...;
            std::uint64_t next_hash = ...;

            // 親ノード ID は自動管理されるため渡す必要なし
            bs.push(next_score, act, next_hash);
        }
    });
    if (!ok) {
        break;
    }
}

std::vector<Action> best_actions = bs.best_actions();
```

### Step 3: 途中ターンの最高点を採用したい場合
- **毎ターンのトップを自動保存**：`bs.step(..., /*track_overall_best=*/true)` または `bs.run(..., /*track_overall_best=*/true)` を指定し、`bs.overall_best_actions()` で取得します。
- **任意条件（ゴール到達など）で保存**：`expand` 内で `if (is_goal(state)) { bs.record_best(state.score); }` を呼び、`bs.overall_best_actions()` で取得します。

### Step 4: `apply` / `undo` のバグ防止アサート
`EulerBeamSearch` がおかしい場合、99% は `undo` が完全な逆操作になっていません。開発時は以下で検証してください。

```cpp
#ifdef LOCAL
    State snapshot = state;
    apply(act);
    undo(act);
    assert(state.score == snapshot.score);
    assert(state.hash == snapshot.hash);
    assert(state.board == snapshot.board);
#endif
```

---

## 8. 実戦サンプルコード集

### 8.1 TSP (2-opt 焼きなまし)

```cpp
#include "heuristic.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

int main() {
    constexpr int n = 100;
    Random init_rng(12345);

    std::vector<long double> x(n), y(n);
    for (int i = 0; i < n; ++i) {
        x[i] = init_rng.next_real(0.0L, 1000.0L);
        y[i] = init_rng.next_real(0.0L, 1000.0L);
    }

    std::vector<long double> dist(static_cast<std::size_t>(n) * n);
    auto d = [&](int i, int j) { return dist[static_cast<long>(i) * n + j]; };
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            long double dx = x[i] - x[j], dy = y[i] - y[j];
            dist[static_cast<long>(i) * n + j] = std::sqrt(dx * dx + dy * dy);
        }
    }

    std::vector<int> order(n);
    for (int i = 0; i < n; ++i) order[i] = i;
    long double cur_len = 0.0L;
    for (int i = 0; i < n; ++i) cur_len += d(order[i], order[(i + 1) % n]);

    long double best_len = cur_len;
    std::vector<int> best_order = order;

    Annealer sa(50.0L, 0.1L, 500);
    while (sa.step()) {
        int l = sa.rng.next_int(1, n - 1);
        int r = sa.rng.next_int(l + 2, n + 1);
        if (l == 1 && r == n) continue;

        int a = order[l - 1], b = order[l];
        int c = order[r - 1], e = order[r == n ? 0 : r];
        long double delta = (d(a, c) + d(b, e)) - (d(a, b) + d(c, e));

        if (sa.accept_min(delta)) {
            std::reverse(order.begin() + l, order.begin() + r);
            cur_len += delta;
            if (cur_len < best_len) {
                best_len = cur_len;
                best_order = order;
            }
        }
    }
    sa.dump();
    std::cout << "best_len = " << best_len << "\n";
}
```

### 8.2 グリッドコイン回収 (`EulerBeamSearch` 版)

```cpp
#include "heuristic.hpp"
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

constexpr int H = 15, W = 15;
constexpr long TURNS = 40, BEAM_WIDTH = 300;
constexpr int DR[4] = {-1, 1, 0, 0}; // dir ^ 1 が逆方向
constexpr int DC[4] = {0, 0, -1, 1};

int main() {
    Random rng(999);
    std::array<long, H * W> coin{};
    std::array<std::uint64_t, H * W> z_coin{}, z_pos{};
    for (int i = 0; i < H * W; ++i) {
        coin[i] = rng.next_long(1, 101);
        z_coin[i] = rng.next_u64();
        z_pos[i] = rng.next_u64();
    }

    struct SingleState {
        int r = 0, c = 0;
        long score = 0;
        std::uint64_t hash = 0;
        std::array<int, H * W> visit_cnt{};
    } state;

    state.visit_cnt[0] = 1;
    state.score = coin[0];
    state.hash = z_pos[0] ^ z_coin[0];

    auto apply = [&](int dir) {
        int oidx = state.r * W + state.c;
        state.r += DR[dir];
        state.c += DC[dir];
        int nidx = state.r * W + state.c;
        state.hash ^= z_pos[oidx] ^ z_pos[nidx];
        if (state.visit_cnt[nidx]++ == 0) {
            state.score += coin[nidx];
            state.hash ^= z_coin[nidx];
        }
    };

    auto undo = [&](int dir) {
        int oidx = state.r * W + state.c;
        if (--state.visit_cnt[oidx] == 0) {
            state.score -= coin[oidx];
            state.hash ^= z_coin[oidx];
        }
        int rev = dir ^ 1;
        state.r += DR[rev];
        state.c += DC[rev];
        int nidx = state.r * W + state.c;
        state.hash ^= z_pos[oidx] ^ z_pos[nidx];
    };

    EulerBeamSearch<int, long> bs;
    bs.reserve_cands(BEAM_WIDTH * 4);
    bs.run(TURNS, BEAM_WIDTH, apply, undo, [&]() {
        int oidx = state.r * W + state.c;
        for (int dir = 0; dir < 4; ++dir) {
            int nr = state.r + DR[dir], nc = state.c + DC[dir];
            if (nr < 0 || nr >= H || nc < 0 || nc >= W) continue;
            int nidx = nr * W + nc;
            long nscore = state.score + (state.visit_cnt[nidx] == 0 ? coin[nidx] : 0);
            std::uint64_t nhash = state.hash ^ z_pos[oidx] ^ z_pos[nidx]
                                             ^ (state.visit_cnt[nidx] == 0 ? z_coin[nidx] : 0ULL);
            bs.push(nscore, dir, nhash);
        }
    });

    bs.dump();
    std::vector<int> ans = bs.best_actions();
    std::cout << "score = " << bs.best_score() << ", steps = " << std::ssize(ans) << "\n";
}
```V

---

## 9. 定数倍高速化チェックリスト

- [ ] **差分計算**: 焼きなまし・ビームサーチともに、変更箇所だけの差分評価を行う（全盤面走査を避ける）。
- [ ] **評価と適用の分離**:
  - `Annealer`: 可能な限り `accept(delta)` で採用判定してから状態を変更し、棄却時は何もしない。
  - `BeamSearch`: `push` の時点で `State` をコピーせず、`advance` で上位 $W$ 個に絞ってから更新する。
- [ ] **`std::unordered_set` を使わない**: 64bit ハッシュの重複判定には `FastHashSet`（または `BeamSearch::push(..., hash)`）を使う。
- [ ] **メモリ事前確保**: ビームサーチでは `reserve_cands(beam_width * branch_factor)` を最初に呼ぶ。
- [ ] **タイマーの引き継ぎ**: 初期解作成や前処理が長いときは `Annealer(..., global_timer)` を使い、TLE を防ぐ。
- [ ] **`-DLOCAL` なし提出**: 提出用ビルド／提出コードでは `-DLOCAL` を定義せず、`dump()` やデバッグ集計のコストを完全にゼロにする。V