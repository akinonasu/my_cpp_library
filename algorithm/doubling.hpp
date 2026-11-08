#include <bit>
#include <concepts>
#include <print>
#include <utility>
#include <vector>

// 1. プライマリテンプレート: モノイド積対応ダブリング (重み付き)
// クラス名の直後の <S, Op> を除去し、基本定義とします
template <typename S = void, typename Op = std::nullptr_t>
class Doubling {
   private:
    int n;
    int log_k;
    std::vector<std::vector<int>> next;
    std::vector<std::vector<S>> data;
    Op op;
    S e;

   public:
    struct JumpResult {
        int to;
        S val;
    };

    struct MaxJumpResult {
        int to;
        S val;
        long step;
    };

    // n: 頂点数, max_k: 想定最大ステップ数 (64bit), op: モノイド演算, e: 単位元
    Doubling(int n, long max_k, Op op, S e) : n(n), log_k(max_k <= 0 ? 1 : std::bit_width(static_cast<unsigned long>(max_k))), next(log_k, std::vector<int>(n, -1)), data(log_k, std::vector<S>(n, e)), op(op), e(e) {}

    // 遷移元 from から遷移先 to への遷移と重みを設定
    void set_next(int from, int to, S weight) {
        next[0][from] = to;
        data[0][from] = weight;
    }

    // テーブルの構築: O(N log K)
    void build() {
        for (int b = 0; b + 1 < log_k; ++b) {
            for (int v = 0; v < n; ++v) {
                int mid = next[b][v];
                if (mid != -1) {
                    next[b + 1][v] = next[b][mid];
                    data[b + 1][v] = op(data[b][v], data[b][mid]);
                }
            }
        }
    }

    // start から step 回遷移した頂点と、通過した重みのモノイド積を返す: O(log K)
    JumpResult jump(int start, long step) const {
        int curr = start;
        S acc = e;
        for (int b = 0; b < log_k; ++b) {
            if ((step >> b) & 1) {
                if (curr == -1) break;
                acc = op(acc, data[b][curr]);
                curr = next[b][curr];
            }
        }
        return {curr, acc};
    }

    // pred(累積積) が true である間ジャンプし続けた結果を返す (二分探索)
    template <typename Pred>
    MaxJumpResult max_jump(int start, Pred pred) const {
        int curr = start;
        S acc = e;
        long step = 0;
        for (int b = log_k - 1; b >= 0; --b) {
            int nxt = next[b][curr];
            if (nxt != -1) {
                S next_acc = op(acc, data[b][curr]);
                if (pred(next_acc)) {
                    acc = next_acc;
                    curr = nxt;
                    step += (1L << b);
                }
            }
        }
        return {curr, acc, step};
    }
};

// 2. 完全特殊化: 重みなしダブリング (遷移先の追跡のみ)
template <>
class Doubling<void, std::nullptr_t> {
   private:
    int n;
    int log_k;
    std::vector<std::vector<int>> next;

   public:
    // n: 頂点数, max_k: 想定される最大ステップ数 (64bit)
    Doubling(int n, long max_k) : n(n), log_k(max_k <= 0 ? 1 : std::bit_width(static_cast<unsigned long>(max_k))), next(log_k, std::vector<int>(n, -1)) {}

    // 遷移元 from から遷移先 to への遷移を設定
    void set_next(int from, int to) { next[0][from] = to; }

    // テーブルの構築: O(N log K)
    void build() {
        for (int b = 0; b + 1 < log_k; ++b) {
            for (int v = 0; v < n; ++v) {
                if (next[b][v] != -1) {
                    next[b + 1][v] = next[b][next[b][v]];
                }
            }
        }
    }

    // start から step 回遷移した頂点を取得 (存在しない場合は -1): O(log K)
    int jump(int start, long step) const {
        int curr = start;
        for (int b = 0; b < log_k; ++b) {
            if ((step >> b) & 1) {
                if (curr == -1) break;
                curr = next[b][curr];
            }
        }
        return curr;
    }

    // pred(v) が true である間ジャンプし続けた先の頂点と累積ステップ数を返す
    template <typename Pred>
    std::pair<int, long> max_jump(int start, Pred pred) const {
        int curr = start;
        long step = 0;
        if (!pred(curr)) return {curr, 0};
        for (int b = log_k - 1; b >= 0; --b) {
            int nxt = next[b][curr];
            if (nxt != -1 && pred(nxt)) {
                curr = nxt;
                step += (1L << b);
            }
        }
        return {curr, step};
    }
};

// テンプレート引数推論補助 (CTAD)
Doubling(int, long) -> Doubling<void, std::nullptr_t>;

template <typename Op, typename S>
Doubling(int, long, Op, S) -> Doubling<S, Op>;
