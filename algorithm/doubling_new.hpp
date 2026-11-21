#pragma once

#include <bit>
#include <cassert>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

// ダブリング（頂点遷移のみ）
// 時間計算量: 構築 O(N log K), クエリ O(log K)
// 空間計算量: O(N log K)
// 事前条件:
//   - 0 <= max_k
//   - 各 i について nxt[i] == -1 または 0 <= nxt[i] < N
struct Doubling {
    Doubling() = default;

    Doubling(const std::vector<long>& nxt, long max_k) : n_(std::ssize(nxt)), max_k_(max_k), log_(max_k == 0 ? 0 : std::bit_width(static_cast<unsigned long>(max_k))), nxt_(log_ * n_, -1) {
        assert(0 <= max_k);
        if (log_ == 0) return;

        for (long i = 0; i < n_; ++i) {
            assert(nxt[i] == -1 || (0 <= nxt[i] && nxt[i] < n_));
            nxt_[i] = nxt[i];
        }
        for (long p = 0; p + 1 < log_; ++p) {
            long offset = p * n_;
            long next_offset = (p + 1) * n_;
            for (long i = 0; i < n_; ++i) {
                long mid = nxt_[offset + i];
                if (mid != -1) {
                    nxt_[next_offset + i] = nxt_[offset + mid];
                }
            }
        }
    }

    // 頂点 u から k 回移動した先の頂点を返す（到達不能なら -1）
    long jump(long u, long k) const {
        assert(0 <= u && u < n_);
        assert(0 <= k && k <= max_k_);
        for (long p = 0; p < log_; ++p) {
            if ((k >> p) & 1) {
                u = nxt_[p * n_ + u];
                if (u == -1) return -1;
            }
        }
        return u;
    }

   private:
    long n_ = 0;
    long max_k_ = 0;
    long log_ = 0;
    std::vector<long> nxt_;
};

// モノイド積付きダブリング
// 頂点 i から nxt[i] への 1 ステップの遷移に重み val[i] が対応する（非可換モノイド対応）
// 時間計算量: 構築 O(N log K), クエリ O(log K)
// 空間計算量: O(N log K)
// 事前条件:
//   - std::ssize(nxt) == std::ssize(val)
//   - 0 <= max_k
//   - 各 i について nxt[i] == -1 または 0 <= nxt[i] < N
//   - (S, op, e) はモノイドをなす（結合律: op(op(a, b), c) == op(a, op(b, c)), 単位元: e）
template <class S, class Op>
struct DoublingMonoid {
    struct StepResult {
        long step;
        long to;
        S prod;
    };

    DoublingMonoid(const std::vector<long>& nxt, const std::vector<S>& val, long max_k, Op op, std::type_identity_t<S> e) : n_(std::ssize(nxt)), max_k_(max_k), log_(max_k == 0 ? 0 : std::bit_width(static_cast<unsigned long>(max_k))), nxt_(log_ * n_, -1), prod_(log_ * n_, e), op_(op), e_(e) {
        assert(std::ssize(nxt) == std::ssize(val));
        assert(0 <= max_k);
        if (log_ == 0) return;

        for (long i = 0; i < n_; ++i) {
            assert(nxt[i] == -1 || (0 <= nxt[i] && nxt[i] < n_));
            nxt_[i] = nxt[i];
            if (nxt[i] != -1) {
                prod_[i] = val[i];
            }
        }
        for (long p = 0; p + 1 < log_; ++p) {
            long offset = p * n_;
            long next_offset = (p + 1) * n_;
            for (long i = 0; i < n_; ++i) {
                long mid = nxt_[offset + i];
                if (mid != -1) {
                    long to = nxt_[offset + mid];
                    nxt_[next_offset + i] = to;
                    if (to != -1) {
                        prod_[next_offset + i] = op_(prod_[offset + i], prod_[offset + mid]);
                    }
                }
            }
        }
    }

    // 頂点 u から k 回移動した先の {到達頂点, モノイド積} を返す
    // 途中で -1 に到達した場合は {-1, e} を返す
    std::pair<long, S> jump(long u, long k) const {
        assert(0 <= u && u < n_);
        assert(0 <= k && k <= max_k_);
        S res = e_;
        for (long p = 0; p < log_; ++p) {
            if ((k >> p) & 1) {
                long idx = p * n_ + u;
                if (nxt_[idx] == -1) return {-1, e_};
                res = op_(res, prod_[idx]);
                u = nxt_[idx];
            }
        }
        return {u, res};
    }

    // 頂点 u から 0 <= step <= max_k の範囲で、途中で -1 に落ちず、
    // かつモノイド積 s が pred(s) を満たす最大の {step, to, prod} を返す（単調性を仮定）
    // 事前条件: pred(e) == true, 0 <= max_k <= max_k_
    template <class Pred>
    StepResult max_step(long u, long max_k, Pred pred) const {
        assert(0 <= u && u < n_);
        assert(0 <= max_k && max_k <= max_k_);
        assert(pred(e_));
        long step = 0;
        S res = e_;
        for (long p = log_ - 1; p >= 0; --p) {
            long width = 1L << p;
            if (step <= max_k - width) {
                long idx = p * n_ + u;
                if (nxt_[idx] != -1) {
                    S nxt_res = op_(res, prod_[idx]);
                    if (pred(nxt_res)) {
                        res = std::move(nxt_res);
                        u = nxt_[idx];
                        step += width;
                    }
                }
            }
        }
        return {step, u, res};
    }

   private:
    long n_ = 0;
    long max_k_ = 0;
    long log_ = 0;
    std::vector<long> nxt_;
    std::vector<S> prod_;
    Op op_;
    S e_;
};
