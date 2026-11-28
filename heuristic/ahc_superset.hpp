#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <utility>
#include <vector>

// ============================================================================
// 1. Timer (高精度タイマー)
// 時間計算量: 各操作 O(1)
// 空間計算量: O(1)
// 事前条件: limit_ms > 0
// ============================================================================
struct Timer {
    Timer() : start_(std::chrono::steady_clock::now()) {}

    void reset() { start_ = std::chrono::steady_clock::now(); }

    long elapsed_us() const {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::microseconds>(now - start_).count();
    }

    long elapsed_ms() const { return elapsed_us() / 1000; }

    long double elapsed_sec() const { return static_cast<long double>(elapsed_us()) * 1e-6L; }

    long double progress(long limit_ms) const {
        assert(limit_ms > 0);
        long double p = static_cast<long double>(elapsed_us()) / (static_cast<long double>(limit_ms) * 1000.0L);
        return std::clamp(p, 0.0L, 1.0L);
    }

    bool is_over(long limit_ms) const { return elapsed_ms() >= limit_ms; }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "[Timer] elapsed=" << elapsed_ms() << "ms\n";
#else
        (void)os;
#endif
    }

   private:
    std::chrono::steady_clock::time_point start_;
};

// ============================================================================
// 2. XorShift (xoshiro256** 超高速・高品質乱数生成器)
// 時間計算量: 生成 O(1), shuffle O(N)
// 空間計算量: O(1)
// 事前条件: next_int(n) は n > 0, next_int(l, r) は l < r
// ============================================================================
struct XorShift {
    explicit XorShift(std::uint64_t seed = 88172645463325252ULL) {
        auto splitmix64 = [](std::uint64_t& x) {
            std::uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            return z ^ (z >> 31);
        };
        for (int i = 0; i < 4; ++i)
            s_[i] = splitmix64(seed);
    }

    std::uint64_t next_u64() {
        const std::uint64_t result = std::rotl(s_[1] * 5ULL, 7) * 9ULL;
        const std::uint64_t t = s_[1] << 17;
        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];
        s_[2] ^= t;
        s_[3] = std::rotl(s_[3], 45);
        return result;
    }

    int next_int(int n) {
        assert(n > 0);
        return static_cast<int>((static_cast<__uint128_t>(next_u64()) * static_cast<std::uint64_t>(n)) >> 64);
    }

    int next_int(int l, int r) {
        assert(l < r);
        return l + next_int(r - l);
    }

    long next_long(long n) {
        assert(n > 0);
        return static_cast<long>((static_cast<__uint128_t>(next_u64()) * static_cast<std::uint64_t>(n)) >> 64);
    }

    long next_long(long l, long r) {
        assert(l < r);
        return l + next_long(r - l);
    }

    long double next_real() {
        constexpr long double inv = 1.0L / 9007199254740992.0L;
        return static_cast<long double>(next_u64() >> 11) * inv;
    }

    long double next_real(long double l, long double r) {
        assert(l < r);
        return l + (r - l) * next_real();
    }

    bool bernoulli(long double p) { return next_real() < p; }

    template <class T>
    void shuffle(std::vector<T>& v) {
        long n = std::ssize(v);
        for (long i = n - 1; i > 0; --i) {
            std::swap(v[i], v[next_long(i + 1)]);
        }
    }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "[XorShift] state=(" << s_[0] << ", " << s_[1] << ", " << s_[2] << ", " << s_[3] << ")\n";
#else
        (void)os;
#endif
    }

   private:
    std::array<std::uint64_t, 4> s_{};
};

// ============================================================================
// 3. ZobristTable (O(1) 差分ハッシュ用テーブル)
// 時間計算量: 構築 O(cells * values), 参照 O(1)
// 空間計算量: O(cells * values)
// 事前条件: 0 <= cell < num_cells, 0 <= val < num_values
// ============================================================================
struct ZobristTable {
    ZobristTable() = default;
    ZobristTable(int num_cells, int num_values, XorShift& rng)
        : cells_(num_cells), values_(num_values), table_(static_cast<long>(num_cells) * num_values) {
        assert(num_cells >= 0 && num_values >= 0);
        for (auto& x : table_)
            x = rng.next_u64();
    }

    std::uint64_t operator()(int cell, int val) const {
        assert(0 <= cell && cell < cells_ && 0 <= val && val < values_);
        return table_[static_cast<long>(cell) * values_ + val];
    }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "[ZobristTable] cells=" << cells_ << " values=" << values_ << "\n";
#else
        (void)os;
#endif
    }

   private:
    int cells_ = 0;
    int values_ = 0;
    std::vector<std::uint64_t> table_;
};

// ============================================================================
// 4. StampArray (O(1) クリア可能な訪問管理配列: 反復 BFS/DFS 用)
// 時間計算量: 構築 O(N), clear / mark / is_marked / try_mark O(1)
// 空間計算量: O(N)
// 事前条件: 0 <= i < n
// ============================================================================
struct StampArray {
    explicit StampArray(long n = 0) : stamp_(n, 0) { assert(n >= 0); }

    void clear() {
        ++cur_;
        if (cur_ == 0) {
            std::fill(stamp_.begin(), stamp_.end(), 0);
            cur_ = 1;
        }
    }

    void mark(long i) {
        assert(0 <= i && i < std::ssize(stamp_));
        stamp_[i] = cur_;
    }

    bool is_marked(long i) const {
        assert(0 <= i && i < std::ssize(stamp_));
        return stamp_[i] == cur_;
    }

    bool try_mark(long i) {
        assert(0 <= i && i < std::ssize(stamp_));
        if (stamp_[i] == cur_) return false;
        stamp_[i] = cur_;
        return true;
    }

    long size() const { return std::ssize(stamp_); }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        long cnt = 0;
        for (int x : stamp_) {
            if (x == cur_) ++cnt;
        }
        os << "[StampArray] marked=" << cnt << " / " << std::ssize(stamp_) << "\n";
#else
        (void)os;
#endif
    }

   private:
    std::vector<int> stamp_;
    int cur_ = 1;
};

// ============================================================================
// 5. FastQueue (ゼロアロケーション固定容量リングキュー: 反復 BFS 用)
// 時間計算量: push / pop / clear O(1)
// 空間計算量: O(max_size)
// 事前条件: 同時格納数が capacity 未満
// ============================================================================
template <class T>
struct FastQueue {
    explicit FastQueue(long max_size = 1024) : buf_(max_size) { assert(max_size > 0); }

    void clear() {
        head_ = 0;
        tail_ = 0;
    }

    void push(const T& x) {
        assert(tail_ < std::ssize(buf_));
        buf_[tail_++] = x;
    }

    T pop() {
        assert(head_ < tail_);
        return buf_[head_++];
    }

    bool empty() const { return head_ == tail_; }
    long size() const { return tail_ - head_; }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "[FastQueue] size=" << (tail_ - head_) << " / cap=" << std::ssize(buf_) << "\n";
#else
        (void)os;
#endif
    }

   private:
    std::vector<T> buf_;
    long head_ = 0;
    long tail_ = 0;
};

// ============================================================================
// 6. FastHashSet (世代管理付きオープンアドレス法ハッシュセット: clear O(1))
// 時間計算量: insert / contains 期待 O(1), clear O(1)
// 空間計算量: O(capacity) (capacity は 2 の冪乗、格納予定数の 4 倍程度を推奨)
// 事前条件: capacity は 2 の冪乗 (>= 4)
// ============================================================================
struct FastHashSet {
    explicit FastHashSet(long capacity = 1 << 18) : mask_(capacity - 1), keys_(capacity, 0), gen_(capacity, 0) {
        assert(capacity >= 4 && (capacity & (capacity - 1)) == 0);
    }

    void clear() {
        ++cur_gen_;
        count_ = 0;
        if (cur_gen_ == 0) {
            std::fill(gen_.begin(), gen_.end(), 0);
            cur_gen_ = 1;
        }
    }

    bool insert(std::uint64_t key) {
        long idx = static_cast<long>(mix(key)) & mask_;
        while (gen_[idx] == cur_gen_) {
            if (keys_[idx] == key) return false;
            idx = (idx + 1) & mask_;
        }
        keys_[idx] = key;
        gen_[idx] = cur_gen_;
        ++count_;
        return true;
    }

    bool contains(std::uint64_t key) const {
        long idx = static_cast<long>(mix(key)) & mask_;
        while (gen_[idx] == cur_gen_) {
            if (keys_[idx] == key) return true;
            idx = (idx + 1) & mask_;
        }
        return false;
    }

    long size() const { return count_; }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "[FastHashSet] size=" << count_ << " / cap=" << (mask_ + 1) << "\n";
#else
        (void)os;
#endif
    }

   private:
    static std::uint64_t mix(std::uint64_t x) {
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27;
        return x;
    }

    long mask_;
    std::vector<std::uint64_t> keys_;
    std::vector<int> gen_;
    int cur_gen_ = 1;
    long count_ = 0;
};

// ============================================================================
// 7. FastSet ([0, n) の整数集合: 追加・削除・存在判定・ランダム取得すべて O(1))
// 時間計算量: 構築 O(N), insert / erase / contains / random_elem O(1)
// 空間計算量: O(N)
// 事前条件: 0 <= x < n
// ============================================================================
struct FastSet {
    explicit FastSet(long n = 0) : pos_(n, -1) {
        assert(n >= 0);
        elems_.reserve(n);
    }

    bool contains(int x) const {
        assert(0 <= x && x < std::ssize(pos_));
        return pos_[x] != -1;
    }

    bool insert(int x) {
        assert(0 <= x && x < std::ssize(pos_));
        if (pos_[x] != -1) return false;
        pos_[x] = static_cast<int>(std::ssize(elems_));
        elems_.push_back(x);
        return true;
    }

    bool erase(int x) {
        assert(0 <= x && x < std::ssize(pos_));
        int idx = pos_[x];
        if (idx == -1) return false;
        int last = elems_.back();
        elems_[idx] = last;
        pos_[last] = idx;
        elems_.pop_back();
        pos_[x] = -1;
        return true;
    }

    void clear() {
        for (int x : elems_)
            pos_[x] = -1;
        elems_.clear();
    }

    int random_elem(XorShift& rng) const {
        assert(!elems_.empty());
        return elems_[rng.next_long(std::ssize(elems_))];
    }

    long size() const { return std::ssize(elems_); }
    bool empty() const { return elems_.empty(); }
    const std::vector<int>& data() const { return elems_; }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "[FastSet] size=" << std::ssize(elems_) << " elems=[";
        for (long i = 0; i < std::ssize(elems_); ++i) {
            if (i > 0) os << ", ";
            os << elems_[i];
        }
        os << "]\n";
#else
        (void)os;
#endif
    }

   private:
    std::vector<int> elems_;
    std::vector<int> pos_;
};

// ============================================================================
// 8. NeighborSelector (焼きなまし用: 指数移動平均による近傍選択重みの自動適応)
// 時間計算量: select O(K), feedback O(1)
// 空間計算量: O(K)
// 事前条件: k >= 1
// ============================================================================
struct NeighborSelector {
    explicit NeighborSelector(int k, long double min_weight = 0.05L, long double decay = 0.995L)
        : k_(k), min_weight_(min_weight), decay_(decay), weight_(k, 1.0L), trials_(k, 0), gains_(k, 0) {
        assert(k >= 1);
    }

    int select(XorShift& rng) {
        long double sum = 0.0L;
        for (int i = 0; i < k_; ++i)
            sum += weight_[i];
        long double r = rng.next_real() * sum;
        for (int i = 0; i < k_; ++i) {
            r -= weight_[i];
            if (r <= 0.0L) {
                ++trials_[i];
                return i;
            }
        }
        ++trials_[k_ - 1];
        return k_ - 1;
    }

    void feedback(int id, bool improved) {
        assert(0 <= id && id < k_);
        if (improved) ++gains_[id];
        long double reward = improved ? 1.0L : 0.0L;
        weight_[id] = std::max(min_weight_, weight_[id] * decay_ + reward * (1.0L - decay_) * 10.0L);
    }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "[NeighborSelector]";
        for (int i = 0; i < k_; ++i) {
            os << " #" << i << "(w=" << weight_[i] << ", try=" << trials_[i] << ", imp=" << gains_[i] << ")";
        }
        os << "\n";
#else
        (void)os;
#endif
    }

   private:
    int k_;
    long double min_weight_;
    long double decay_;
    std::vector<long double> weight_;
    std::vector<long> trials_;
    std::vector<long> gains_;
};

// ============================================================================
// 9. Annealer (時間間引き・O(1) 対数テーブル・指数温度スケジュール 焼きなまし)
// 時間計算量: step / accept ならし O(1)
// 空間計算量: O(1) (テーブル 65536 要素)
// 事前条件: limit_ms > 0, start_temp > 0, end_temp > 0, check_mask は 2^k - 1
// ============================================================================
struct Annealer {
    static constexpr int kLogTableSize = 1 << 16;

    Annealer(long limit_ms, long double start_temp, long double end_temp, bool maximize = true, long check_mask = 1023)
        : limit_ms_(limit_ms),
          start_temp_(start_temp),
          end_temp_(end_temp),
          log_ratio_(std::log(end_temp / start_temp)),
          temp_(start_temp),
          maximize_(maximize),
          check_mask_(check_mask) {
        assert(limit_ms > 0);
        assert(start_temp > 0.0L && end_temp > 0.0L);
        assert((check_mask & (check_mask + 1)) == 0);
        for (int i = 0; i < kLogTableSize; ++i) {
            log_table_[i] = std::log((static_cast<long double>(i) + 0.5L) / static_cast<long double>(kLogTableSize));
        }
    }

    bool step() {
        if ((iter_ & check_mask_) == 0) {
            long double p = timer_.progress(limit_ms_);
            if (p >= 1.0L) return false;
            temp_ = start_temp_ * std::exp(log_ratio_ * p);
        }
        ++iter_;
        return true;
    }

    bool accept(long double delta, XorShift& rng) {
        long double gain = maximize_ ? delta : -delta;
        if (gain >= 0.0L) {
            ++accepted_;
            if (gain > 0.0L) ++improved_;
            return true;
        }
        long double threshold = temp_ * log_table_[rng.next_u64() & (kLogTableSize - 1)];
        if (gain >= threshold) {
            ++accepted_;
            return true;
        }
        return false;
    }

    long double temperature() const { return temp_; }
    long iterations() const { return iter_; }
    long accepted() const { return accepted_; }
    long improved() const { return improved_; }
    long elapsed_ms() const { return timer_.elapsed_ms(); }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "[Annealer] iter=" << iter_ << " accepted=" << accepted_ << " improved=" << improved_ << " temp=" << temp_
           << " elapsed=" << timer_.elapsed_ms() << "ms\n";
#else
        (void)os;
#endif
    }

   private:
    Timer timer_;
    std::array<long double, kLogTableSize> log_table_{};
    long limit_ms_;
    long double start_temp_;
    long double end_temp_;
    long double log_ratio_;
    long double temp_;
    bool maximize_;
    long check_mask_;
    long iter_ = 0;
    long accepted_ = 0;
    long improved_ = 0;
};

// ============================================================================
// 10. EulerTourBeam (状態コピーゼロ・差分更新型オイラーツアービームサーチ)
// 時間計算量: 1 ターンあたり O(W * (分岐数 + 差分更新コスト))
// 空間計算量: O(ターン数 * W)
// 事前条件: max_turns >= 1, beam_width >= 1
// ============================================================================
template <class State, class Op, class Score = long>
struct EulerTourBeam {
    struct Candidate {
        int parent = -1;
        Op op{};
        Score score{};
        std::uint64_t hash = 0;
    };

    EulerTourBeam(int beam_width, long hash_capacity = 1 << 18) : beam_width_(beam_width), seen_(hash_capacity) { assert(beam_width >= 1); }

    std::pair<Score, std::vector<Op>> run(State& state, int max_turns, long limit_ms = -1) {
        assert(max_turns >= 1);
        Timer timer;
        nodes_.clear();
        nodes_.push_back(Node{-1, 0, Op{}, Score{}});

        std::vector<int> leaves = {0};
        std::vector<Candidate> cands;
        cands.reserve(static_cast<long>(beam_width_) * 8);

        int cur_node = 0;

        auto move_to = [&](int target) {
            up_path_.clear();
            down_path_.clear();
            int u = cur_node;
            int v = target;
            while (u != v) {
                if (nodes_[u].depth > nodes_[v].depth) {
                    up_path_.push_back(u);
                    u = nodes_[u].parent;
                } else if (nodes_[v].depth > nodes_[u].depth) {
                    down_path_.push_back(v);
                    v = nodes_[v].parent;
                } else {
                    up_path_.push_back(u);
                    down_path_.push_back(v);
                    u = nodes_[u].parent;
                    v = nodes_[v].parent;
                }
            }
            for (int x : up_path_)
                state.undo(nodes_[x].op);
            for (long i = std::ssize(down_path_) - 1; i >= 0; --i) {
                state.apply(nodes_[down_path_[i]].op);
            }
            cur_node = target;
        };

        int best_node = 0;
        Score best_score = std::numeric_limits<Score>::lowest();

        for (int turn = 0; turn < max_turns; ++turn) {
            if (limit_ms > 0 && timer.is_over(limit_ms)) break;
            if (leaves.empty()) break;

            cands.clear();
            for (int leaf : leaves) {
                move_to(leaf);
                state.expand([&](const Op& op, Score score, std::uint64_t hash) { cands.push_back(Candidate{leaf, op, score, hash}); });
            }
            if (cands.empty()) break;

            auto comp = [](const Candidate& a, const Candidate& b) { return a.score > b.score; };
            if (std::ssize(cands) > static_cast<long>(beam_width_) * 2) {
                std::nth_element(cands.begin(), cands.begin() + beam_width_ * 2, cands.end(), comp);
                cands.resize(beam_width_ * 2);
            }
            std::sort(cands.begin(), cands.end(), comp);

            seen_.clear();
            std::vector<int> next_leaves;
            next_leaves.reserve(beam_width_);

            for (const auto& c : cands) {
                if (!seen_.insert(c.hash)) continue;
                int id = static_cast<int>(std::ssize(nodes_));
                nodes_.push_back(Node{c.parent, nodes_[c.parent].depth + 1, c.op, c.score});
                next_leaves.push_back(id);
                if (static_cast<int>(std::ssize(next_leaves)) == beam_width_) break;
            }

            if (next_leaves.empty()) break;
            std::sort(next_leaves.begin(), next_leaves.end(), [&](int a, int b) { return nodes_[a].parent < nodes_[b].parent; });

            for (int id : next_leaves) {
                if (nodes_[id].score > best_score) {
                    best_score = nodes_[id].score;
                    best_node = id;
                }
            }
            leaves = std::move(next_leaves);
            turns_done_ = turn + 1;
        }

        move_to(0);

        std::vector<Op> ops;
        for (int u = best_node; u != 0; u = nodes_[u].parent) {
            ops.push_back(nodes_[u].op);
        }
        std::reverse(ops.begin(), ops.end());
        best_score_ = best_score;
        return {best_score, ops};
    }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "[EulerTourBeam] turns=" << turns_done_ << " total_nodes=" << std::ssize(nodes_) << " best_score=" << best_score_ << "\n";
#else
        (void)os;
#endif
    }

   private:
    struct Node {
        int parent = -1;
        int depth = 0;
        Op op{};
        Score score{};
    };

    int beam_width_;
    FastHashSet seen_;
    std::vector<Node> nodes_;
    std::vector<int> up_path_;
    std::vector<int> down_path_;
    int turns_done_ = 0;
    Score best_score_{};
};

// ============================================================================
// 11. ChokudaiSearch (時間いっぱい自動でビーム幅を拡張する Chokudai サーチ)
//
// State が実装すべきインターフェース:
//   - Score eval() const
//   - std::uint64_t hash() const
//   - template <class F> void expand(F push_next) const
//       遷移先の (next_state, op) ごとに push_next(next_state, op) を呼ぶ
//
// 時間計算量: 制限時間 limit_ms まで反復
// 空間計算量: O(展開ノード総数)
// 事前条件: max_turns >= 1, width_per_iter >= 1, limit_ms > 0
// ============================================================================
template <class State, class Op, class Score = long>
struct ChokudaiSearch {
    explicit ChokudaiSearch(int width_per_iter = 1, long hash_capacity = 1 << 18) : width_(width_per_iter), seen_(hash_capacity) {
        assert(width_per_iter >= 1);
    }

    std::pair<Score, std::vector<Op>> run(const State& init_state, int max_turns, long limit_ms) {
        assert(max_turns >= 1 && limit_ms > 0);
        Timer timer;
        seen_.clear();
        nodes_.clear();

        auto cmp = [](const Entry& a, const Entry& b) { return a.score < b.score; };
        std::vector<std::priority_queue<Entry, std::vector<Entry>, decltype(cmp)>> beams(
            max_turns + 1, std::priority_queue<Entry, std::vector<Entry>, decltype(cmp)>(cmp));

        Score init_score = init_state.eval();
        nodes_.push_back(Node{-1, Op{}, init_state});
        seen_.insert(init_state.hash());
        beams[0].push(Entry{0, init_score});

        int best_node = 0;
        Score best_score = init_score;
        iters_ = 0;

        while (!timer.is_over(limit_ms)) {
            bool progressed = false;
            for (int t = 0; t < max_turns; ++t) {
                if (timer.is_over(limit_ms)) break;
                for (int w = 0; w < width_; ++w) {
                    if (beams[t].empty()) break;
                    Entry cur = beams[t].top();
                    beams[t].pop();
                    progressed = true;

                    int parent_id = cur.id;
                    // nodes_ の再確保に備えてコピーを取らず参照前に退避
                    State cur_state = nodes_[parent_id].state;
                    cur_state.expand([&](const State& nxt, const Op& op) {
                        if (!seen_.insert(nxt.hash())) return;
                        Score sc = nxt.eval();
                        int nid = static_cast<int>(std::ssize(nodes_));
                        nodes_.push_back(Node{parent_id, op, nxt});
                        beams[t + 1].push(Entry{nid, sc});
                        if (sc > best_score) {
                            best_score = sc;
                            best_node = nid;
                        }
                    });
                }
            }
            ++iters_;
            if (!progressed) break;
        }

        std::vector<Op> ops;
        for (int u = best_node; u != 0; u = nodes_[u].parent) {
            ops.push_back(nodes_[u].op);
        }
        std::reverse(ops.begin(), ops.end());
        best_score_ = best_score;
        return {best_score, ops};
    }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "[ChokudaiSearch] iters=" << iters_ << " nodes=" << std::ssize(nodes_) << " best_score=" << best_score_ << "\n";
#else
        (void)os;
#endif
    }

   private:
    struct Node {
        int parent = -1;
        Op op{};
        State state{};
    };
    struct Entry {
        int id = 0;
        Score score{};
    };

    int width_;
    FastHashSet seen_;
    std::vector<Node> nodes_;
    long iters_ = 0;
    Score best_score_{};
};

// ============================================================================
// 12. TourOptimizer (順列・巡回経路向け 2-opt / Or-opt 高速最適化)
// 時間計算量: 2-opt 1 パス O(N^2), Or-opt 1 パス O(N^2)
// 空間計算量: O(N)
// 事前条件: dist(u, v) は頂点 u -> v のコスト（long double 変換可能）
// ============================================================================
struct TourOptimizer {
    // 2-opt による局所最適化 (closed=true で閉路、false で始点・終点固定の開路)
    template <class DistFn>
    static void two_opt(std::vector<int>& tour, DistFn dist, bool closed = true, long max_iters = 100000) {
        long n = std::ssize(tour);
        if (n < 4) return;
        for (long iter = 0; iter < max_iters; ++iter) {
            bool improved = false;
            long lim_i = closed ? n - 1 : n - 2;
            for (long i = 0; i < lim_i; ++i) {
                int a = tour[i];
                int b = tour[i + 1];
                long lim_j = closed ? (i == 0 ? n - 1 : n) : n - 1;
                for (long j = i + 2; j < lim_j; ++j) {
                    int c = tour[j];
                    int d = tour[(j + 1) % n];
                    if (dist(a, c) + dist(b, d) + 1e-12L < dist(a, b) + dist(c, d)) {
                        std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                        improved = true;
                        b = tour[i + 1];
                    }
                }
            }
            if (!improved) break;
        }
    }

    // Or-opt (長さ 1..max_seg の連続区間を別位置に挿入) による局所最適化
    template <class DistFn>
    static void or_opt(std::vector<int>& tour, DistFn dist, bool closed = true, int max_seg = 3, long max_iters = 10000) {
        long n = std::ssize(tour);
        if (n < 4) return;
        for (long iter = 0; iter < max_iters; ++iter) {
            bool improved = false;
            for (int seg = 1; seg <= max_seg && seg < n - 1; ++seg) {
                long start_l = closed ? 0 : 1;
                long end_l = closed ? n - seg : n - seg - 1;
                for (long l = start_l; l <= end_l; ++l) {
                    long r = l + seg - 1;
                    int prev_l = tour[(l - 1 + n) % n];
                    int first = tour[l];
                    int last = tour[r];
                    int next_r = tour[(r + 1) % n];
                    long double remove_gain = (dist(prev_l, first) + dist(last, next_r)) - dist(prev_l, next_r);

                    long best_pos = -1;
                    long double best_delta = -1e-12L;
                    long ins_end = closed ? n : n - 1;
                    for (long k = 0; k < ins_end; ++k) {
                        if (k >= l - 1 && k <= r) continue;
                        int u = tour[k];
                        int v = tour[(k + 1) % n];
                        long double add_cost = (dist(u, first) + dist(last, v)) - dist(u, v);
                        long double delta = add_cost - remove_gain;
                        if (delta < best_delta) {
                            best_delta = delta;
                            best_pos = k;
                        }
                    }
                    if (best_pos != -1) {
                        std::vector<int> block(tour.begin() + l, tour.begin() + r + 1);
                        tour.erase(tour.begin() + l, tour.begin() + r + 1);
                        long insert_at = (best_pos < l) ? (best_pos + 1) : (best_pos - seg + 1);
                        tour.insert(tour.begin() + insert_at, block.begin(), block.end());
                        improved = true;
                        break;
                    }
                }
                if (improved) break;
            }
            if (!improved) break;
        }
    }

    template <class DistFn>
    static long double total_cost(const std::vector<int>& tour, DistFn dist, bool closed = true) {
        long n = std::ssize(tour);
        if (n <= 1) return 0.0L;
        long double sum = 0.0L;
        for (long i = 0; i + 1 < n; ++i)
            sum += dist(tour[i], tour[i + 1]);
        if (closed) sum += dist(tour[n - 1], tour[0]);
        return sum;
    }
};

// ============================================================================
// 13. Pos & Grid (キャッシュ効率重視の 1 次元フラット 2D グリッド)
// 時間計算量: 構築 O(HW), 要素アクセス O(1)
// 空間計算量: O(HW)
// 事前条件: h >= 0, w >= 0, アクセス時は 0 <= r < h, 0 <= c < w
// ============================================================================
struct Pos {
    int r = 0;
    int c = 0;

    Pos operator+(const Pos& o) const { return {r + o.r, c + o.c}; }
    Pos operator-(const Pos& o) const { return {r - o.r, c - o.c}; }
    bool operator==(const Pos& o) const = default;

    int manhattan(const Pos& o) const { return std::abs(r - o.r) + std::abs(c - o.c); }
};

inline constexpr std::array<Pos, 4> kDir4 = {{{-1, 0}, {0, 1}, {1, 0}, {0, -1}}};
inline constexpr std::array<char, 4> kDir4Char = {'U', 'R', 'D', 'L'};
inline constexpr std::array<Pos, 8> kDir8 = {{{-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}}};

template <class T>
struct Grid {
    Grid() = default;
    Grid(int h, int w, const T& init = T{}) : h_(h), w_(w), data_(static_cast<long>(h) * w, init) { assert(h >= 0 && w >= 0); }

    bool in_bounds(int r, int c) const { return 0 <= r && r < h_ && 0 <= c && c < w_; }
    bool in_bounds(Pos p) const { return in_bounds(p.r, p.c); }
    int id(int r, int c) const {
        assert(in_bounds(r, c));
        return r * w_ + c;
    }
    int id(Pos p) const { return id(p.r, p.c); }
    Pos pos(int idx) const {
        assert(0 <= idx && idx < h_ * w_);
        return {idx / w_, idx % w_};
    }

    T& operator()(int r, int c) {
        assert(in_bounds(r, c));
        return data_[static_cast<long>(r) * w_ + c];
    }
    const T& operator()(int r, int c) const {
        assert(in_bounds(r, c));
        return data_[static_cast<long>(r) * w_ + c];
    }
    T& operator[](Pos p) { return (*this)(p.r, p.c); }
    const T& operator[](Pos p) const { return (*this)(p.r, p.c); }

    void fill(const T& val) { std::fill(data_.begin(), data_.end(), val); }
    int height() const { return h_; }
    int width() const { return w_; }
    long size() const { return std::ssize(data_); }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "[Grid] " << h_ << "x" << w_ << "\n";
        for (int r = 0; r < h_; ++r) {
            for (int c = 0; c < w_; ++c) {
                if (c > 0) os << ' ';
                os << (*this)(r, c);
            }
            os << '\n';
        }
#else
        (void)os;
#endif
    }

   private:
    int h_ = 0;
    int w_ = 0;
    std::vector<T> data_;
};
