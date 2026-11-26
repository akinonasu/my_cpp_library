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
#include <utility>
#include <vector>

// 64bit 整数ハッシュ混合関数 (SplitMix64)
// 計算量: 時間 O(1), 空間 O(1)
inline std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// 経過時間計測タイマー
// 計算量: 各操作 時間 O(1), 空間 O(1)
struct Timer {
    Timer() : start_(std::chrono::steady_clock::now()) {}

    void reset() { start_ = std::chrono::steady_clock::now(); }

    long elapsed_ns() const {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::nanoseconds>(now - start_).count();
    }

    long elapsed_us() const {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::microseconds>(now - start_).count();
    }

    long elapsed_ms() const {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(now - start_).count();
    }

    long double elapsed_sec() const { return static_cast<long double>(elapsed_ns()) * 1e-9L; }

    // 事前条件: limit_ms > 0
    long double progress(long limit_ms) const {
        assert(limit_ms > 0);
        return static_cast<long double>(elapsed_ns()) / (static_cast<long double>(limit_ms) * 1'000'000.0L);
    }

    bool is_over(long limit_ms) const { return elapsed_ns() >= limit_ms * 1'000'000L; }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "Timer: elapsed = " << elapsed_ms() << " ms (" << elapsed_us() << " us)\n";
#else
        (void)os;
#endif
    }

   private:
    std::chrono::steady_clock::time_point start_;
};

// 高速擬似乱数生成器 (xoroshiro128++)
// 計算量: 生成 時間 O(1), 空間 O(1) (状態サイズ 16 bytes)
struct Random {
    explicit Random(std::uint64_t seed = 42) { set_seed(seed); }

    void set_seed(std::uint64_t seed) {
        s0_ = splitmix64(seed);
        s1_ = splitmix64(s0_);
        if (s0_ == 0 && s1_ == 0) {
            s0_ = 0x9e3779b97f4a7c15ULL;
        }
    }

    std::uint64_t next_u64() {
        std::uint64_t s0 = s0_;
        std::uint64_t s1 = s1_;
        std::uint64_t res = std::rotl(s0 + s1, 17) + s0;
        s1 ^= s0;
        s0_ = std::rotl(s0, 49) ^ s1 ^ (s1 << 21);
        s1_ = std::rotl(s1, 28);
        return res;
    }

    std::uint32_t next_u32() { return static_cast<std::uint32_t>(next_u64() >> 32); }

    // [0, n) の一様乱数 (事前条件: n > 0)
    int next_int(int n) {
        assert(n > 0);
        return static_cast<int>((static_cast<std::uint64_t>(next_u32()) * static_cast<std::uint32_t>(n)) >> 32);
    }

    // [l, r) の一様乱数 (事前条件: l < r)
    int next_int(int l, int r) {
        assert(l < r);
        std::uint64_t span = static_cast<std::uint64_t>(static_cast<long>(r) - static_cast<long>(l));
        return static_cast<int>(static_cast<long>(l) + static_cast<long>((static_cast<std::uint64_t>(next_u32()) * span) >> 32));
    }

    // [0, n) の一様乱数 (事前条件: n > 0)
    long next_long(long n) {
        assert(n > 0);
        return static_cast<long>((static_cast<__uint128_t>(next_u64()) * static_cast<std::uint64_t>(n)) >> 64);
    }

    // [l, r) の一様乱数 (事前条件: l < r)
    long next_long(long l, long r) {
        assert(l < r);
        std::uint64_t span = static_cast<std::uint64_t>(r) - static_cast<std::uint64_t>(l);
        return static_cast<long>(static_cast<std::uint64_t>(l) +
                                 static_cast<std::uint64_t>((static_cast<__uint128_t>(next_u64()) * span) >> 64));
    }

    // [0, 1) の一様乱数
    long double next_real() { return static_cast<long double>(next_u64()) * 0x1.0p-64L; }

    // [l, r) の一様乱数 (事前条件: l < r)
    long double next_real(long double l, long double r) {
        assert(l < r);
        return l + (r - l) * next_real();
    }

    bool next_bool() { return static_cast<bool>(next_u64() >> 63); }

    bool next_bool(long double p) { return next_real() < p; }

    template <class T>
    void shuffle(std::vector<T>& v) {
        long n = std::ssize(v);
        for (long i = n - 1; i > 0; --i) {
            long j = next_long(i + 1);
            std::swap(v[i], v[j]);
        }
    }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "Random: s0=" << s0_ << ", s1=" << s1_ << "\n";
#else
        (void)os;
#endif
    }

   private:
    std::uint64_t s0_;
    std::uint64_t s1_;
};

// O(1) クリア可能な高速オープンアドレス法 64bit ハッシュ集合
// 計算量: insert / contains 平均時間 O(1), clear 時間 O(1), 空間 O(N)
struct FastHashSet {
    explicit FastHashSet(int init_bits = 10)
        : bits_(init_bits),
          cap_(1L << init_bits),
          mask_(static_cast<std::uint64_t>(cap_ - 1)),
          shift_(64 - init_bits),
          size_(0),
          cur_gen_(1),
          table_(static_cast<std::size_t>(cap_), Entry_{0, 0}) {
        assert(init_bits >= 4 && init_bits <= 28);
    }

    void clear() {
        if (cur_gen_ == std::numeric_limits<int>::max()) {
            for (Entry_& e : table_) {
                e.gen = 0;
            }
            cur_gen_ = 1;
        } else {
            ++cur_gen_;
        }
        size_ = 0;
    }

    bool insert(std::uint64_t key) {
        std::uint64_t idx = (key * 0x9e3779b97f4a7c15ULL) >> shift_;
        while (table_[idx].gen == cur_gen_) {
            if (table_[idx].key == key) {
                return false;
            }
            idx = (idx + 1) & mask_;
        }
        table_[idx] = Entry_{key, cur_gen_};
        ++size_;
        if (size_ * 10 >= cap_ * 7) {
            rehash_(bits_ + 1);
        }
        return true;
    }

    bool contains(std::uint64_t key) const {
        std::uint64_t idx = (key * 0x9e3779b97f4a7c15ULL) >> shift_;
        while (table_[idx].gen == cur_gen_) {
            if (table_[idx].key == key) {
                return true;
            }
            idx = (idx + 1) & mask_;
        }
        return false;
    }

    long size() const { return size_; }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "FastHashSet: size=" << size_ << ", cap=" << cap_ << ", gen=" << cur_gen_ << "\n";
#else
        (void)os;
#endif
    }

   private:
    struct Entry_ {
        std::uint64_t key;
        int gen;
    };

    int bits_;
    long cap_;
    std::uint64_t mask_;
    int shift_;
    long size_;
    int cur_gen_;
    std::vector<Entry_> table_;

    void rehash_(int new_bits) {
        long new_cap = 1L << new_bits;
        std::uint64_t new_mask = static_cast<std::uint64_t>(new_cap - 1);
        int new_shift = 64 - new_bits;
        std::vector<Entry_> new_table(static_cast<std::size_t>(new_cap), Entry_{0, 0});
        for (const Entry_& e : table_) {
            if (e.gen == cur_gen_) {
                std::uint64_t idx = (e.key * 0x9e3779b97f4a7c15ULL) >> new_shift;
                while (new_table[idx].gen == 1) {
                    idx = (idx + 1) & new_mask;
                }
                new_table[idx] = Entry_{e.key, 1};
            }
        }
        table_.swap(new_table);
        bits_ = new_bits;
        cap_ = new_cap;
        mask_ = new_mask;
        shift_ = new_shift;
        cur_gen_ = 1;
    }
};

namespace internal {

inline const long double* get_sa_log_table() {
    static const std::array<long double, 4096> table = [] {
        std::array<long double, 4096> t{};
        t[0] = 0.0L;
        for (long i = 1; i < 4096; ++i) {
            t[i] = std::log((static_cast<long double>(i) + 0.5L) / 4096.0L);
        }
        return t;
    }();
    return table.data();
}

// ビームサーチ重複統合用 O(1) クリア高速ハッシュマップ (key -> int candidate_index)
// 計算量: insert_or_get 平均時間 O(1), clear 時間 O(1), 空間 O(N)
struct FastHashMap {
    explicit FastHashMap(int init_bits = 10)
        : bits_(init_bits),
          cap_(1L << init_bits),
          mask_(static_cast<std::uint64_t>(cap_ - 1)),
          shift_(64 - init_bits),
          size_(0),
          cur_gen_(1),
          table_(static_cast<std::size_t>(cap_), Entry_{0, 0, 0}) {}

    void clear() {
        if (cur_gen_ == std::numeric_limits<int>::max()) {
            for (Entry_& e : table_) {
                e.gen = 0;
            }
            cur_gen_ = 1;
        } else {
            ++cur_gen_;
        }
        size_ = 0;
    }

    int insert_or_get(std::uint64_t key, int new_val) {
        std::uint64_t idx = (key * 0x9e3779b97f4a7c15ULL) >> shift_;
        while (table_[idx].gen == cur_gen_) {
            if (table_[idx].key == key) {
                return table_[idx].val;
            }
            idx = (idx + 1) & mask_;
        }
        table_[idx] = Entry_{key, cur_gen_, new_val};
        ++size_;
        if (size_ * 10 >= cap_ * 7) {
            rehash_(bits_ + 1);
        }
        return -1;
    }

    long size() const { return size_; }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "FastHashMap: size=" << size_ << ", cap=" << cap_ << ", gen=" << cur_gen_ << "\n";
#else
        (void)os;
#endif
    }

   private:
    struct Entry_ {
        std::uint64_t key;
        int gen;
        int val;
    };

    int bits_;
    long cap_;
    std::uint64_t mask_;
    int shift_;
    long size_;
    int cur_gen_;
    std::vector<Entry_> table_;

    void rehash_(int new_bits) {
        long new_cap = 1L << new_bits;
        std::uint64_t new_mask = static_cast<std::uint64_t>(new_cap - 1);
        int new_shift = 64 - new_bits;
        std::vector<Entry_> new_table(static_cast<std::size_t>(new_cap), Entry_{0, 0, 0});
        for (const Entry_& e : table_) {
            if (e.gen == cur_gen_) {
                std::uint64_t idx = (e.key * 0x9e3779b97f4a7c15ULL) >> new_shift;
                while (new_table[idx].gen == 1) {
                    idx = (idx + 1) & new_mask;
                }
                new_table[idx] = Entry_{e.key, 1, e.val};
            }
        }
        table_.swap(new_table);
        bits_ = new_bits;
        cap_ = new_cap;
        mask_ = new_mask;
        shift_ = new_shift;
        cur_gen_ = 1;
    }
};

}  // namespace internal

enum class TempSchedule { Exp, Linear };

// 高速焼きなまし法コントローラ
// 計算量: step() / accept() 償却時間 O(1), 空間 O(1) (log テーブル 64KB は全インスタンス共有)
// 事前条件: Exp では start_temp > 0 && end_temp > 0、Linear では start_temp >= 0 && end_temp >= 0
struct Annealer {
    Random rng;
    Timer timer;

    Annealer(long double start_temp, long double end_temp, long time_limit_ms, TempSchedule schedule = TempSchedule::Exp,
             std::uint64_t seed = 42)
        : Annealer(start_temp, end_temp, time_limit_ms, Timer{}, schedule, seed) {}

    Annealer(long double start_temp, long double end_temp, long time_limit_ms, const Timer& base_timer,
             TempSchedule schedule = TempSchedule::Exp, std::uint64_t seed = 42)
        : rng(seed),
          timer(base_timer),
          start_temp_(start_temp),
          end_temp_(end_temp),
          temp_(start_temp),
          log_start_temp_(0.0L),
          log_temp_ratio_(0.0L),
          progress_(0.0L),
          time_limit_ns_(time_limit_ms * 1'000'000L),
          max_iters_(-1),
          max_check_stride_(16384),
          iters_(0),
          prev_stride_(0),
          check_stride_(1),
          check_rem_(1),
          last_check_ns_(timer.elapsed_ns()),
          schedule_(schedule),
          finished_(false),
          log_table_(internal::get_sa_log_table())
#ifdef LOCAL
          ,
          improved_cnt_(0),
          bad_accept_cnt_(0),
          reject_cnt_(0)
#endif
    {
        assert(time_limit_ms > 0);
        init_schedule_();
    }

    static Annealer by_iters(long double start_temp, long double end_temp, long max_iters, TempSchedule schedule = TempSchedule::Exp,
                             std::uint64_t seed = 42) {
        assert(max_iters > 0);
        Annealer sa(start_temp, end_temp, 1, Timer{}, schedule, seed);
        sa.time_limit_ns_ = -1;
        sa.max_iters_ = max_iters;
        sa.max_check_stride_ = std::max(1L, max_iters / 1000);
        return sa;
    }

    bool step() {
        if (--check_rem_ > 0) {
            return true;
        }
        return update_();
    }

    template <class T>
    bool accept(T delta) {
        return accept_max(delta);
    }

    template <class T, class U>
    bool accept(T new_score, U cur_score) {
        return accept_max(new_score, cur_score);
    }

    template <class T>
    bool accept_max(T delta) {
        if (delta >= 0) {
#ifdef LOCAL
            ++improved_cnt_;
#endif
            return true;
        }
        return accept_negative_(static_cast<long double>(delta));
    }

    template <class T, class U>
    bool accept_max(T new_score, U cur_score) {
        if (new_score >= cur_score) {
#ifdef LOCAL
            ++improved_cnt_;
#endif
            return true;
        }
        return accept_negative_(static_cast<long double>(new_score) - static_cast<long double>(cur_score));
    }

    template <class T>
    bool accept_min(T delta) {
        if (delta <= 0) {
#ifdef LOCAL
            ++improved_cnt_;
#endif
            return true;
        }
        return accept_negative_(-static_cast<long double>(delta));
    }

    template <class T, class U>
    bool accept_min(T new_score, U cur_score) {
        if (new_score <= cur_score) {
#ifdef LOCAL
            ++improved_cnt_;
#endif
            return true;
        }
        return accept_negative_(static_cast<long double>(cur_score) - static_cast<long double>(new_score));
    }

    long double temp() const { return temp_; }

    long double progress() const { return progress_; }

    long iters() const { return iters_ + (prev_stride_ - check_rem_); }

    long elapsed_ms() const { return timer.elapsed_ms(); }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        long total = iters();
        long accepted = improved_cnt_ + bad_accept_cnt_;
        long double accept_rate = total > 0 ? 100.0L * static_cast<long double>(accepted) / static_cast<long double>(total) : 0.0L;
        os << "Annealer: iters=" << total << ", elapsed=" << timer.elapsed_ms() << " ms"
           << ", temp=" << temp_ << ", progress=" << progress_ << ", improved=" << improved_cnt_ << ", bad_accept=" << bad_accept_cnt_
           << ", reject=" << reject_cnt_ << ", accept_rate=" << accept_rate << "%\n";
#else
        (void)os;
#endif
    }

   private:
    long double start_temp_;
    long double end_temp_;
    long double temp_;
    long double log_start_temp_;
    long double log_temp_ratio_;
    long double progress_;
    long time_limit_ns_;
    long max_iters_;
    long max_check_stride_;
    long iters_;
    long prev_stride_;
    long check_stride_;
    long check_rem_;
    long last_check_ns_;
    TempSchedule schedule_;
    bool finished_;
    const long double* log_table_;
#ifdef LOCAL
    long improved_cnt_;
    long bad_accept_cnt_;
    long reject_cnt_;
#endif

    void init_schedule_() {
        if (schedule_ == TempSchedule::Exp) {
            assert(start_temp_ > 0.0L && end_temp_ > 0.0L);
            log_start_temp_ = std::log(start_temp_);
            log_temp_ratio_ = std::log(end_temp_ / start_temp_);
        } else {
            assert(start_temp_ >= 0.0L && end_temp_ >= 0.0L);
        }
    }

    long double next_log_() {
        std::uint64_t r = rng.next_u64();
        long idx = static_cast<long>(r >> 52);
        if (idx != 0) {
            return log_table_[idx];
        }
        return std::log((static_cast<long double>(r) + 0.5L) * 0x1.0p-64L);
    }

    bool accept_negative_(long double neg_delta) {
        bool ok = neg_delta >= temp_ * next_log_();
#ifdef LOCAL
        if (ok) {
            ++bad_accept_cnt_;
        } else {
            ++reject_cnt_;
        }
#endif
        return ok;
    }

    void update_temp_(long double p) {
        if (p < 0.0L) {
            p = 0.0L;
        } else if (p > 1.0L) {
            p = 1.0L;
        }
        progress_ = p;
        if (schedule_ == TempSchedule::Exp) {
            temp_ = std::exp(log_start_temp_ + p * log_temp_ratio_);
        } else {
            temp_ = start_temp_ + (end_temp_ - start_temp_) * p;
        }
    }

    bool update_() {
        if (finished_) {
            return false;
        }
        iters_ += prev_stride_;

        if (max_iters_ > 0) {
            if (iters_ >= max_iters_) {
                prev_stride_ = 0;
                check_rem_ = 0;
                finished_ = true;
                update_temp_(1.0L);
                return false;
            }
            long rem = max_iters_ - iters_;
            check_stride_ = std::min(max_check_stride_, rem);
            prev_stride_ = check_stride_;
            check_rem_ = check_stride_;
            update_temp_(static_cast<long double>(iters_) / static_cast<long double>(max_iters_));
            return true;
        }

        long elapsed_ns = timer.elapsed_ns();
        if (elapsed_ns >= time_limit_ns_) {
            prev_stride_ = 0;
            check_rem_ = 0;
            finished_ = true;
            update_temp_(1.0L);
            return false;
        }

        long dt_ns = elapsed_ns - last_check_ns_;
        last_check_ns_ = elapsed_ns;
        if (prev_stride_ > 0) {
            if (dt_ns < 100'000 && check_stride_ < max_check_stride_) {
                check_stride_ <<= 1;
            } else if (dt_ns > 400'000 && check_stride_ > 1) {
                check_stride_ >>= 1;
            }
            long rem_ns = time_limit_ns_ - elapsed_ns;
            long est_step_ns = std::max(1L, dt_ns / prev_stride_);
            long safe_stride = std::max(1L, rem_ns / (est_step_ns * 2));
            if (check_stride_ > safe_stride) {
                check_stride_ = safe_stride;
            }
        }
        prev_stride_ = check_stride_;
        check_rem_ = check_stride_;
        update_temp_(static_cast<long double>(elapsed_ns) / static_cast<long double>(time_limit_ns_));
        return true;
    }
};

// 2段階展開・コピー型ビームサーチ（軽量状態向け）
// 計算量 (1ターン): 時間 O(M + W log W) (M: 候補数, W: ビーム幅), 空間 O(W * 深さ + M)
// 事前条件: State はコピー/ムーブ可能、Action はコピー可能であること
template <class State, class Action, class Score = long, bool Maximize = true>
struct BeamSearch {
    explicit BeamSearch(const State& init_state, Score init_score = Score{}) : beam_size_(1), turn_(0) { reset(init_state, init_score); }

    void reset(const State& init_state, Score init_score = Score{}) {
        beam_size_ = 1;
        turn_ = 0;
        states_.clear();
        states_.push_back(init_state);
        next_states_.clear();
        scores_.clear();
        scores_.push_back(init_score);
        next_scores_.clear();
        node_ids_.clear();
        node_ids_.push_back(-1);
        next_node_ids_.clear();
        nodes_.clear();
        cands_.clear();
        dedup_.clear();
    }

    void reserve_cands(long n) {
        assert(n >= 0);
        cands_.reserve(static_cast<std::size_t>(n));
    }

    long size() const { return beam_size_; }

    long turn() const { return turn_; }

    const State& state(long i) const {
        assert(0 <= i && i < beam_size_);
        return states_[i];
    }

    Score score(long i) const {
        assert(0 <= i && i < beam_size_);
        return scores_[i];
    }

    void push(Score score, long parent, const Action& action) {
        assert(0 <= parent && parent < beam_size_);
        cands_.push_back(Cand_{score, parent, action});
    }

    void push(Score score, long parent, const Action& action, std::uint64_t hash) {
        assert(0 <= parent && parent < beam_size_);
        int existing = dedup_.insert_or_get(hash, static_cast<int>(std::ssize(cands_)));
        if (existing == -1) {
            cands_.push_back(Cand_{score, parent, action});
        } else if (better_(score, cands_[existing].score)) {
            cands_[existing] = Cand_{score, parent, action};
        }
    }

    template <class ApplyFn>
    bool advance(long next_beam_width, ApplyFn apply, bool sort_beam = true) {
        assert(next_beam_width > 0);
        long total = std::ssize(cands_);
        if (total == 0) {
            dedup_.clear();
            return false;
        }

        long keep = std::min(next_beam_width, total);
        if (keep < total) {
            std::nth_element(cands_.begin(), cands_.begin() + keep, cands_.end(),
                             [](const Cand_& a, const Cand_& b) { return better_(a.score, b.score); });
        }
        if (sort_beam && keep > 1) {
            std::sort(cands_.begin(), cands_.begin() + keep, [](const Cand_& a, const Cand_& b) { return better_(a.score, b.score); });
        }

        parent_rem_.assign(static_cast<std::size_t>(beam_size_), 0);
        for (long i = 0; i < keep; ++i) {
            ++parent_rem_[cands_[i].parent];
        }

        next_scores_.clear();
        next_node_ids_.clear();
        for (long i = 0; i < keep; ++i) {
            long p = cands_[i].parent;
            --parent_rem_[p];
            if (i < std::ssize(next_states_)) {
                if (parent_rem_[p] == 0) {
                    std::swap(next_states_[i], states_[p]);
                } else {
                    next_states_[i] = states_[p];
                }
                apply(next_states_[i], cands_[i].action);
            } else {
                if (parent_rem_[p] == 0) {
                    next_states_.push_back(std::move(states_[p]));
                } else {
                    next_states_.push_back(states_[p]);
                }
                apply(next_states_.back(), cands_[i].action);
            }
            next_scores_.push_back(cands_[i].score);
            long nid = std::ssize(nodes_);
            nodes_.push_back(Node_{node_ids_[p], cands_[i].action});
            next_node_ids_.push_back(nid);
        }

        beam_size_ = keep;
        ++turn_;
        states_.swap(next_states_);
        scores_.swap(next_scores_);
        node_ids_.swap(next_node_ids_);
        cands_.clear();
        dedup_.clear();
        return true;
    }

    template <class ExpandFn, class ApplyFn>
    long run(long max_turns, long beam_width, ExpandFn expand, ApplyFn apply, bool sort_beam = true) {
        assert(max_turns >= 0 && beam_width > 0);
        for (long t = 0; t < max_turns; ++t) {
            for (long i = 0; i < beam_size_; ++i) {
                expand(i, states_[i]);
            }
            if (!advance(beam_width, apply, sort_beam)) {
                return t;
            }
        }
        return max_turns;
    }

    long best_index() const {
        assert(beam_size_ > 0);
        long best_i = 0;
        for (long i = 1; i < beam_size_; ++i) {
            if (better_(scores_[i], scores_[best_i])) {
                best_i = i;
            }
        }
        return best_i;
    }

    const State& best_state() const { return states_[best_index()]; }

    Score best_score() const { return scores_[best_index()]; }

    std::vector<Action> restore_actions(long idx = -1) const {
        if (idx == -1) {
            idx = best_index();
        }
        assert(0 <= idx && idx < beam_size_);
        std::vector<Action> actions;
        long nid = node_ids_[idx];
        while (nid != -1) {
            actions.push_back(nodes_[nid].action);
            nid = nodes_[nid].parent;
        }
        std::reverse(actions.begin(), actions.end());
        return actions;
    }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "BeamSearch: turn=" << turn_ << ", beam_size=" << beam_size_ << ", best_score=" << best_score()
           << ", total_nodes=" << std::ssize(nodes_) << "\n";
#else
        (void)os;
#endif
    }

   private:
    struct Cand_ {
        Score score;
        long parent;
        Action action;
    };

    struct Node_ {
        long parent;
        Action action;
    };

    long beam_size_;
    long turn_;
    std::vector<State> states_;
    std::vector<State> next_states_;
    std::vector<Score> scores_;
    std::vector<Score> next_scores_;
    std::vector<long> node_ids_;
    std::vector<long> next_node_ids_;
    std::vector<Node_> nodes_;
    std::vector<Cand_> cands_;
    std::vector<int> parent_rem_;
    internal::FastHashMap dedup_;

    static bool better_(Score a, Score b) {
        if constexpr (Maximize) {
            return a > b;
        } else {
            return a < b;
        }
    }
};

// 単一状態ロールバック型 (Euler Tour) ビームサーチ（盤面などの巨大状態向け）
// 計算量 (1ターン): 時間 O(木エッジ数 + M + W), 空間 O(W + 深さ + M)
// 事前条件: Action はデフォルト構築・コピー可能であり、apply(act) と undo(act) が逆操作であること
template <class Action, class Score = long, bool Maximize = true>
struct EulerBeamSearch {
    EulerBeamSearch()
        : root_(0),
          cur_node_(0),
          last_best_node_(0),
          has_last_best_(false),
          has_overall_best_(false),
          turn_(0),
          last_best_score_{},
          overall_best_score_{} {
        reset();
    }

    void reset() {
        nodes_.clear();
        free_ids_.clear();
        leaves_.clear();
        next_leaves_.clear();
        cands_.clear();
        committed_actions_.clear();
        overall_best_actions_.clear();
        dedup_.clear();
        nodes_.push_back(Node_{Action{}, -1, -1, -1, -1});
        root_ = 0;
        cur_node_ = 0;
        last_best_node_ = 0;
        has_last_best_ = false;
        has_overall_best_ = false;
        turn_ = 0;
        leaves_.push_back(0);
    }

    void reserve_cands(long n) {
        assert(n >= 0);
        cands_.reserve(static_cast<std::size_t>(n));
    }

    long turn() const { return turn_; }

    long size() const { return std::ssize(leaves_); }

    void push(Score score, const Action& action) { cands_.push_back(Cand_{score, cur_node_, action}); }

    void push(Score score, const Action& action, std::uint64_t hash) {
        int existing = dedup_.insert_or_get(hash, static_cast<int>(std::ssize(cands_)));
        if (existing == -1) {
            cands_.push_back(Cand_{score, cur_node_, action});
        } else if (better_(score, cands_[existing].score)) {
            cands_[existing] = Cand_{score, cur_node_, action};
        }
    }

    void record_best(Score score) {
        if (!has_overall_best_ || better_(score, overall_best_score_)) {
            has_overall_best_ = true;
            overall_best_score_ = score;
            overall_best_actions_ = path_from_root_(cur_node_);
        }
    }

    template <class ApplyFn, class UndoFn, class ExpandFn>
    bool step(long beam_width, ApplyFn apply, UndoFn undo, ExpandFn expand, bool track_overall_best = false) {
        assert(beam_width > 0);
        cands_.clear();
        dedup_.clear();

        if (nodes_[root_].child == -1) {
            cur_node_ = root_;
            expand();
        } else {
            int v = nodes_[root_].child;
            apply(nodes_[v].action);
            while (v != root_) {
                if (nodes_[v].child != -1) {
                    v = nodes_[v].child;
                    apply(nodes_[v].action);
                } else {
                    cur_node_ = v;
                    expand();
                    while (v != root_) {
                        undo(nodes_[v].action);
                        if (nodes_[v].next_sib != -1) {
                            v = nodes_[v].next_sib;
                            apply(nodes_[v].action);
                            break;
                        }
                        v = nodes_[v].parent;
                    }
                }
            }
        }

        long total = std::ssize(cands_);
        if (total == 0) {
            return false;
        }

        long keep = std::min(beam_width, total);
        if (keep < total) {
            std::nth_element(cands_.begin(), cands_.begin() + keep, cands_.end(),
                             [](const Cand_& a, const Cand_& b) { return better_(a.score, b.score); });
        }

        next_leaves_.clear();
        long best_i = 0;
        for (long i = 0; i < keep; ++i) {
            if (i > 0 && better_(cands_[i].score, cands_[best_i].score)) {
                best_i = i;
            }
            int nid = add_child_(cands_[i].parent, cands_[i].action);
            next_leaves_.push_back(nid);
            if (i == best_i) {
                last_best_node_ = nid;
            }
        }
        has_last_best_ = true;
        last_best_score_ = cands_[best_i].score;

        for (int u : leaves_) {
            if (nodes_[u].child == -1) {
                erase_branch_(u);
            }
        }
        leaves_.swap(next_leaves_);

        while (nodes_[root_].child != -1 && nodes_[nodes_[root_].child].next_sib == -1) {
            int c = nodes_[root_].child;
            apply(nodes_[c].action);
            committed_actions_.push_back(nodes_[c].action);
            free_ids_.push_back(root_);
            root_ = c;
            nodes_[root_].parent = -1;
        }

        if (track_overall_best) {
            if (!has_overall_best_ || better_(last_best_score_, overall_best_score_)) {
                has_overall_best_ = true;
                overall_best_score_ = last_best_score_;
                overall_best_actions_ = path_from_root_(last_best_node_);
            }
        }

        ++turn_;
        return true;
    }

    template <class ApplyFn, class UndoFn, class ExpandFn>
    long run(long max_turns, long beam_width, ApplyFn apply, UndoFn undo, ExpandFn expand, bool track_overall_best = false) {
        assert(max_turns >= 0 && beam_width > 0);
        for (long t = 0; t < max_turns; ++t) {
            if (!step(beam_width, apply, undo, expand, track_overall_best)) {
                return t;
            }
        }
        return max_turns;
    }

    bool has_best() const { return has_last_best_; }

    Score best_score() const {
        assert(has_last_best_);
        return last_best_score_;
    }

    std::vector<Action> best_actions() const {
        if (!has_last_best_) {
            return {};
        }
        return path_from_root_(last_best_node_);
    }

    bool has_overall_best() const { return has_overall_best_; }

    Score overall_best_score() const {
        assert(has_overall_best_);
        return overall_best_score_;
    }

    std::vector<Action> overall_best_actions() const {
        if (!has_overall_best_) {
            return {};
        }
        return overall_best_actions_;
    }

    void dump(std::ostream& os = std::cerr) const {
#ifdef LOCAL
        os << "EulerBeamSearch: turn=" << turn_ << ", leaves=" << std::ssize(leaves_)
           << ", active_nodes=" << (std::ssize(nodes_) - std::ssize(free_ids_)) << ", committed=" << std::ssize(committed_actions_);
        if (has_last_best_) {
            os << ", last_best_score=" << last_best_score_;
        }
        if (has_overall_best_) {
            os << ", overall_best_score=" << overall_best_score_;
        }
        os << "\n";
#else
        (void)os;
#endif
    }

   private:
    struct Node_ {
        Action action;
        int parent;
        int child;
        int next_sib;
        int prev_sib;
    };

    struct Cand_ {
        Score score;
        int parent;
        Action action;
    };

    int root_;
    int cur_node_;
    int last_best_node_;
    bool has_last_best_;
    bool has_overall_best_;
    long turn_;
    Score last_best_score_;
    Score overall_best_score_;
    std::vector<Node_> nodes_;
    std::vector<int> free_ids_;
    std::vector<int> leaves_;
    std::vector<int> next_leaves_;
    std::vector<Cand_> cands_;
    std::vector<Action> committed_actions_;
    std::vector<Action> overall_best_actions_;
    internal::FastHashMap dedup_;

    static bool better_(Score a, Score b) {
        if constexpr (Maximize) {
            return a > b;
        } else {
            return a < b;
        }
    }

    int add_child_(int p, const Action& act) {
        int id;
        if (!free_ids_.empty()) {
            id = free_ids_.back();
            free_ids_.pop_back();
            nodes_[id] = Node_{act, p, -1, nodes_[p].child, -1};
        } else {
            id = static_cast<int>(std::ssize(nodes_));
            nodes_.push_back(Node_{act, p, -1, nodes_[p].child, -1});
        }
        int old_first = nodes_[p].child;
        if (old_first != -1) {
            nodes_[old_first].prev_sib = id;
        }
        nodes_[p].child = id;
        return id;
    }

    void erase_branch_(int u) {
        while (u != root_ && nodes_[u].child == -1) {
            int p = nodes_[u].parent;
            int prev = nodes_[u].prev_sib;
            int next = nodes_[u].next_sib;
            if (prev != -1) {
                nodes_[prev].next_sib = next;
            } else {
                nodes_[p].child = next;
            }
            if (next != -1) {
                nodes_[next].prev_sib = prev;
            }
            free_ids_.push_back(u);
            u = p;
        }
    }

    std::vector<Action> path_from_root_(int u) const {
        std::vector<Action> tail;
        while (u != root_) {
            tail.push_back(nodes_[u].action);
            u = nodes_[u].parent;
        }
        std::vector<Action> res;
        res.reserve(committed_actions_.size() + tail.size());
        res.insert(res.end(), committed_actions_.begin(), committed_actions_.end());
        res.insert(res.end(), tail.rbegin(), tail.rend());
        return res;
    }
};