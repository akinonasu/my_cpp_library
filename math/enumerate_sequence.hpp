#include <algorithm>
#include <cassert>
#include <concepts>
#include <limits>
#include <print>
#include <type_traits>
#include <utility>
#include <vector>

static_assert(sizeof(long) == 8, "64-bit long is required.");

namespace seq {

namespace detail {

using i128 = __int128;

constexpr i128 abs_i128(i128 x) noexcept {
    return x < 0 ? -x : x;
}

constexpr i128 floor_div_i128(i128 a, i128 b) noexcept {
    assert(b != 0);
    i128 q = a / b;
    i128 r = a % b;
    if (r != 0 && ((r < 0) != (b < 0))) --q;
    return q;
}

constexpr i128 ceil_div_i128(i128 a, i128 b) noexcept {
    assert(b != 0);
    i128 q = a / b;
    i128 r = a % b;
    if (r != 0 && ((r < 0) == (b < 0))) ++q;
    return q;
}

struct IndexBounds {
    bool ok;
    i128 k_min;
    i128 k_max;
};

// d != 0 のとき、k >= k_start かつ lo <= a + k * d <= hi を満たす k の範囲を求める
constexpr IndexBounds ap_index_bounds(long a, long d, long lo, long hi, long k_start) noexcept {
    assert(d != 0);
    assert(lo <= hi);
    assert(k_start >= 0);

    i128 lo_k, hi_k;
    if (d > 0) {
        lo_k = ceil_div_i128(static_cast<i128>(lo) - a, d);
        hi_k = floor_div_i128(static_cast<i128>(hi) - a, d);
    } else {
        lo_k = ceil_div_i128(static_cast<i128>(hi) - a, d);
        hi_k = floor_div_i128(static_cast<i128>(lo) - a, d);
    }
    if (lo_k < k_start) lo_k = k_start;
    if (lo_k > hi_k) return {false, 0, -1};
    return {true, lo_k, hi_k};
}

}  // namespace detail

template <class Func>
concept Callback = std::invocable<Func, long> && (std::same_as<std::invoke_result_t<Func, long>, void> || std::same_as<std::invoke_result_t<Func, long>, bool>);

template <class Pred>
concept LongPredicate = std::predicate<Pred, long>;

struct ApConfig {
    long a = 0;
    long d = 1;
    long lo = std::numeric_limits<long>::min();
    long hi = std::numeric_limits<long>::max();
    long k_start = 0;
    long limit = -1;        // -1: 無制限
    bool distinct = false;  // true なら d == 0 でも最大1項のみ
};

struct GpConfig {
    long a = 1;
    long r = 2;
    long lo = std::numeric_limits<long>::min();
    long hi = std::numeric_limits<long>::max();
    long k_start = 0;
    long limit = -1;        // -1: 無制限
    bool distinct = false;  // true なら重複値（a=0, r=0, r=1, r=-1）を1度のみ出力
};

// ==========================================
// 等差数列 (Arithmetic Progression)
// ==========================================

template <Callback Func>
void for_each_ap(const ApConfig& cfg, Func&& func) {
    assert(cfg.lo <= cfg.hi);
    assert(cfg.k_start >= 0);
    assert(cfg.limit >= -1);
    if (cfg.limit == 0) return;

    constexpr bool returns_bool = std::same_as<std::invoke_result_t<Func, long>, bool>;
    auto emit = [&](long val) -> bool {
        if constexpr (returns_bool) {
            return func(val);
        } else {
            func(val);
            return true;
        }
    };

    if (cfg.d == 0) {
        if (cfg.a < cfg.lo || cfg.a > cfg.hi) return;
        if (cfg.distinct) {
            emit(cfg.a);
            return;
        }
        if constexpr (!returns_bool) {
            assert(cfg.limit != -1 && "Infinite arithmetic progression with d == 0");
        }
        for (long cnt = 0; cfg.limit == -1 || cnt < cfg.limit; ++cnt) {
            if (!emit(cfg.a)) return;
        }
        return;
    }

    auto [ok, k_min, k_max] = detail::ap_index_bounds(cfg.a, cfg.d, cfg.lo, cfg.hi, cfg.k_start);
    if (!ok) return;

    detail::i128 cur = static_cast<detail::i128>(cfg.a) + k_min * static_cast<detail::i128>(cfg.d);
    long cnt = 0;
    for (detail::i128 k = k_min; k <= k_max; ++k) {
        if (cfg.limit != -1 && cnt >= cfg.limit) break;
        if (!emit(static_cast<long>(cur))) break;
        ++cnt;
        cur += cfg.d;
    }
}

template <LongPredicate WhilePred, LongPredicate FilterPred>
std::vector<long> ap_enum_while(const ApConfig& cfg, WhilePred&& cond, FilterPred&& filter) {
    assert(cfg.lo <= cfg.hi);
    assert(cfg.k_start >= 0);
    assert(cfg.limit >= -1);

    std::vector<long> res;
    if (cfg.limit == 0) return res;

    if (cfg.d == 0) {
        if (cfg.a < cfg.lo || cfg.a > cfg.hi) return res;
        if (!cond(cfg.a) || !filter(cfg.a)) return res;
        if (cfg.distinct) {
            res.push_back(cfg.a);
            return res;
        }
        assert(cfg.limit != -1 && "Infinite arithmetic progression with d == 0");
        res.assign(static_cast<std::size_t>(cfg.limit), cfg.a);
        return res;
    }

    ApConfig scan = cfg;
    scan.limit = -1;  // 採用数ベースで limit 判定するため走査自体は抑えない
    for_each_ap(scan, [&](long x) -> bool {
        if (!cond(x)) return false;
        if (filter(x)) {
            res.push_back(x);
            if (cfg.limit != -1 && static_cast<long>(res.size()) >= cfg.limit) return false;
        }
        return true;
    });
    return res;
}

template <LongPredicate WhilePred>
std::vector<long> ap_enum_while(const ApConfig& cfg, WhilePred&& cond) {
    return ap_enum_while(cfg, std::forward<WhilePred>(cond), [](long) { return true; });
}

template <LongPredicate FilterPred>
std::vector<long> ap_enum(const ApConfig& cfg, FilterPred&& filter) {
    return ap_enum_while(cfg, [](long) { return true; }, std::forward<FilterPred>(filter));
}

inline std::vector<long> ap_enum(const ApConfig& cfg) {
    return ap_enum_while(cfg, [](long) { return true; }, [](long) { return true; });
}

// [lo, hi] に含まれる等差数列を列挙（d == 0 で範囲内の場合は assert）
inline std::vector<long> ap_range(long a, long d, long lo, long hi) {
    assert(lo <= hi);
    assert(d != 0 || a < lo || a > hi);
    return ap_enum(ApConfig{.a = a, .d = d, .lo = lo, .hi = hi});
}

// 初項 a, 公差 d から最大 n 項を列挙
inline std::vector<long> ap_take(long a, long d, long n, bool stop_on_overflow = true) {
    assert(n >= 0);
    std::vector<long> res;
    if (n == 0) return res;
    res.reserve(static_cast<std::size_t>(std::min<long>(n, 1'000'000)));
    long cur = a;
    res.push_back(cur);
    for (long i = 1; i < n; ++i) {
        long nxt = 0;
        if (__builtin_add_overflow(cur, d, &nxt)) {
            assert(stop_on_overflow && "Arithmetic progression overflowed long");
            break;
        }
        cur = nxt;
        res.push_back(cur);
    }
    return res;
}

// x ≡ rem (mod step) かつ lo <= x <= hi を満たす整数を昇順列挙 (step > 0)
inline std::vector<long> ap_mod_range(long rem, long step, long lo, long hi) {
    assert(lo <= hi);
    assert(step > 0);
    detail::i128 s = step;
    detail::i128 r_mod = ((static_cast<detail::i128>(rem) % s) + s) % s;
    detail::i128 l_mod = ((static_cast<detail::i128>(lo) % s) + s) % s;
    detail::i128 offset = (r_mod - l_mod + s) % s;
    detail::i128 first = static_cast<detail::i128>(lo) + offset;
    if (first > hi) return {};
    return ap_range(static_cast<long>(first), step, lo, hi);
}

// --- AP の O(1) クエリ ---

// 個数を返す（d == 0 かつ範囲内で無限の場合は -1）
inline long ap_count_in_range(long a, long d, long lo, long hi) {
    assert(lo <= hi);
    if (d == 0) return (lo <= a && a <= hi) ? -1 : 0;
    auto [ok, k_min, k_max] = detail::ap_index_bounds(a, d, lo, hi, 0);
    if (!ok) return 0;
    detail::i128 cnt = k_max - k_min + 1;
    assert(cnt <= std::numeric_limits<long>::max() && "Count overflows long");
    return static_cast<long>(cnt);
}

// 範囲に入る最小の k >= 0（存在しない場合は -1）
inline long ap_first_index(long a, long d, long lo, long hi) {
    assert(lo <= hi);
    if (d == 0) return (lo <= a && a <= hi) ? 0 : -1;
    auto [ok, k_min, k_max] = detail::ap_index_bounds(a, d, lo, hi, 0);
    if (!ok) return -1;
    assert(k_min <= std::numeric_limits<long>::max());
    return static_cast<long>(k_min);
}

// 範囲に入る最大の k >= 0（存在しない、または無限の場合は -1）
inline long ap_last_index(long a, long d, long lo, long hi) {
    assert(lo <= hi);
    if (d == 0) return -1;
    auto [ok, k_min, k_max] = detail::ap_index_bounds(a, d, lo, hi, 0);
    if (!ok) return -1;
    assert(k_max <= std::numeric_limits<long>::max());
    return static_cast<long>(k_max);
}

// 範囲に入る最初の値 {ok, val}
inline std::pair<bool, long> ap_first_value(long a, long d, long lo, long hi) {
    assert(lo <= hi);
    if (d == 0) {
        if (lo <= a && a <= hi) return {true, a};
        return {false, 0};
    }
    auto [ok, k_min, k_max] = detail::ap_index_bounds(a, d, lo, hi, 0);
    if (!ok) return {false, 0};
    detail::i128 val = static_cast<detail::i128>(a) + k_min * static_cast<detail::i128>(d);
    return {true, static_cast<long>(val)};
}

// 範囲に入る最後の値 {ok, val}（無限の場合は {false, 0}）
inline std::pair<bool, long> ap_last_value(long a, long d, long lo, long hi) {
    assert(lo <= hi);
    if (d == 0) return {false, 0};
    auto [ok, k_min, k_max] = detail::ap_index_bounds(a, d, lo, hi, 0);
    if (!ok) return {false, 0};
    detail::i128 val = static_cast<detail::i128>(a) + k_max * static_cast<detail::i128>(d);
    return {true, static_cast<long>(val)};
}

// ==========================================
// 等比数列 (Geometric Progression)
// ==========================================

template <Callback Func>
void for_each_gp(const GpConfig& cfg, Func&& func) {
    assert(cfg.lo <= cfg.hi);
    assert(cfg.k_start >= 0);
    assert(cfg.limit >= -1);
    if (cfg.limit == 0) return;

    constexpr bool returns_bool = std::same_as<std::invoke_result_t<Func, long>, bool>;
    long emitted = 0;
    auto emit = [&](long val) -> bool {
        if (cfg.limit != -1 && emitted >= cfg.limit) return false;
        ++emitted;
        if constexpr (returns_bool) {
            return func(val);
        } else {
            func(val);
            return true;
        }
    };

    auto emit_const = [&](long val) {
        if (val < cfg.lo || val > cfg.hi) return;
        if (cfg.distinct) {
            emit(val);
            return;
        }
        if constexpr (!returns_bool) {
            assert(cfg.limit != -1 && "Infinite geometric progression");
        }
        while (cfg.limit == -1 || emitted < cfg.limit) {
            if (!emit(val)) return;
        }
    };

    if (cfg.a == 0) {
        emit_const(0);
        return;
    }

    if (cfg.r == 0) {
        if (cfg.k_start == 0 && cfg.lo <= cfg.a && cfg.a <= cfg.hi) {
            if (!emit(cfg.a)) return;
        }
        emit_const(0);
        return;
    }

    if (cfg.r == 1) {
        emit_const(cfg.a);
        return;
    }

    if (cfg.r == -1) {
        detail::i128 v0 = (cfg.k_start % 2 == 0) ? static_cast<detail::i128>(cfg.a) : -static_cast<detail::i128>(cfg.a);
        detail::i128 v1 = -v0;
        bool in0 = (cfg.lo <= v0 && v0 <= cfg.hi);
        bool in1 = (cfg.lo <= v1 && v1 <= cfg.hi);
        if (!in0 && !in1) return;

        if (cfg.distinct) {
            if (in0 && !emit(static_cast<long>(v0))) return;
            if (in1) emit(static_cast<long>(v1));
            return;
        }

        if constexpr (!returns_bool) {
            assert(cfg.limit != -1 && "Infinite geometric progression with r == -1");
        }
        while (cfg.limit == -1 || emitted < cfg.limit) {
            if (in0 && !emit(static_cast<long>(v0))) return;
            if (in1 && !emit(static_cast<long>(v1))) return;
        }
        return;
    }

    // |r| >= 2 かつ a != 0
    detail::i128 max_abs = std::max(detail::abs_i128(cfg.lo), detail::abs_i128(cfg.hi));
    detail::i128 cur = cfg.a;
    for (long k = 0; k < cfg.k_start; ++k) {
        if (detail::abs_i128(cur) > max_abs) return;
        cur *= cfg.r;
    }
    while (detail::abs_i128(cur) <= max_abs) {
        if (cfg.lo <= cur && cur <= cfg.hi) {
            if (!emit(static_cast<long>(cur))) return;
        }
        cur *= cfg.r;
    }
}

template <LongPredicate WhilePred, LongPredicate FilterPred>
std::vector<long> gp_enum_while(const GpConfig& cfg, WhilePred&& cond, FilterPred&& filter) {
    assert(cfg.lo <= cfg.hi);
    assert(cfg.k_start >= 0);
    assert(cfg.limit >= -1);

    std::vector<long> res;
    if (cfg.limit == 0) return res;

    auto append_const = [&](long val) {
        if (val < cfg.lo || val > cfg.hi) return;
        if (!cond(val) || !filter(val)) return;
        if (cfg.distinct) {
            res.push_back(val);
            return;
        }
        assert(cfg.limit != -1 && "Infinite geometric progression");
        while (static_cast<long>(res.size()) < cfg.limit) {
            res.push_back(val);
        }
    };

    if (cfg.a == 0) {
        append_const(0);
        return res;
    }

    if (cfg.r == 0) {
        if (cfg.k_start == 0 && cfg.lo <= cfg.a && cfg.a <= cfg.hi) {
            if (!cond(cfg.a)) return res;
            if (filter(cfg.a)) {
                res.push_back(cfg.a);
                if (cfg.limit != -1 && static_cast<long>(res.size()) >= cfg.limit) return res;
            }
        }
        append_const(0);
        return res;
    }

    if (cfg.r == 1) {
        append_const(cfg.a);
        return res;
    }

    if (cfg.r == -1) {
        detail::i128 v0 = (cfg.k_start % 2 == 0) ? static_cast<detail::i128>(cfg.a) : -static_cast<detail::i128>(cfg.a);
        detail::i128 v1 = -v0;
        bool in0 = (cfg.lo <= v0 && v0 <= cfg.hi);
        bool in1 = (cfg.lo <= v1 && v1 <= cfg.hi);
        if (!in0 && !in1) return res;

        bool added_in_cycle = false;
        auto step_one = [&](long val) -> bool {
            if (!cond(val)) return false;
            if (filter(val)) {
                res.push_back(val);
                added_in_cycle = true;
                if (cfg.limit != -1 && static_cast<long>(res.size()) >= cfg.limit) return false;
            }
            return true;
        };

        if (in0 && !step_one(static_cast<long>(v0))) return res;
        if (in1 && !step_one(static_cast<long>(v1))) return res;
        if (cfg.distinct || !added_in_cycle) return res;

        assert(cfg.limit != -1 && "Infinite geometric progression with r == -1");
        while (true) {
            if (in0 && !step_one(static_cast<long>(v0))) return res;
            if (in1 && !step_one(static_cast<long>(v1))) return res;
        }
    }

    // |r| >= 2
    GpConfig scan = cfg;
    scan.limit = -1;
    for_each_gp(scan, [&](long x) -> bool {
        if (!cond(x)) return false;
        if (filter(x)) {
            res.push_back(x);
            if (cfg.limit != -1 && static_cast<long>(res.size()) >= cfg.limit) return false;
        }
        return true;
    });
    return res;
}

template <LongPredicate WhilePred>
std::vector<long> gp_enum_while(const GpConfig& cfg, WhilePred&& cond) {
    return gp_enum_while(cfg, std::forward<WhilePred>(cond), [](long) { return true; });
}

template <LongPredicate FilterPred>
std::vector<long> gp_enum(const GpConfig& cfg, FilterPred&& filter) {
    return gp_enum_while(cfg, [](long) { return true; }, std::forward<FilterPred>(filter));
}

inline std::vector<long> gp_enum(const GpConfig& cfg) {
    return gp_enum_while(cfg, [](long) { return true; }, [](long) { return true; });
}

// [lo, hi] に含まれる等比数列を列挙（デフォルトで重複なしのため無限にならない）
inline std::vector<long> gp_range(long a, long r, long lo, long hi, bool distinct = true) {
    assert(lo <= hi);
    return gp_enum(GpConfig{.a = a, .r = r, .lo = lo, .hi = hi, .distinct = distinct});
}

// 初項 a, 公比 r から最大 n 項を列挙
inline std::vector<long> gp_take(long a, long r, long n, bool stop_on_overflow = true) {
    assert(n >= 0);
    std::vector<long> res;
    if (n == 0) return res;
    res.reserve(static_cast<std::size_t>(std::min<long>(n, 1'000'000)));
    long cur = a;
    res.push_back(cur);
    for (long i = 1; i < n; ++i) {
        long nxt = 0;
        if (__builtin_mul_overflow(cur, r, &nxt)) {
            assert(stop_on_overflow && "Geometric progression overflowed long");
            break;
        }
        cur = nxt;
        res.push_back(cur);
    }
    return res;
}

// --- GP のクエリ ---

// 個数を返す（distinct == false で無限の場合は -1）
inline long gp_count_in_range(long a, long r, long lo, long hi, bool distinct = true) {
    assert(lo <= hi);
    if (!distinct) {
        if (a == 0) return (lo <= 0 && 0 <= hi) ? -1 : 0;
        if (r == 0) {
            if (lo <= 0 && 0 <= hi) return -1;
            return (lo <= a && a <= hi) ? 1 : 0;
        }
        if (r == 1) return (lo <= a && a <= hi) ? -1 : 0;
        if (r == -1) {
            detail::i128 na = -static_cast<detail::i128>(a);
            if ((lo <= a && a <= hi) || (lo <= na && na <= hi)) return -1;
            return 0;
        }
    }
    return static_cast<long>(gp_range(a, r, lo, hi, distinct).size());
}

// 範囲に入る最初の値 {ok, val}
inline std::pair<bool, long> gp_first_value(long a, long r, long lo, long hi) {
    assert(lo <= hi);
    std::pair<bool, long> ans = {false, 0};
    for_each_gp(GpConfig{.a = a, .r = r, .lo = lo, .hi = hi, .limit = 1}, [&](long x) -> bool {
        ans = {true, x};
        return false;
    });
    return ans;
}

}  // namespace seq
