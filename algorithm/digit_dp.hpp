#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>
#include <string_view>
#include <utility>
#include <vector>

namespace internal {

template <class T>
void reset_table(T& x) {
    if constexpr (requires {
                      typename T::key_type;
                      x.clear();
                  }) {
        x.clear();
    } else if constexpr (requires {
                             x.begin();
                             x.end();
                         }) {
        for (auto& e : x)
            reset_table(e);
    } else if constexpr (requires {
                             x.first;
                             x.second;
                         }) {
        reset_table(x.first);
        reset_table(x.second);
    } else {
        x = T{};
    }
}

inline std::vector<int> to_digits(long n, int base) {
    assert(n >= 0);
    assert(base >= 2);
    if (n == 0) return {0};
    std::vector<int> res;
    while (n > 0) {
        res.push_back(static_cast<int>(n % base));
        n /= base;
    }
    std::reverse(res.begin(), res.end());
    return res;
}

inline std::vector<int> to_digits(std::string_view s, int base) {
    assert(!s.empty());
    assert(2 <= base && base <= 10);
    std::vector<int> res;
    res.reserve(s.size());
    for (char c : s) {
        int d = c - '0';
        assert(0 <= d && d < base);
        res.push_back(d);
    }
    return res;
}

inline void pad_digits(std::vector<int>& dl, std::vector<int>& dr) {
    if (dl.size() < dr.size()) {
        std::vector<int> tmp(dr.size() - dl.size(), 0);
        tmp.insert(tmp.end(), dl.begin(), dl.end());
        dl = std::move(tmp);
    } else if (dr.size() < dl.size()) {
        std::vector<int> tmp(dl.size() - dr.size(), 0);
        tmp.insert(tmp.end(), dr.begin(), dr.end());
        dr = std::move(tmp);
    }
}

template <class F, class Table>
void call_trans(F& f, const Table& cur, Table& nxt, int d, long pos) {
    if constexpr (requires { f(cur, nxt, d, pos); }) {
        f(cur, nxt, d, pos);
    } else {
        f(cur, nxt, d);
    }
}

template <class F, class Table>
void call_start(F& f, Table& nxt, int d, long pos) {
    if constexpr (requires { f(nxt, d, pos); }) {
        f(nxt, d, pos);
    } else {
        f(nxt, d);
    }
}

}  // namespace internal

template <class Table>
struct DigitDPResult {
    Table less;   // 0 <= x < N
    Table equal;  // x == N
};

template <class Table>
struct DigitDPNZResult {
    Table less;       // 1 <= x < N
    Table equal;      // x == N (N >= 1)
    bool zero_less;   // 0 < N か (N >= 1 なら true)
    bool zero_equal;  // N == 0 か
};

template <class Table>
struct DigitDPRangeResult {
    Table dp[2][2];  // dp[tight_l][tight_r]
};

template <class Table>
struct DigitDPRangeNZResult {
    Table dp[2][2];   // dp[tight_l][tight_r] (max(1, L) <= x <= R)
    bool zero[2][2];  // 0 が [L, R] に含まれるか
};

namespace internal {

template <class Table, class Start, class Trans>
DigitDPNZResult<Table> digit_dp_nz_impl(const std::vector<int>& digits, const Table& empty, Start start, Trans trans, int base) {
    assert(!digits.empty());
    assert(base >= 2);
    for (int d : digits)
        assert(0 <= d && d < base);

    Table cur[2] = {empty, empty};
    Table nxt[2] = {empty, empty};
    bool lz[2] = {false, true};
    bool nlz[2] = {false, false};

    long n = std::ssize(digits);
    for (long pos = 0; pos < n; ++pos) {
        int lim_tight = digits[pos];
        for (int tight = 0; tight < 2; ++tight) {
            int lim = tight ? lim_tight : base - 1;
            for (int d = 0; d <= lim; ++d) {
                int nt = tight && (d == lim);
                call_trans(trans, cur[tight], nxt[nt], d, pos);
            }
        }
        for (int tight = 0; tight < 2; ++tight) {
            if (!lz[tight]) continue;
            int lim = tight ? lim_tight : base - 1;
            nlz[tight && (0 == lim)] = true;
            for (int d = 1; d <= lim; ++d) {
                int nt = tight && (d == lim);
                call_start(start, nxt[nt], d, pos);
            }
        }
        std::swap(cur[0], nxt[0]);
        std::swap(cur[1], nxt[1]);
        nxt[0] = empty;
        nxt[1] = empty;
        lz[0] = nlz[0];
        lz[1] = nlz[1];
        nlz[0] = nlz[1] = false;
    }
    return {std::move(cur[0]), std::move(cur[1]), lz[0], lz[1]};
}

template <class Table, class Start, class Trans>
DigitDPRangeNZResult<Table> digit_dp_range_nz_impl(const std::vector<int>& dl, const std::vector<int>& dr, const Table& empty, Start start,
                                                   Trans trans, int base) {
    assert(!dl.empty() && dl.size() == dr.size());
    assert(base >= 2);
    for (int d : dl)
        assert(0 <= d && d < base);
    for (int d : dr)
        assert(0 <= d && d < base);
    assert(dl <= dr);

    Table cur[2][2] = {{empty, empty}, {empty, empty}};
    Table nxt[2][2] = {{empty, empty}, {empty, empty}};
    bool lz[2][2] = {{false, false}, {false, true}};
    bool nlz[2][2] = {{false, false}, {false, false}};

    long n = std::ssize(dl);
    for (long pos = 0; pos < n; ++pos) {
        for (int tl = 0; tl < 2; ++tl) {
            for (int tr = 0; tr < 2; ++tr) {
                int lo = tl ? dl[pos] : 0;
                int hi = tr ? dr[pos] : base - 1;
                for (int d = lo; d <= hi; ++d) {
                    int ntl = tl && (d == lo);
                    int ntr = tr && (d == hi);
                    call_trans(trans, cur[tl][tr], nxt[ntl][ntr], d, pos);
                }
            }
        }
        for (int tl = 0; tl < 2; ++tl) {
            for (int tr = 0; tr < 2; ++tr) {
                if (!lz[tl][tr]) continue;
                int lo = tl ? dl[pos] : 0;
                int hi = tr ? dr[pos] : base - 1;
                if (lo == 0 && 0 <= hi) {
                    int ntl = tl && (0 == lo);
                    int ntr = tr && (0 == hi);
                    nlz[ntl][ntr] = true;
                }
                for (int d = std::max(1, lo); d <= hi; ++d) {
                    int ntl = tl && (d == lo);
                    int ntr = tr && (d == hi);
                    call_start(start, nxt[ntl][ntr], d, pos);
                }
            }
        }
        for (int tl = 0; tl < 2; ++tl) {
            for (int tr = 0; tr < 2; ++tr) {
                std::swap(cur[tl][tr], nxt[tl][tr]);
                nxt[tl][tr] = empty;
                lz[tl][tr] = nlz[tl][tr];
                nlz[tl][tr] = false;
            }
        }
    }

    DigitDPRangeNZResult<Table> res;
    for (int tl = 0; tl < 2; ++tl) {
        for (int tr = 0; tr < 2; ++tr) {
            res.dp[tl][tr] = std::move(cur[tl][tr]);
            res.zero[tl][tr] = lz[tl][tr];
        }
    }
    return res;
}

}  // namespace internal

// 汎用桁DP (0 <= x <= N, leading zero を区別しない)
// 計算量: 時間 O(L * base * T_trans), 空間 O(テーブルサイズ) (L = digits.size())
// 事前条件: !digits.empty(), base >= 2, 各桁 d について 0 <= d < base
template <class Table, class F>
DigitDPResult<Table> digit_dp(const std::vector<int>& digits, const Table& init, const Table& empty, F trans, int base = 10) {
    assert(!digits.empty());
    assert(base >= 2);
    for (int d : digits)
        assert(0 <= d && d < base);

    Table cur[2] = {empty, init};
    Table nxt[2] = {empty, empty};

    long n = std::ssize(digits);
    for (long pos = 0; pos < n; ++pos) {
        for (int tight = 0; tight < 2; ++tight) {
            int lim = tight ? digits[pos] : base - 1;
            for (int d = 0; d <= lim; ++d) {
                int nt = tight && (d == lim);
                internal::call_trans(trans, cur[tight], nxt[nt], d, pos);
            }
        }
        std::swap(cur[0], nxt[0]);
        std::swap(cur[1], nxt[1]);
        nxt[0] = empty;
        nxt[1] = empty;
    }
    return {std::move(cur[0]), std::move(cur[1])};
}

template <class Table, class F>
DigitDPResult<Table> digit_dp(const std::vector<int>& digits, const Table& init, F trans, int base = 10) {
    Table empty = init;
    internal::reset_table(empty);
    return digit_dp(digits, init, empty, trans, base);
}

template <class Table, class F>
DigitDPResult<Table> digit_dp(std::string_view s, const Table& init, const Table& empty, F trans, int base = 10) {
    return digit_dp(internal::to_digits(s, base), init, empty, trans, base);
}

template <class Table, class F>
DigitDPResult<Table> digit_dp(std::string_view s, const Table& init, F trans, int base = 10) {
    return digit_dp(internal::to_digits(s, base), init, trans, base);
}

template <class Table, class F>
DigitDPResult<Table> digit_dp(long n, const Table& init, const Table& empty, F trans, int base = 10) {
    return digit_dp(internal::to_digits(n, base), init, empty, trans, base);
}

template <class Table, class F>
DigitDPResult<Table> digit_dp(long n, const Table& init, F trans, int base = 10) {
    return digit_dp(internal::to_digits(n, base), init, trans, base);
}

// Leading zero を区別する汎用桁DP (1 <= x <= N を遷移, 0 は zero_less / zero_equal で判定)
// 計算量: 時間 O(L * base * T_trans), 空間 O(テーブルサイズ) (L = digits.size())
// 事前条件: !digits.empty(), base >= 2, 各桁 d について 0 <= d < base
template <class Table, class Start, class Trans>
DigitDPNZResult<Table> digit_dp_nz(const std::vector<int>& digits, const Table& empty, Start start, Trans trans, int base = 10) {
    return internal::digit_dp_nz_impl(digits, empty, start, trans, base);
}

template <class Table, class Trans>
DigitDPNZResult<Table> digit_dp_nz(const std::vector<int>& digits, const Table& init, const Table& empty, Trans trans, int base = 10) {
    auto start = [&](Table& nxt, int d, long pos) { internal::call_trans(trans, init, nxt, d, pos); };
    return internal::digit_dp_nz_impl(digits, empty, start, trans, base);
}

template <class Table, class Trans>
DigitDPNZResult<Table> digit_dp_nz(const std::vector<int>& digits, const Table& init, Trans trans, int base = 10) {
    Table empty = init;
    internal::reset_table(empty);
    return digit_dp_nz(digits, init, empty, trans, base);
}

template <class Table, class Start, class Trans>
DigitDPNZResult<Table> digit_dp_nz(std::string_view s, const Table& empty, Start start, Trans trans, int base = 10) {
    return digit_dp_nz(internal::to_digits(s, base), empty, start, trans, base);
}

template <class Table, class Trans>
DigitDPNZResult<Table> digit_dp_nz(std::string_view s, const Table& init, const Table& empty, Trans trans, int base = 10) {
    return digit_dp_nz(internal::to_digits(s, base), init, empty, trans, base);
}

template <class Table, class Trans>
DigitDPNZResult<Table> digit_dp_nz(std::string_view s, const Table& init, Trans trans, int base = 10) {
    return digit_dp_nz(internal::to_digits(s, base), init, trans, base);
}

template <class Table, class Start, class Trans>
DigitDPNZResult<Table> digit_dp_nz(long n, const Table& empty, Start start, Trans trans, int base = 10) {
    return digit_dp_nz(internal::to_digits(n, base), empty, start, trans, base);
}

template <class Table, class Trans>
DigitDPNZResult<Table> digit_dp_nz(long n, const Table& init, const Table& empty, Trans trans, int base = 10) {
    return digit_dp_nz(internal::to_digits(n, base), init, empty, trans, base);
}

template <class Table, class Trans>
DigitDPNZResult<Table> digit_dp_nz(long n, const Table& init, Trans trans, int base = 10) {
    return digit_dp_nz(internal::to_digits(n, base), init, trans, base);
}

// 区間汎用桁DP (L <= x <= R, max/min 等減算できない問題用)
// 計算量: 時間 O(L_max * base * T_trans), 空間 O(テーブルサイズ)
// 事前条件: 0 <= L <= R, base >= 2
template <class Table, class F>
DigitDPRangeResult<Table> digit_dp_range(std::vector<int> dl, std::vector<int> dr, const Table& init, const Table& empty, F trans,
                                         int base = 10) {
    assert(!dl.empty() && !dr.empty());
    assert(base >= 2);
    for (int d : dl)
        assert(0 <= d && d < base);
    for (int d : dr)
        assert(0 <= d && d < base);
    internal::pad_digits(dl, dr);
    assert(dl <= dr);

    Table cur[2][2] = {{empty, empty}, {empty, init}};
    Table nxt[2][2] = {{empty, empty}, {empty, empty}};

    long n = std::ssize(dl);
    for (long pos = 0; pos < n; ++pos) {
        for (int tl = 0; tl < 2; ++tl) {
            for (int tr = 0; tr < 2; ++tr) {
                int lo = tl ? dl[pos] : 0;
                int hi = tr ? dr[pos] : base - 1;
                for (int d = lo; d <= hi; ++d) {
                    int ntl = tl && (d == lo);
                    int ntr = tr && (d == hi);
                    internal::call_trans(trans, cur[tl][tr], nxt[ntl][ntr], d, pos);
                }
            }
        }
        for (int tl = 0; tl < 2; ++tl) {
            for (int tr = 0; tr < 2; ++tr) {
                std::swap(cur[tl][tr], nxt[tl][tr]);
                nxt[tl][tr] = empty;
            }
        }
    }

    DigitDPRangeResult<Table> res;
    for (int tl = 0; tl < 2; ++tl) {
        for (int tr = 0; tr < 2; ++tr) {
            res.dp[tl][tr] = std::move(cur[tl][tr]);
        }
    }
    return res;
}

template <class Table, class F>
DigitDPRangeResult<Table> digit_dp_range(std::vector<int> dl, std::vector<int> dr, const Table& init, F trans, int base = 10) {
    Table empty = init;
    internal::reset_table(empty);
    return digit_dp_range(std::move(dl), std::move(dr), init, empty, trans, base);
}

template <class Table, class F>
DigitDPRangeResult<Table> digit_dp_range(std::string_view l, std::string_view r, const Table& init, const Table& empty, F trans,
                                         int base = 10) {
    return digit_dp_range(internal::to_digits(l, base), internal::to_digits(r, base), init, empty, trans, base);
}

template <class Table, class F>
DigitDPRangeResult<Table> digit_dp_range(std::string_view l, std::string_view r, const Table& init, F trans, int base = 10) {
    return digit_dp_range(internal::to_digits(l, base), internal::to_digits(r, base), init, trans, base);
}

template <class Table, class F>
DigitDPRangeResult<Table> digit_dp_range(long l, long r, const Table& init, const Table& empty, F trans, int base = 10) {
    return digit_dp_range(internal::to_digits(l, base), internal::to_digits(r, base), init, empty, trans, base);
}

template <class Table, class F>
DigitDPRangeResult<Table> digit_dp_range(long l, long r, const Table& init, F trans, int base = 10) {
    return digit_dp_range(internal::to_digits(l, base), internal::to_digits(r, base), init, trans, base);
}

// 区間汎用桁DP (L <= x <= R, Leading zero 区別版)
// 計算量: 時間 O(L_max * base * T_trans), 空間 O(テーブルサイズ)
// 事前条件: 0 <= L <= R, base >= 2
template <class Table, class Start, class Trans>
DigitDPRangeNZResult<Table> digit_dp_range_nz(std::vector<int> dl, std::vector<int> dr, const Table& empty, Start start, Trans trans,
                                              int base = 10) {
    assert(!dl.empty() && !dr.empty());
    assert(base >= 2);
    for (int d : dl)
        assert(0 <= d && d < base);
    for (int d : dr)
        assert(0 <= d && d < base);
    internal::pad_digits(dl, dr);
    assert(dl <= dr);
    return internal::digit_dp_range_nz_impl(dl, dr, empty, start, trans, base);
}

template <class Table, class Trans>
DigitDPRangeNZResult<Table> digit_dp_range_nz(std::vector<int> dl, std::vector<int> dr, const Table& init, const Table& empty, Trans trans,
                                              int base = 10) {
    auto start = [&](Table& nxt, int d, long pos) { internal::call_trans(trans, init, nxt, d, pos); };
    return digit_dp_range_nz(std::move(dl), std::move(dr), empty, start, trans, base);
}

template <class Table, class Trans>
DigitDPRangeNZResult<Table> digit_dp_range_nz(std::vector<int> dl, std::vector<int> dr, const Table& init, Trans trans, int base = 10) {
    Table empty = init;
    internal::reset_table(empty);
    return digit_dp_range_nz(std::move(dl), std::move(dr), init, empty, trans, base);
}

template <class Table, class Start, class Trans>
DigitDPRangeNZResult<Table> digit_dp_range_nz(std::string_view l, std::string_view r, const Table& empty, Start start, Trans trans,
                                              int base = 10) {
    return digit_dp_range_nz(internal::to_digits(l, base), internal::to_digits(r, base), empty, start, trans, base);
}

template <class Table, class Trans>
DigitDPRangeNZResult<Table> digit_dp_range_nz(std::string_view l, std::string_view r, const Table& init, const Table& empty, Trans trans,
                                              int base = 10) {
    return digit_dp_range_nz(internal::to_digits(l, base), internal::to_digits(r, base), init, empty, trans, base);
}

template <class Table, class Trans>
DigitDPRangeNZResult<Table> digit_dp_range_nz(std::string_view l, std::string_view r, const Table& init, Trans trans, int base = 10) {
    return digit_dp_range_nz(internal::to_digits(l, base), internal::to_digits(r, base), init, trans, base);
}

template <class Table, class Start, class Trans>
DigitDPRangeNZResult<Table> digit_dp_range_nz(long l, long r, const Table& empty, Start start, Trans trans, int base = 10) {
    return digit_dp_range_nz(internal::to_digits(l, base), internal::to_digits(r, base), empty, start, trans, base);
}

template <class Table, class Trans>
DigitDPRangeNZResult<Table> digit_dp_range_nz(long l, long r, const Table& init, const Table& empty, Trans trans, int base = 10) {
    return digit_dp_range_nz(internal::to_digits(l, base), internal::to_digits(r, base), init, empty, trans, base);
}

template <class Table, class Trans>
DigitDPRangeNZResult<Table> digit_dp_range_nz(long l, long r, const Table& init, Trans trans, int base = 10) {
    return digit_dp_range_nz(internal::to_digits(l, base), internal::to_digits(r, base), init, trans, base);
}
