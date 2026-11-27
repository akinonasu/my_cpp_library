#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

namespace internal {

#ifdef LOCAL
inline void debug_print(std::ostream& os, __int128_t x) {
    if (x == 0) {
        os << '0';
        return;
    }
    if (x < 0) {
        os << '-';
        x = -x;
    }
    std::string s;
    while (x > 0) {
        s.push_back(static_cast<char>('0' + x % 10));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    os << s;
}

template <class T>
void debug_print(std::ostream& os, const T& x) {
    os << x;
}
#endif

template <class T>
struct Fenwick {
    Fenwick() = default;
    explicit Fenwick(long n) : n_(n), data_(n, T{}) {}

    void add(long i, T x) {
        assert(0 <= i && i < n_);
        for (++i; i <= n_; i += i & -i) {
            data_[i - 1] += x;
        }
    }

    T sum(long r) const {
        assert(0 <= r && r <= n_);
        T s{};
        for (; r > 0; r -= r & -r) {
            s += data_[r - 1];
        }
        return s;
    }

    T sum(long l, long r) const {
        assert(0 <= l && l <= r && r <= n_);
        return sum(r) - sum(l);
    }

   private:
    long n_ = 0;
    std::vector<T> data_;
};

struct CoordCompress {
    std::vector<long> xs;

    void add(long x) { xs.push_back(x); }

    void build() {
        std::sort(xs.begin(), xs.end());
        xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    }

    long id(long x) const { return std::lower_bound(xs.begin(), xs.end(), x) - xs.begin(); }

    long size() const { return std::ssize(xs); }
};

}  // namespace internal

// 静的 2D 点群に対する矩形和 [xl, xr) x [yl, yr)
// 計算量: 時間 O((N + Q) log N), 空間 O(N + Q) (N: 点数, Q: クエリ数)
// 事前条件: T は +, -, T{} を持つ可換群。xl <= xr かつ yl <= yr。
template <class T = long>
struct RectangleSum {
    void add_point(long x, long y, T w = T{1}) { points_.push_back({x, y, w}); }

    long add_query(long xl, long xr, long yl, long yr) {
        assert(xl <= xr && yl <= yr);
        long id = std::ssize(queries_);
        queries_.push_back({xl, xr, yl, yr});
        return id;
    }

    std::vector<T> solve() const {
        long q = std::ssize(queries_);
        std::vector<T> ans(q, T{});
        if (points_.empty() || queries_.empty()) return ans;

        internal::CoordCompress cy;
        for (const auto& p : points_)
            cy.add(p.y);
        cy.build();

        struct Event {
            long x;
            long id;
            int sign;
            long yl_id;
            long yr_id;
        };
        std::vector<Event> events;
        events.reserve(2 * q);
        for (long i = 0; i < q; ++i) {
            const auto& qr = queries_[i];
            if (qr.xl == qr.xr || qr.yl == qr.yr) continue;
            long yl_id = cy.id(qr.yl);
            long yr_id = cy.id(qr.yr);
            if (yl_id == yr_id) continue;
            events.push_back({qr.xl, i, -1, yl_id, yr_id});
            events.push_back({qr.xr, i, 1, yl_id, yr_id});
        }
        std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) { return a.x < b.x; });

        struct P {
            long x;
            long y_id;
            T w;
        };
        std::vector<P> pts;
        pts.reserve(points_.size());
        for (const auto& p : points_) {
            pts.push_back({p.x, cy.id(p.y), p.w});
        }
        std::sort(pts.begin(), pts.end(), [](const P& a, const P& b) { return a.x < b.x; });

        internal::Fenwick<T> fw(cy.size());
        long idx = 0;
        long n = std::ssize(pts);
        for (const auto& ev : events) {
            while (idx < n && pts[idx].x < ev.x) {
                fw.add(pts[idx].y_id, pts[idx].w);
                ++idx;
            }
            T val = fw.sum(ev.yl_id, ev.yr_id);
            if (ev.sign == 1)
                ans[ev.id] += val;
            else
                ans[ev.id] -= val;
        }
        return ans;
    }

    void dump() const {
#ifdef LOCAL
        std::cerr << "RectangleSum(points=" << points_.size() << ", queries=" << queries_.size() << ")\n";
        for (long i = 0; i < std::ssize(points_); ++i) {
            std::cerr << "  pt[" << i << "]: (" << points_[i].x << ", " << points_[i].y << ") w=";
            internal::debug_print(std::cerr, points_[i].w);
            std::cerr << "\n";
        }
        for (long i = 0; i < std::ssize(queries_); ++i) {
            std::cerr << "  qr[" << i << "]: [" << queries_[i].xl << ", " << queries_[i].xr << ") x [" << queries_[i].yl << ", "
                      << queries_[i].yr << ")\n";
        }
#endif
    }

   private:
    struct Point {
        long x, y;
        T w;
    };
    struct Query {
        long xl, xr, yl, yr;
    };
    std::vector<Point> points_;
    std::vector<Query> queries_;
};

// 静的 2D 矩形加算 [xl, xr) x [yl, yr) に対する点取得
// 計算量: 時間 O((N + Q) log Q), 空間 O(N + Q) (N: 矩形数, Q: クエリ数)
// 事前条件: T は +, -, T{} を持つ可換群。xl <= xr かつ yl <= yr。
template <class T = long>
struct RectangleAddPointGet {
    void add_rect(long xl, long xr, long yl, long yr, T w = T{1}) {
        assert(xl <= xr && yl <= yr);
        if (xl == xr || yl == yr) return;
        rects_.push_back({xl, xr, yl, yr, w});
    }

    long add_query(long x, long y) {
        long id = std::ssize(queries_);
        queries_.push_back({x, y});
        return id;
    }

    std::vector<T> solve() const {
        long q = std::ssize(queries_);
        std::vector<T> ans(q, T{});
        if (rects_.empty() || queries_.empty()) return ans;

        internal::CoordCompress cy;
        for (const auto& qr : queries_)
            cy.add(qr.y);
        cy.build();
        long m = cy.size();

        struct Event {
            long x;
            long yl_id;
            long yr_id;
            T w;
        };
        std::vector<Event> events;
        events.reserve(2 * rects_.size());
        for (const auto& r : rects_) {
            long yl_id = cy.id(r.yl);
            long yr_id = cy.id(r.yr);
            if (yl_id == yr_id) continue;
            events.push_back({r.xl, yl_id, yr_id, r.w});
            events.push_back({r.xr, yl_id, yr_id, -r.w});
        }
        std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) { return a.x < b.x; });

        struct Q {
            long x;
            long y_id;
            long id;
        };
        std::vector<Q> qs;
        qs.reserve(q);
        for (long i = 0; i < q; ++i) {
            qs.push_back({queries_[i].x, cy.id(queries_[i].y), i});
        }
        std::sort(qs.begin(), qs.end(), [](const Q& a, const Q& b) { return a.x < b.x; });

        internal::Fenwick<T> fw(m);
        long idx = 0;
        long n_ev = std::ssize(events);
        for (const auto& qr : qs) {
            while (idx < n_ev && events[idx].x <= qr.x) {
                fw.add(events[idx].yl_id, events[idx].w);
                if (events[idx].yr_id < m) {
                    fw.add(events[idx].yr_id, -events[idx].w);
                }
                ++idx;
            }
            ans[qr.id] = fw.sum(qr.y_id + 1);
        }
        return ans;
    }

    void dump() const {
#ifdef LOCAL
        std::cerr << "RectangleAddPointGet(rects=" << rects_.size() << ", queries=" << queries_.size() << ")\n";
        for (long i = 0; i < std::ssize(rects_); ++i) {
            std::cerr << "  rect[" << i << "]: [" << rects_[i].xl << ", " << rects_[i].xr << ") x [" << rects_[i].yl << ", " << rects_[i].yr
                      << ") w=";
            internal::debug_print(std::cerr, rects_[i].w);
            std::cerr << "\n";
        }
        for (long i = 0; i < std::ssize(queries_); ++i) {
            std::cerr << "  qr[" << i << "]: (" << queries_[i].x << ", " << queries_[i].y << ")\n";
        }
#endif
    }

   private:
    struct Rect {
        long xl, xr, yl, yr;
        T w;
    };
    struct Query {
        long x, y;
    };
    std::vector<Rect> rects_;
    std::vector<Query> queries_;
};

// 重み付き矩形群とクエリ矩形 [xl, xr) x [yl, yr) の共通面積の重み付き和
// 計算量: 時間 O((N + Q) log N), 空間 O(N + Q) (N: 矩形数, Q: クエリ数)
// 事前条件: xl <= xr かつ yl <= yr。最終結果が型 T に収まること。
template <class T = long>
struct RectangleIntersectionAreaSum {
    void add_rect(long xl, long xr, long yl, long yr, T w = T{1}) {
        assert(xl <= xr && yl <= yr);
        if (xl == xr || yl == yr) return;
        rects_.push_back({xl, xr, yl, yr, w});
    }

    long add_query(long xl, long xr, long yl, long yr) {
        assert(xl <= xr && yl <= yr);
        long id = std::ssize(queries_);
        queries_.push_back({xl, xr, yl, yr});
        return id;
    }

    std::vector<T> solve() const {
        long q = std::ssize(queries_);
        std::vector<T> ans(q, T{});
        if (rects_.empty() || queries_.empty()) return ans;

        internal::CoordCompress cy;
        for (const auto& r : rects_) {
            cy.add(r.yl);
            cy.add(r.yr);
        }
        cy.build();

        struct Coeff {
            __int128_t xy = 0, x = 0, y = 0, c = 0;
            Coeff& operator+=(const Coeff& o) {
                xy += o.xy;
                x += o.x;
                y += o.y;
                c += o.c;
                return *this;
            }
            Coeff& operator-=(const Coeff& o) {
                xy -= o.xy;
                x -= o.x;
                y -= o.y;
                c -= o.c;
                return *this;
            }
        };

        struct Pt {
            long x;
            long y_id;
            Coeff val;
        };
        std::vector<Pt> pts;
        pts.reserve(4 * rects_.size());
        for (const auto& r : rects_) {
            __int128_t w = static_cast<__int128_t>(r.w);
            __int128_t x1 = r.xl, x2 = r.xr, y1 = r.yl, y2 = r.yr;
            long y1_id = cy.id(r.yl);
            long y2_id = cy.id(r.yr);
            pts.push_back({r.xl, y1_id, {w, -w * y1, -w * x1, w * x1 * y1}});
            pts.push_back({r.xr, y1_id, {-w, w * y1, w * x2, -w * x2 * y1}});
            pts.push_back({r.xl, y2_id, {-w, w * y2, w * x1, -w * x1 * y2}});
            pts.push_back({r.xr, y2_id, {w, -w * y2, -w * x2, w * x2 * y2}});
        }
        std::sort(pts.begin(), pts.end(), [](const Pt& a, const Pt& b) { return a.x < b.x; });

        struct Event {
            long x;
            long y;
            long y_id;
            long qid;
            int sign;
        };
        std::vector<Event> events;
        events.reserve(4 * q);
        for (long i = 0; i < q; ++i) {
            const auto& qr = queries_[i];
            if (qr.xl == qr.xr || qr.yl == qr.yr) continue;
            long yl_id = cy.id(qr.yl);
            long yr_id = cy.id(qr.yr);
            events.push_back({qr.xr, qr.yr, yr_id, i, 1});
            events.push_back({qr.xl, qr.yr, yr_id, i, -1});
            events.push_back({qr.xr, qr.yl, yl_id, i, -1});
            events.push_back({qr.xl, qr.yl, yl_id, i, 1});
        }
        std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) { return a.x < b.x; });

        std::vector<__int128_t> res(q, 0);
        internal::Fenwick<Coeff> fw(cy.size());
        long idx = 0;
        long n_pts = std::ssize(pts);
        for (const auto& ev : events) {
            while (idx < n_pts && pts[idx].x < ev.x) {
                fw.add(pts[idx].y_id, pts[idx].val);
                ++idx;
            }
            Coeff s = fw.sum(ev.y_id);
            __int128_t vx = ev.x, vy = ev.y;
            __int128_t val = s.xy * vx * vy + s.x * vx + s.y * vy + s.c;
            res[ev.qid] += ev.sign * val;
        }

        for (long i = 0; i < q; ++i) {
            ans[i] = static_cast<T>(res[i]);
        }
        return ans;
    }

    void dump() const {
#ifdef LOCAL
        std::cerr << "RectangleIntersectionAreaSum(rects=" << rects_.size() << ", queries=" << queries_.size() << ")\n";
        for (long i = 0; i < std::ssize(rects_); ++i) {
            std::cerr << "  rect[" << i << "]: [" << rects_[i].xl << ", " << rects_[i].xr << ") x [" << rects_[i].yl << ", " << rects_[i].yr
                      << ") w=";
            internal::debug_print(std::cerr, rects_[i].w);
            std::cerr << "\n";
        }
        for (long i = 0; i < std::ssize(queries_); ++i) {
            std::cerr << "  qr[" << i << "]: [" << queries_[i].xl << ", " << queries_[i].xr << ") x [" << queries_[i].yl << ", "
                      << queries_[i].yr << ")\n";
        }
#endif
    }

   private:
    struct Rect {
        long xl, xr, yl, yr;
        T w;
    };
    struct Query {
        long xl, xr, yl, yr;
    };
    std::vector<Rect> rects_;
    std::vector<Query> queries_;
};

// オフライン動的点加算・矩形和 [xl, xr) x [yl, yr)
// 計算量: 時間 O((N + Q) log^2 N), 空間 O(N log N) (N: 点加算数, Q: 矩形和数)
// 事前条件: T は +, -, T{} を持つ可換群。xl <= xr かつ yl <= yr。
template <class T = long>
struct PointAddRectangleSum {
    void add_point(long x, long y, T w = T{1}) { ops_.push_back({0, x, 0, y, 0, w}); }

    long add_query(long xl, long xr, long yl, long yr) {
        assert(xl <= xr && yl <= yr);
        long id = q_cnt_++;
        ops_.push_back({1, xl, xr, yl, yr, T{}});
        return id;
    }

    std::vector<T> solve() const {
        std::vector<T> ans;
        ans.reserve(q_cnt_);
        if (q_cnt_ == 0) return ans;

        internal::CoordCompress cx;
        for (const auto& op : ops_) {
            if (op.type == 0) cx.add(op.x1);
        }
        cx.build();
        long nx = cx.size();
        if (nx == 0) {
            ans.assign(q_cnt_, T{});
            return ans;
        }

        std::vector<std::vector<long>> ys(nx + 1);
        for (const auto& op : ops_) {
            if (op.type == 0) {
                long xi = cx.id(op.x1);
                for (long i = xi + 1; i <= nx; i += i & -i) {
                    ys[i].push_back(op.y1);
                }
            }
        }

        std::vector<internal::Fenwick<T>> fws(nx + 1);
        for (long i = 1; i <= nx; ++i) {
            std::sort(ys[i].begin(), ys[i].end());
            ys[i].erase(std::unique(ys[i].begin(), ys[i].end()), ys[i].end());
            fws[i] = internal::Fenwick<T>(std::ssize(ys[i]));
        }

        auto get_prefix = [&](long x, long yl, long yr) {
            T res{};
            long xi = cx.id(x);
            for (long i = xi; i > 0; i -= i & -i) {
                long l = std::lower_bound(ys[i].begin(), ys[i].end(), yl) - ys[i].begin();
                long r = std::lower_bound(ys[i].begin(), ys[i].end(), yr) - ys[i].begin();
                res += fws[i].sum(l, r);
            }
            return res;
        };

        for (const auto& op : ops_) {
            if (op.type == 0) {
                long xi = cx.id(op.x1);
                for (long i = xi + 1; i <= nx; i += i & -i) {
                    long yi = std::lower_bound(ys[i].begin(), ys[i].end(), op.y1) - ys[i].begin();
                    fws[i].add(yi, op.w);
                }
            } else {
                if (op.x1 == op.x2 || op.y1 == op.y2) {
                    ans.push_back(T{});
                } else {
                    ans.push_back(get_prefix(op.x2, op.y1, op.y2) - get_prefix(op.x1, op.y1, op.y2));
                }
            }
        }
        return ans;
    }

    void dump() const {
#ifdef LOCAL
        std::cerr << "PointAddRectangleSum(ops=" << ops_.size() << ", queries=" << q_cnt_ << ")\n";
        for (long i = 0; i < std::ssize(ops_); ++i) {
            const auto& op = ops_[i];
            if (op.type == 0) {
                std::cerr << "  op[" << i << "]: add (" << op.x1 << ", " << op.y1 << ") w=";
                internal::debug_print(std::cerr, op.w);
                std::cerr << "\n";
            } else {
                std::cerr << "  op[" << i << "]: sum [" << op.x1 << ", " << op.x2 << ") x [" << op.y1 << ", " << op.y2 << ")\n";
            }
        }
#endif
    }

   private:
    struct Op {
        int type;  // 0: add_point, 1: add_query
        long x1, x2, y1, y2;
        T w;
    };
    std::vector<Op> ops_;
    long q_cnt_ = 0;
};

// 複数の矩形 [xl, xr) x [yl, yr) の和集合の面積
// 計算量: 時間 O(N log N), 空間 O(N) (N: 矩形数)
// 事前条件: xl <= xr かつ yl <= yr。面積が型 T に収まること。
template <class T = long>
struct RectangleUnion {
    void add_rect(long xl, long xr, long yl, long yr) {
        assert(xl <= xr && yl <= yr);
        if (xl == xr || yl == yr) return;
        rects_.push_back({xl, xr, yl, yr});
    }

    T area() const {
        if (rects_.empty()) return T{};

        internal::CoordCompress cy;
        for (const auto& r : rects_) {
            cy.add(r.yl);
            cy.add(r.yr);
        }
        cy.build();
        long m = cy.size() - 1;

        struct Event {
            long x;
            long l;
            long r;
            int val;
        };
        std::vector<Event> events;
        events.reserve(2 * rects_.size());
        for (const auto& r : rects_) {
            long l = cy.id(r.yl);
            long r_id = cy.id(r.yr);
            events.push_back({r.xl, l, r_id, 1});
            events.push_back({r.xr, l, r_id, -1});
        }
        std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) { return a.x < b.x; });

        std::vector<int> cnt(4 * m + 4, 0);
        std::vector<T> len(4 * m + 4, T{});
        auto add = [&](auto& self, long ql, long qr, int val, long k, long l, long r) -> void {
            if (qr <= l || r <= ql) return;
            if (ql <= l && r <= qr) {
                cnt[k] += val;
            } else {
                long mid = (l + r) / 2;
                self(self, ql, qr, val, 2 * k, l, mid);
                self(self, ql, qr, val, 2 * k + 1, mid, r);
            }
            if (cnt[k] > 0) {
                len[k] = static_cast<T>(cy.xs[r]) - static_cast<T>(cy.xs[l]);
            } else if (r - l == 1) {
                len[k] = T{};
            } else {
                len[k] = len[2 * k] + len[2 * k + 1];
            }
        };

        T ans{};
        long prev_x = events[0].x;
        for (const auto& ev : events) {
            if (ev.x > prev_x) {
                ans += (static_cast<T>(ev.x) - static_cast<T>(prev_x)) * len[1];
                prev_x = ev.x;
            }
            add(add, ev.l, ev.r, ev.val, 1, 0, m);
        }
        return ans;
    }

    void dump() const {
#ifdef LOCAL
        std::cerr << "RectangleUnion(rects=" << rects_.size() << ")\n";
        for (long i = 0; i < std::ssize(rects_); ++i) {
            std::cerr << "  rect[" << i << "]: [" << rects_[i].xl << ", " << rects_[i].xr << ") x [" << rects_[i].yl << ", " << rects_[i].yr
                      << ")\n";
        }
#endif
    }

   private:
    struct Rect {
        long xl, xr, yl, yr;
    };
    std::vector<Rect> rects_;
};

// 重み付き矩形 [xl, xr) x [yl, yr) の重なり和の最大値
// 計算量: 時間 O(N log N), 空間 O(N) (N: 矩形数)
// 事前条件: T は +, <, T{} を持つ型。xl <= xr かつ yl <= yr。
template <class T = long>
struct RectangleMaxOverlap {
    void add_rect(long xl, long xr, long yl, long yr, T w = T{1}) {
        assert(xl <= xr && yl <= yr);
        if (xl == xr || yl == yr) return;
        rects_.push_back({xl, xr, yl, yr, w});
    }

    T max_overlap() const {
        if (rects_.empty()) return T{};

        internal::CoordCompress cy;
        for (const auto& r : rects_) {
            cy.add(r.yl);
            cy.add(r.yr);
        }
        cy.build();
        long m = cy.size() - 1;

        struct Event {
            long x;
            long l;
            long r;
            T w;
        };
        std::vector<Event> events;
        events.reserve(2 * rects_.size());
        for (const auto& r : rects_) {
            long l = cy.id(r.yl);
            long r_id = cy.id(r.yr);
            events.push_back({r.xl, l, r_id, r.w});
            events.push_back({r.xr, l, r_id, -r.w});
        }
        std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) { return a.x < b.x; });

        std::vector<T> mx(4 * m + 4, T{}), lazy(4 * m + 4, T{});
        auto push = [&](long k) {
            if (lazy[k] != T{}) {
                for (long c : {2 * k, 2 * k + 1}) {
                    mx[c] += lazy[k];
                    lazy[c] += lazy[k];
                }
                lazy[k] = T{};
            }
        };
        auto add = [&](auto& self, long ql, long qr, T val, long k, long l, long r) -> void {
            if (qr <= l || r <= ql) return;
            if (ql <= l && r <= qr) {
                mx[k] += val;
                lazy[k] += val;
                return;
            }
            push(k);
            long mid = (l + r) / 2;
            self(self, ql, qr, val, 2 * k, l, mid);
            self(self, ql, qr, val, 2 * k + 1, mid, r);
            mx[k] = std::max(mx[2 * k], mx[2 * k + 1]);
        };

        T ans{};
        long n_ev = std::ssize(events);
        for (long i = 0; i < n_ev;) {
            long x = events[i].x;
            while (i < n_ev && events[i].x == x) {
                add(add, events[i].l, events[i].r, events[i].w, 1, 0, m);
                ++i;
            }
            ans = std::max(ans, mx[1]);
        }
        return ans;
    }

    void dump() const {
#ifdef LOCAL
        std::cerr << "RectangleMaxOverlap(rects=" << rects_.size() << ")\n";
        for (long i = 0; i < std::ssize(rects_); ++i) {
            std::cerr << "  rect[" << i << "]: [" << rects_[i].xl << ", " << rects_[i].xr << ") x [" << rects_[i].yl << ", " << rects_[i].yr
                      << ") w=";
            internal::debug_print(std::cerr, rects_[i].w);
            std::cerr << "\n";
        }
#endif
    }

   private:
    struct Rect {
        long xl, xr, yl, yr;
        T w;
    };
    std::vector<Rect> rects_;
};
