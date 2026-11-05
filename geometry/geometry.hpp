#include <algorithm>
#include <cassert>
#include <compare>
#include <numeric>
#include <print>

// --- 2次元平面上の点・ベクトル ---
struct Point {
    int x, y;

    constexpr Point(int x = 0, int y = 0) : x(x), y(y) {}

    constexpr Point operator+(const Point& p) const { return {x + p.x, y + p.y}; }
    constexpr Point operator-(const Point& p) const { return {x - p.x, y - p.y}; }
    constexpr Point operator*(int k) const { return {x * k, y * k}; }

    constexpr auto operator<=>(const Point&) const = default;
};

// 内積 (dot product) -> 64bit
constexpr long dot(Point a, Point b) {
    return (long)a.x * b.x + (long)a.y * b.y;
}

// 外積 (cross product) -> 64bit
constexpr long cross(Point a, Point b) {
    return (long)a.x * b.y - (long)a.y * b.x;
}

// --- 直線: ax + by + c = 0 ---
struct Line {
    // 係数は座標の2乗オーダーになり得るため 64bit long を使用
    long a, b, c;

    Line() : a(0), b(0), c(0) {}
    Line(long a, long b, long c) : a(a), b(b), c(c) { normalize(); }

    // 2点を通る直線の方程式を構築
    Line(Point p1, Point p2) {
        assert(p1 != p2 && "2点は異なる必要があります");
        // (y2 - y1)x - (x2 - x1)y + (x2*y1 - x1*y2) = 0
        a = (long)p2.y - p1.y;
        b = (long)p1.x - p2.x;
        c = (long)p2.x * p1.y - (long)p1.x * p2.y;
        normalize();
    }

    // 方程式を一意な形に正規化 (gcdで割り、先頭の非ゼロ係数を正にする)
    // これにより、同一の直線は必ず同じ (a, b, c) を持つため map や set で管理可能
    void normalize() {
        long g = std::gcd(std::gcd(std::abs(a), std::abs(b)), std::abs(c));
        if (g != 0) {
            a /= g;
            b /= g;
            c /= g;
        }
        if (a < 0 || (a == 0 && b < 0)) {
            a = -a;
            b = -b;
            c = -c;
        }
    }

    // 点 p を代入した値: ax + by + c
    long eval(Point p) const { return a * p.x + b * p.y + c; }

    // 点 p が直線上にあるか
    bool contains(Point p) const { return eval(p) == 0; }

    // 直線同士が平行か
    bool is_parallel(const Line& l) const { return a * l.b - b * l.a == 0; }

    // 直線同士が直交か
    bool is_orthogonal(const Line& l) const { return a * l.a + b * l.b == 0; }

    auto operator<=>(const Line&) const = default;
};

// --- CCW (点の進行方向の判定) ---
// a -> b から見て c がどの位置にあるか
constexpr int CCW_COUNTER_CLOCKWISE = 1;  // 反時計回り (左折)
constexpr int CCW_CLOCKWISE = -1;         // 時計回り (右折)
constexpr int CCW_ONLINE_BACK = 2;        // c -> a -> b (同一直線上, 後方)
constexpr int CCW_ONLINE_FRONT = -2;      // a -> b -> c (同一直線上, 前方)
constexpr int CCW_ON_SEGMENT = 0;         // a -> c -> b (線分ab上)

inline int ccw(Point a, Point b, Point c) {
    long cp = cross(b - a, c - a);
    if (cp > 0) return CCW_COUNTER_CLOCKWISE;
    if (cp < 0) return CCW_CLOCKWISE;
    if (dot(b - a, c - a) < 0) return CCW_ONLINE_BACK;
    if (dot(a - b, c - b) < 0) return CCW_ONLINE_FRONT;
    return CCW_ON_SEGMENT;
}

// 線分 ab と線分 cd が交差するか判定 (端点の接触も交差に含む)
inline bool intersect(Point a, Point b, Point c, Point d) {
    return (ccw(a, b, c) * ccw(a, b, d) <= 0) && (ccw(c, d, a) * ccw(c, d, b) <= 0);
}
