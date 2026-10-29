#include <print>
// 拡張ユークリッドの互助法の結果を格納する構造体
struct ExtGCDResult {
    long gcd;
    long x;
    long y;
};

/**
 * @brief 拡張ユークリッドの互助法
 * a * x + b * y = gcd(a, b) を満たす整数解 (x, y) および gcd(a, b) を求めます。
 * gcd は常に非負になるよう正規化されます。
 *
 * @param a 整数1
 * @param b 整数2
 * @return ExtGCDResult {gcd, x, y}
 */
constexpr ExtGCDResult ext_gcd(long a, long b) {
    long x0 = 1, x1 = 0;
    long y0 = 0, y1 = 1;
    while (b != 0) {
        long q = a / b;
        long r = a % b;
        a = b;
        b = r;

        long next_x = x0 - q * x1;
        x0 = x1;
        x1 = next_x;

        long next_y = y0 - q * y1;
        y0 = y1;
        y1 = next_y;
    }
    // gcd を正の値に正規化
    if (a < 0) {
        a = -a;
        x0 = -x0;
        y0 = -y0;
    }
    return {a, x0, y0};
}

/**
 * @brief mod m における a の乗法逆元を求めます (a * x ≡ 1 (mod m))
 * gcd(a, m) == 1 である必要があります。
 *
 * @param a
 * @param m 法
 * @return long 逆元 (0 <= x < m)。存在しない場合は -1
 */
constexpr long mod_inv(long a, long m) {
    auto [g, x, y] = ext_gcd(a, m);
    if (g != 1) return -1;  // 互いに素でない場合は逆元が存在しない
    long res = x % m;
    if (res < 0) res += m;
    return res;
}
