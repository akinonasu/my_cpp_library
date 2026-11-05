#include <chrono>
#include <concepts>
#include <print>
#include <ranges>
#include <vector>

class FastRandom {
   public:
    using result_type = unsigned long;

    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return ~0UL; }

   private:
    unsigned long s[4];

    static constexpr unsigned long rotl(const unsigned long x, int k) noexcept { return (x << k) | (x >> (64 - k)); }

    // シード初期化用のSplitMix64
    static unsigned long splitmix64(unsigned long& x) noexcept {
        unsigned long z = (x += 0x9e3779b97f4a7c15UL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9UL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebUL;
        return z ^ (z >> 31);
    }

    // Lemire's nearly divisionless method: [0, n)
    unsigned long next_bound(unsigned long n) noexcept {
        __uint128_t m = static_cast<__uint128_t>(next()) * n;
        unsigned long l = static_cast<unsigned long>(m);
        if (l < n) {
            unsigned long t = -n % n;
            while (l < t) {
                m = static_cast<__uint128_t>(next()) * n;
                l = static_cast<unsigned long>(m);
            }
        }
        return static_cast<unsigned long>(m >> 64);
    }

   public:
    // 時刻とスタックアドレスを組み合わせてハックを防ぐ安全なシード初期化
    FastRandom() noexcept {
        unsigned long seed_val = std::chrono::steady_clock::now().time_since_epoch().count();
        seed_val ^= reinterpret_cast<unsigned long>(&seed_val);
        seed(seed_val);
    }

    explicit FastRandom(unsigned long seed_val) noexcept { seed(seed_val); }

    void seed(unsigned long seed_val) noexcept {
        s[0] = splitmix64(seed_val);
        s[1] = splitmix64(seed_val);
        s[2] = splitmix64(seed_val);
        s[3] = splitmix64(seed_val);
    }

    // 64bit乱数生値 (Xoshiro256**)
    unsigned long next() noexcept {
        const unsigned long result = rotl(s[1] * 5, 7) * 9;
        const unsigned long t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }

    unsigned long operator()() noexcept { return next(); }

    // [0, r] の整数
    template <std::integral T>
    T operator()(T r) noexcept {
        return operator()(T(0), r);
    }

    // [l, r] の整数
    template <std::integral T>
    T operator()(T l, T r) noexcept {
        unsigned long length = static_cast<unsigned long>(r - l) + 1UL;
        if (length == 0) return static_cast<T>(next());  // 全探索範囲のオーバーフロー対策
        return static_cast<T>(l + next_bound(length));
    }

    // [l, r) の浮動小数点数 (53bit精度)
    double real(double l = 0.0, double r = 1.0) noexcept {
        double val = (next() >> 11) * (1.0 / 9007199254740992.0);
        return l + val * (r - l);
    }

    // コンテナ・配列のシャッフル
    template <std::ranges::random_access_range R>
    void shuffle(R&& range) noexcept {
        int n = static_cast<int>(std::ranges::size(range));
        for (int i = n - 1; i > 0; --i) {
            int j = operator()(0, i);
            std::swap(range[i], range[j]);
        }
    }

    // コンテナからランダムに1要素取得
    template <std::ranges::random_access_range R>
    auto choice(const R& range) noexcept -> decltype(range[0]) {
        return range[operator()(0, static_cast<int>(std::ranges::size(range)) - 1)];
    }
};

// グローバルインスタンス (競プロ用)
inline FastRandom rng;
