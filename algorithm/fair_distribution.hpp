#include <algorithm>
#include <cassert>
#include <print>
#include <vector>

// 均等分配（Round-Robin 分配）を管理するクラス
struct FairDistribution {
    long m;  // 配る物の総数
    long n;  // 配る対象の人数
    long q;  // 全員に配られる基本の個数 (m / n)
    long r;  // 余り。先頭 r 人が (q + 1) 個もらう (m % n)

    FairDistribution(long total_items, long num_recipients) : m(total_items), n(num_recipients) {
        assert(n > 0);
        assert(m >= 0);
        q = m / n;
        r = m % n;
    }

    // i 番目の人 (0-indexed: 0 <= i < n) がもらう個数を取得: O(1)
    long get(long i) const {
        assert(0 <= i && i < n);
        return q + (i < r ? 1 : 0);
    }

    // 配列風アクセス: O(1)
    long operator[](long i) const { return get(i); }

    // 先頭 k 人 (0 <= k <= n) がもらう合計個数を取得: O(1)
    long prefix_sum(long k) const {
        assert(0 <= k && k <= n);
        return k * q + std::min(k, r);
    }

    // 区間 [l, r) の人がもらう合計個数を取得: O(1)
    long range_sum(long l, long r_idx) const {
        assert(0 <= l && l <= r_idx && r_idx <= n);
        return prefix_sum(r_idx) - prefix_sum(l);
    }

    // 全員の個数を std::vector として取得: O(n)
    // 32bit整数で収まる場合は to_vector<int>() と指定可能
    template <typename T = long>
    std::vector<T> to_vector() const {
        std::vector<T> res(n, static_cast<T>(q));
        std::fill_n(res.begin(), r, static_cast<T>(q + 1));
        return res;
    }
};
