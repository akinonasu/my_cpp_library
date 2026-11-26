#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "heuristic.hpp"

void solve_tsp() {
    constexpr int n = 100;
    Random init_rng(12345);

    std::vector<long double> x(n), y(n);
    for (int i = 0; i < n; ++i) {
        x[i] = init_rng.next_real(0.0L, 1000.0L);
        y[i] = init_rng.next_real(0.0L, 1000.0L);
    }

    // 2点間距離を事前計算
    std::vector<long double> dist(static_cast<std::size_t>(n) * n);
    auto d = [&](int i, int j) -> long double { return dist[static_cast<long>(i) * n + j]; };
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            long double dx = x[i] - x[j];
            long double dy = y[i] - y[j];
            dist[static_cast<long>(i) * n + j] = std::sqrt(dx * dx + dy * dy);
        }
    }

    // 初期解 (0 -> 1 -> ... -> n-1 -> 0)
    std::vector<int> order(n);
    for (int i = 0; i < n; ++i) {
        order[i] = i;
    }
    long double cur_len = 0.0L;
    for (int i = 0; i < n; ++i) {
        cur_len += d(order[i], order[(i + 1) % n]);
    }

    std::vector<int> best_order = order;
    long double best_len = cur_len;

    std::cout << "TSP init length: " << cur_len << "\n";

    // 開始温度 50.0 -> 終了温度 0.1、制限時間 500 ms (指数スケジュール)
    Annealer sa(50.0L, 0.1L, 500, TempSchedule::Exp, 42);

    while (sa.step()) {
        // 0 を固定し、半開区間 [l, r) (長さ 2 以上) を反転する 2-opt 近傍
        int l = sa.rng.next_int(1, n - 1);
        int r = sa.rng.next_int(l + 2, n + 1);
        if (l == 1 && r == n) {
            continue;
        }

        int a = order[l - 1];
        int b = order[l];
        int c = order[r - 1];
        int e = order[r == n ? 0 : r];

        long double delta = (d(a, c) + d(b, e)) - (d(a, b) + d(c, e));

        // 最小化問題なので accept_min を使用
        if (sa.accept_min(delta)) {
            std::reverse(order.begin() + l, order.begin() + r);
            cur_len += delta;
            if (cur_len < best_len) {
                best_len = cur_len;
                best_order = order;
            }
        }
    }

    sa.dump();  // -DLOCAL 時のみ反復数・採択率などを出力
    std::cout << "TSP best length: " << best_len << "\n";
}

int main() {
    solve_tsp();
    std::cout << "----\n";
}