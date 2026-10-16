#pragma once
#include <algorithm>
#include <cassert>
#include <limits>
#include <utility>
#include <vector>

// Rectangular Hungarian algorithm: every row gets a distinct column.
template <class T>
std::pair<T, std::vector<int>> assignment(const std::vector<std::vector<T>>& cost) {
    assert(cost.size() < static_cast<std::size_t>(std::numeric_limits<int>::max()));
    const int n = static_cast<int>(cost.size());
    if (n == 0) return {T{}, {}};
    assert(cost[0].size() < static_cast<std::size_t>(std::numeric_limits<int>::max()));
    const int m = static_cast<int>(cost[0].size());
    assert(n <= m);
    for (const auto& row : cost) {
        assert(row.size() == cost[0].size());
        (void)row;
    }
    std::vector<T> row_potential(n + 1), col_potential(m + 1), slack(m + 1);
    std::vector<int> matched(m + 1), previous(m + 1);
    std::vector<bool> used(m + 1);
    for (int row = 1; row <= n; ++row) {
        matched[0] = row;
        std::fill(slack.begin(), slack.end(), std::numeric_limits<T>::max());
        std::fill(used.begin(), used.end(), false);
        int col = 0;
        do {
            used[col] = true;
            const int active = matched[col];
            T delta = std::numeric_limits<T>::max();
            int next = 0;
            for (int j = 1; j <= m; ++j)
                if (!used[j]) {
                    const T reduced = cost[active - 1][j - 1] - row_potential[active] - col_potential[j];
                    if (reduced < slack[j]) {
                        slack[j] = reduced;
                        previous[j] = col;
                    }
                    if (slack[j] < delta) {
                        delta = slack[j];
                        next = j;
                    }
                }
            assert(next != 0);
            for (int j = 0; j <= m; ++j) {
                if (used[j]) {
                    row_potential[matched[j]] += delta;
                    col_potential[j] -= delta;
                } else
                    slack[j] -= delta;
            }
            col = next;
        } while (matched[col] != 0);
        do {
            const int prev = previous[col];
            matched[col] = matched[prev];
            col = prev;
        } while (col != 0);
    }
    std::vector<int> match(n);
    for (int j = 1; j <= m; ++j)
        if (matched[j]) match[matched[j] - 1] = j - 1;
    T total{};
    for (int i = 0; i < n; ++i)
        total += cost[i][match[i]];
    return {total, std::move(match)};
}