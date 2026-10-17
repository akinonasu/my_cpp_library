#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <limits>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

template <class Cost = long long>
struct Hungarian {
    static_assert(std::numeric_limits<Cost>::is_integer && std::numeric_limits<Cost>::is_signed);
    Cost cost = 0;
    std::vector<int> row_to_col, col_to_row;
    template <class Matrix>
    explicit Hungarian(const Matrix& a, bool maximize = false) {
        int n = int(a.size()), m = n ? int(a[0].size()) : 0;
        using T = std::remove_cvref_t<decltype(a[0][0])>;
        std::vector<const T*> rows(n);
        for (int i = 0; i < n; ++i) {
            assert(int(a[i].size()) == m);
            rows[i] = std::data(a[i]);
        }
        if (maximize)
            solve<true>(n, m, rows);
        else
            solve<false>(n, m, rows);
    }
    template <class T>
    Hungarian(int n, int m, const T* a, int stride = -1, bool maximize = false) {
        assert(0 <= n && n <= m);
        if (stride == -1) stride = m;
        assert(stride >= m);
        std::vector<const T*> rows(n);
        for (int i = 0; i < n; ++i)
            rows[i] = a + std::size_t(i) * stride;
        if (maximize)
            solve<true>(n, m, rows);
        else
            solve<false>(n, m, rows);
    }

   private:
    template <bool maximize, class T>
    void solve(int n, int m, const std::vector<const T*>& a) {
        assert(n <= m);
        row_to_col.assign(n, -1);
        col_to_row.assign(m, -1);
        if (!n) return;
        auto value = [&](int i, int j) -> Cost {
            if constexpr (maximize)
                return -Cost(a[i][j]);
            else
                return Cost(a[i][j]);
        };
        auto finish = [&] {
            for (int i = 0; i < n; ++i)
                cost += Cost(a[i][row_to_col[i]]);
        };
        if (n == m) {
            using Wide = std::conditional_t<(sizeof(T) < 8), long long, __int128_t>;
            int order = 3;
            for (int i = 1; i < n && order; ++i)
                for (int j = 1; j < m; ++j) {
                    Wide x = Wide(a[i - 1][j - 1]) + Wide(a[i][j]);
                    Wide y = Wide(a[i - 1][j]) + Wide(a[i][j - 1]);
                    if (x > y) order &= 2;
                    if (x < y) order &= 1;
                    if (!order) break;
                }
            if (order) {
                bool reverse = maximize ? bool(order & 1) : !(order & 1);
                for (int i = 0; i < n; ++i) {
                    int j = reverse ? n - 1 - i : i;
                    row_to_col[i] = j;
                    col_to_row[j] = i;
                }
                finish();
                return;
            }
        }
        std::vector<Cost> price(m), dist(m);
        std::vector<int> columns(m), pred(m);
        if (n == m) {
            std::vector<T> minimum(a[0], a[0] + m);
            std::vector<int> owner(m), count(n);
            for (int i = 1; i < n; ++i)
                for (int j = 0; j < m; ++j) {
                    T x = a[i][j];
                    bool better = maximize ? x > minimum[j] : x < minimum[j];
                    minimum[j] = better ? x : minimum[j];
                    owner[j] = better ? i : owner[j];
                }
            for (int j = 0; j < m; ++j)
                price[j] = maximize ? -Cost(minimum[j]) : Cost(minimum[j]);
            for (int j = m - 1; j >= 0; --j) {
                int i = owner[j];
                if (!count[i]++) row_to_col[i] = j, col_to_row[j] = i;
            }
            for (int i = 0; i < n; ++i)
                if (count[i] == 1) {
                    Cost delta = std::numeric_limits<Cost>::max();
                    int matched = row_to_col[i];
                    for (int j = 0; j < m; ++j)
                        if (j != matched) delta = std::min(delta, Cost(value(i, j) - price[j]));
                    price[matched] -= delta;
                }
            for (int pass = 0; pass < 2; ++pass)
                for (int i = 0; i < n; ++i)
                    if (row_to_col[i] == -1) {
                        int best_col = 0;
                        Cost best = value(i, 0) - price[0], second = std::numeric_limits<Cost>::max();
                        for (int j = 1; j < m; ++j) {
                            Cost x = value(i, j) - price[j];
                            if (x < best)
                                second = best, best = x, best_col = j;
                            else if (x == best && col_to_row[best_col] != -1 && col_to_row[j] == -1)
                                second = best, best_col = j;
                            else
                                second = std::min(second, x);
                        }
                        if (best < second) price[best_col] -= second - best;
                        int old = col_to_row[best_col];
                        if (old != -1) row_to_col[old] = -1;
                        row_to_col[i] = best_col;
                        col_to_row[best_col] = i;
                    }
        } else {
            for (int i = 0; i < n; ++i) {
                int best_col = 0;
                Cost best = value(i, 0);
                for (int j = 1; j < m; ++j) {
                    Cost x = value(i, j);
                    if (x < best || (x == best && col_to_row[best_col] != -1 && col_to_row[j] == -1)) best = x, best_col = j;
                }
                if (col_to_row[best_col] == -1) row_to_col[i] = best_col, col_to_row[best_col] = i;
            }
        }
        std::iota(columns.begin(), columns.end(), 0);
        for (int root = 0; root < n; ++root)
            if (row_to_col[root] == -1) {
                for (int j = 0; j < m; ++j)
                    dist[j] = value(root, j) - price[j], pred[j] = root;
                int scanned = 0, queued = 0, target = -1;
                Cost level = 0;
                while (target == -1) {
                    if (scanned == queued) {
                        level = dist[columns[scanned]];
                        for (int k = scanned; k < m; ++k) {
                            int j = columns[k];
                            if (dist[j] > level) continue;
                            if (dist[j] < level) level = dist[j], queued = scanned, target = -1;
                            std::swap(columns[k], columns[queued++]);
                            if (col_to_row[j] == -1) target = j;
                        }
                        if (target != -1) break;
                    }
                    int j0 = columns[scanned++], i = col_to_row[j0];
                    Cost offset = value(i, j0) - price[j0] - level;
                    for (int k = queued; k < m; ++k) {
                        int j = columns[k];
                        Cost candidate = value(i, j) - price[j] - offset;
                        if (candidate >= dist[j]) continue;
                        dist[j] = candidate;
                        pred[j] = i;
                        if (candidate == level) {
                            if (col_to_row[j] == -1) {
                                target = j;
                                break;
                            }
                            std::swap(columns[k], columns[queued++]);
                        }
                    }
                }
                for (int k = 0; k < scanned; ++k) {
                    int j = columns[k];
                    price[j] += dist[j] - level;
                }
                while (target != -1) {
                    int i = pred[target], next = row_to_col[i];
                    row_to_col[i] = target;
                    col_to_row[target] = i;
                    target = next;
                }
            }
        finish();
    }
};
