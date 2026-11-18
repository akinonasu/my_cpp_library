#pragma once

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

struct TileRun {
    long start;  // 最初のタイルの左端座標
    long count;  // start, start + K, ..., start + (count - 1) * K と連続配置する個数
};

// 区間列 [l, r] をマージし、ソート済みかつ互いに素な区間列を返す
template <typename T>
inline std::vector<std::pair<T, T>> merge_intervals(std::vector<std::pair<T, T>> intervals) {
    for (const auto& [l, r] : intervals) {
        assert(l <= r);
    }
    if (intervals.empty()) {
        return {};
    }
    std::sort(intervals.begin(), intervals.end());

    std::vector<std::pair<T, T>> merged;
    merged.reserve(intervals.size());
    merged.push_back(intervals[0]);
    for (std::size_t i = 1; i < intervals.size(); ++i) {
        const auto& [l, r] = intervals[i];
        if (l <= merged.back().second) {
            merged.back().second = std::max(merged.back().second, r);
        } else {
            merged.emplace_back(l, r);
        }
    }
    return merged;
}

// 既存の被覆区間 existing がある状態で、長さ K のタイルを追加して [L, R] を覆う配置計画を返す
// （L == R のときは空を返す。各既存区間が長さ K 未満でタイルが跨ぐケースにも対応）
inline std::vector<TileRun> build_tiles_to_cover(long L, long R, long K, const std::vector<std::pair<long, long>>& existing = {}) {
    assert(L <= R);
    assert(K > 0);
    for (const auto& [l, r] : existing) {
        assert(l <= r);
    }
    if (L == R) {
        return {};
    }

    std::vector<std::pair<long, long>> clipped;
    clipped.reserve(existing.size());
    for (const auto& [l, r] : existing) {
        long cl = std::max(l, L);
        long cr = std::min(r, R);
        if (cl < cr) {
            clipped.emplace_back(cl, cr);
        }
    }
    auto merged = merge_intervals(std::move(clipped));

    std::vector<TileRun> runs;
    long cur = L;
    std::size_t idx = 0;

    while (cur < R) {
        while (idx < merged.size() && merged[idx].second <= cur) {
            ++idx;
        }
        if (idx < merged.size() && merged[idx].first <= cur) {
            cur = merged[idx].second;
            ++idx;
            continue;
        }

        long target = (idx < merged.size()) ? merged[idx].first : R;
        long gap = target - cur;
        long cnt = (gap - 1) / K + 1;
        runs.push_back({cur, cnt});

        long rem = R - cur;
        long max_cnt = (rem - 1) / K + 1;
        if (cnt == max_cnt) {
            break;
        }
        cur += cnt * K;
    }
    return runs;
}

// 既存の被覆区間 existing がある状態で、長さ K のタイルを追加して [L, R] を覆う最小個数を返す
inline long min_tiles_to_cover(long L, long R, long K, const std::vector<std::pair<long, long>>& existing = {}) {
    auto runs = build_tiles_to_cover(L, R, K, existing);
    long total = 0;
    for (const auto& run : runs) {
        total += run.count;
    }
    return total;
}

// 候補区間リスト intervals から最小個数を選んで [L, R] を覆う
// 戻り値: {成功したか, 選んだ区間の元のインデックス(0-indexed)}
inline std::pair<bool, std::vector<int>> min_interval_cover(long L, long R, const std::vector<std::pair<long, long>>& intervals) {
    assert(L <= R);
    for (const auto& [l, r] : intervals) {
        assert(l <= r);
    }
    if (L == R) {
        return {true, {}};
    }

    struct Item {
        long l, r;
        int id;
    };
    int n = static_cast<int>(intervals.size());
    std::vector<Item> items;
    items.reserve(n);
    for (int i = 0; i < n; ++i) {
        items.push_back({intervals[i].first, intervals[i].second, i});
    }
    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        if (a.l != b.l) return a.l < b.l;
        return a.r > b.r;
    });

    std::vector<int> chosen;
    long cur = L;
    int idx = 0;
    while (cur < R) {
        long best_r = cur;
        int best_id = -1;
        while (idx < n && items[idx].l <= cur) {
            if (items[idx].r > best_r) {
                best_r = items[idx].r;
                best_id = items[idx].id;
            }
            ++idx;
        }
        if (best_id == -1) {
            return {false, {}};
        }
        chosen.push_back(best_id);
        cur = best_r;
    }
    return {true, chosen};
}

// 候補区間リスト intervals から最小個数を選んで [L, R] を覆うときの最小個数（不可能なら -1）
inline int min_interval_cover_count(long L, long R, const std::vector<std::pair<long, long>>& intervals) {
    auto [ok, chosen] = min_interval_cover(L, R, intervals);
    if (!ok) {
        return -1;
    }
    return static_cast<int>(chosen.size());
}
