#include <array>
#include <cassert>
#include <iostream>
#include <limits>
#include <print>
#include <string>
#include <string_view>
#include <vector>

// S の空でない部分文字列のうち、T を（連続するとは限らない）部分列として含むものの個数を求める
inline long count_substrings_containing_subsequence(std::string_view s, std::string_view t) {
    assert(s.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()));
    assert(t.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()));

    for (char c : s) {
        assert('a' <= c && c <= 'z');
    }
    for (char c : t) {
        assert('a' <= c && c <= 'z');
    }

    int n = static_cast<int>(s.size());
    int m = static_cast<int>(t.size());

    if (n == 0 || m > n) {
        return 0L;
    }
    if (m == 0) {
        return static_cast<long>(n) * (n + 1) / 2;
    }

    // 各英小文字が T のどの位置に現れるかを降順で保持する
    std::array<std::vector<int>, 26> pos;
    for (int j = m - 1; j >= 0; --j) {
        pos[t[j] - 'a'].push_back(j);
    }

    // dp[j]: 現在の右端 r に対して、T[0..j] を部分列として含む最大の左端 l（不可能なら -1）
    std::vector<int> dp(m, -1);
    long ans = 0;

    for (int r = 0; r < n; ++r) {
        for (int j : pos[s[r] - 'a']) {
            dp[j] = (j == 0) ? r : dp[j - 1];
        }
        ans += static_cast<long>(dp[m - 1] + 1);
    }

    return ans;
}
