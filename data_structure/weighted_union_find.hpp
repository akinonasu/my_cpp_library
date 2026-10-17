#include <cassert>
#include <cstdint>
#include <utility>
#include <vector>

template <class T = std::int64_t>
struct WeightedUnionFind {
    std::vector<int> parent;
    std::vector<T> delta;
    explicit WeightedUnionFind(int n) {
        assert(n >= 0);
        parent.assign(n, -1);
        delta.resize(n);
    }
    int leader(int x) {
        assert(0 <= x && x < static_cast<int>(parent.size()));
        if (parent[x] < 0) return x;
        int r = leader(parent[x]);
        delta[x] += delta[parent[x]];
        return parent[x] = r;
    }
    T weight(int x) {
        leader(x);
        return delta[x];
    }
    bool same(int x, int y) { return leader(x) == leader(y); }
    int size(int x) { return -parent[leader(x)]; }
    T diff(int x, int y) {
        assert(same(x, y));
        return weight(y) - weight(x);
    }
    bool merge(int x, int y, T w) {
        int rx = leader(x), ry = leader(y);
        w += delta[x] - delta[y];
        if (rx == ry) return w == T(0);
        if (parent[rx] > parent[ry]) std::swap(rx, ry), w = -w;
        parent[rx] += parent[ry];
        parent[ry] = rx;
        delta[ry] = w;
        return true;
    }
};
