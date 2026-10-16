#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <limits>
#include <string_view>
#include <vector>

// Fixed byte alphabet [First, First + Alphabet). Insert all patterns before build.
template <int Alphabet = 26, unsigned char First = 'a'>
class AhoCorasick {
    static_assert(Alphabet > 0 && static_cast<int>(First) + Alphabet <= 256);
    struct Node {
        std::array<int, Alphabet> next;
        int parent, link = 0;
        explicit Node(int p) : parent(p) { next.fill(-1); }
    };
    std::vector<Node> nodes_;
    std::vector<int> order_;
    bool built_ = false;

    static int symbol(char c) {
        const int x = static_cast<unsigned char>(c) - static_cast<int>(First);
        assert(0 <= x && x < Alphabet);
        return x;
    }

   public:
    explicit AhoCorasick(int reserve_nodes = 0) {
        assert(reserve_nodes >= 0);
        nodes_.reserve(static_cast<std::size_t>(reserve_nodes));
        nodes_.emplace_back(-1);
    }
    int size() const { return static_cast<int>(nodes_.size()); }

    int add(std::string_view pattern) {
        assert(!built_);
        int v = 0;
        for (char c : pattern) {
            const int x = symbol(c);
            if (nodes_[v].next[x] == -1) {
                assert(nodes_.size() < static_cast<std::size_t>(std::numeric_limits<int>::max()));
                const int u = size();
                nodes_.emplace_back(v);
                nodes_[v].next[x] = u;
            }
            v = nodes_[v].next[x];
        }
        return v;
    }

    void build() {
        if (built_) return;
        order_.reserve(nodes_.size());
        order_.push_back(0);
        for (int x = 0; x < Alphabet; ++x) {
            const int child = nodes_[0].next[x];
            if (child == -1)
                nodes_[0].next[x] = 0;
            else
                order_.push_back(child);
        }
        for (std::size_t head = 1; head < order_.size(); ++head) {
            const int v = order_[head], failure = nodes_[v].link;
            for (int x = 0; x < Alphabet; ++x) {
                const int child = nodes_[v].next[x];
                if (child == -1)
                    nodes_[v].next[x] = nodes_[failure].next[x];
                else {
                    nodes_[child].link = nodes_[failure].next[x];
                    order_.push_back(child);
                }
            }
        }
        built_ = true;
    }

    int parent(int v) const {
        assert(0 <= v && v < size());
        return nodes_[v].parent;
    }
    int link(int v) const {
        assert(built_ && 0 <= v && v < size());
        return nodes_[v].link;
    }
    int next(int v, char c) const {
        assert(built_ && 0 <= v && v < size());
        return nodes_[v].next[symbol(c)];
    }

    // Overlapping occurrences for every trie node, including |text|+1 at root.
    std::vector<long long> count(std::string_view text) const {
        assert(built_);
        assert(text.size() < static_cast<std::size_t>(std::numeric_limits<long long>::max()));
        std::vector<long long> result(nodes_.size());
        int v = 0;
        for (char c : text) {
            v = next(v, c);
            ++result[v];
        }
        for (std::size_t i = order_.size(); i > 1; --i) {
            const int u = order_[i - 1];
            result[nodes_[u].link] += result[u];
        }
        ++result[0];
        return result;
    }
};
