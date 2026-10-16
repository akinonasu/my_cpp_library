#pragma once
#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

// Multiset of unsigned B-bit integers. xor_min returns the XOR value.
template <class UInt = std::uint32_t, int Bits = std::numeric_limits<UInt>::digits>
class BinaryTrie {
    static_assert(std::is_integral_v<UInt> && std::is_unsigned_v<UInt> && !std::is_same_v<UInt, bool>);
    static_assert(1 <= Bits && Bits <= std::numeric_limits<UInt>::digits);
    struct Node {
        int child[2] = {-1, -1};
        int count = 0;
    };
    std::vector<Node> nodes_;
    std::vector<int> free_;
    int root_ = -1;
    static void check_value([[maybe_unused]] UInt value) {
        if constexpr (Bits < std::numeric_limits<UInt>::digits) assert((value >> Bits) == 0);
    }
    int allocate() {
        if (!free_.empty()) {
            int v = free_.back();
            free_.pop_back();
            nodes_[v] = Node{};
            return v;
        }
        assert(nodes_.size() < static_cast<std::size_t>(std::numeric_limits<int>::max()));
        nodes_.push_back(Node{});
        return static_cast<int>(nodes_.size()) - 1;
    }
    int count_node(int v) const { return v == -1 ? 0 : nodes_[v].count; }

   public:
    BinaryTrie() = default;
    BinaryTrie(const BinaryTrie&) = default;
    BinaryTrie& operator=(const BinaryTrie&) = default;
    BinaryTrie(BinaryTrie&& other) noexcept
        : nodes_(std::move(other.nodes_)), free_(std::move(other.free_)), root_(std::exchange(other.root_, -1)) {
        other.nodes_.clear();
        other.free_.clear();
    }
    BinaryTrie& operator=(BinaryTrie&& other) noexcept {
        if (this != &other) {
            nodes_ = std::move(other.nodes_);
            free_ = std::move(other.free_);
            root_ = std::exchange(other.root_, -1);
            other.nodes_.clear();
            other.free_.clear();
        }
        return *this;
    }
    int size() const { return count_node(root_); }
    bool empty() const { return root_ == -1; }
    void insert(UInt value) {
        check_value(value);
        assert(size() < std::numeric_limits<int>::max());
        if (root_ == -1) root_ = allocate();
        int v = root_;
        ++nodes_[v].count;
        for (int bit = Bits - 1; bit >= 0; --bit) {
            int d = (value >> bit) & UInt{1};
            if (nodes_[v].child[d] == -1) {
                int next = allocate();
                nodes_[v].child[d] = next;
            }
            v = nodes_[v].child[d];
            ++nodes_[v].count;
        }
    }
    int count(UInt value) const {
        check_value(value);
        int v = root_;
        for (int bit = Bits - 1; bit >= 0 && v != -1; --bit)
            v = nodes_[v].child[(value >> bit) & UInt{1}];
        return count_node(v);
    }
    bool erase(UInt value) {
        check_value(value);
        if (!count(value)) return false;
        std::array<int, Bits + 1> path{};
        int v = root_;
        path[0] = v;
        for (int bit = Bits - 1; bit >= 0; --bit) {
            v = nodes_[v].child[(value >> bit) & UInt{1}];
            path[Bits - bit] = v;
        }
        for (int depth = Bits; depth >= 0; --depth) {
            v = path[depth];
            if (--nodes_[v].count == 0) {
                free_.push_back(v);
                if (depth == 0)
                    root_ = -1;
                else
                    nodes_[path[depth - 1]].child[(value >> (Bits - depth)) & UInt{1}] = -1;
            }
        }
        return true;
    }
    int rank(UInt value) const {
        check_value(value);
        int answer = 0, v = root_;
        for (int bit = Bits - 1; bit >= 0 && v != -1; --bit) {
            int d = (value >> bit) & UInt{1};
            if (d) answer += count_node(nodes_[v].child[0]);
            v = nodes_[v].child[d];
        }
        return answer;
    }
    std::optional<UInt> kth(int k) const {
        if (k < 0 || k >= size()) return std::nullopt;
        int v = root_;
        UInt answer = 0;
        for (int bit = Bits - 1; bit >= 0; --bit) {
            int left = count_node(nodes_[v].child[0]);
            int d = k >= left;
            if (d) {
                k -= left;
                answer |= UInt{1} << bit;
            }
            v = nodes_[v].child[d];
        }
        return answer;
    }
    std::optional<UInt> xor_min(UInt value) const {
        check_value(value);
        if (empty()) return std::nullopt;
        int v = root_;
        UInt answer = 0;
        for (int bit = Bits - 1; bit >= 0; --bit) {
            int d = (value >> bit) & UInt{1};
            if (nodes_[v].child[d] == -1) {
                d ^= 1;
                answer |= UInt{1} << bit;
            }
            v = nodes_[v].child[d];
        }
        return answer;
    }
    void clear() {
        root_ = -1;
        std::vector<Node>().swap(nodes_);
        std::vector<int>().swap(free_);
    }
};
