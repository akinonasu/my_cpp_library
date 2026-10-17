#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

template <class S, auto op, auto e, class Coord = std::int64_t>
class DynamicSegtree {
    static_assert(std::is_integral_v<Coord> && !std::is_same_v<Coord, bool> && sizeof(Coord) <= 8);
    static constexpr unsigned B = 8, limit = (1u << 26) - 1;
    using U = std::make_unsigned_t<Coord>;
    static constexpr unsigned tag = 1u << 31;
    static constexpr U sign = std::is_signed_v<Coord> ? U(1) << (std::numeric_limits<U>::digits - 1) : U(0);
    template <std::size_t... I>
    static auto units(std::index_sequence<I...>) {
        return std::array<S, sizeof...(I)>{(void(I), e())...};
    }
    struct Leaf {
        std::array<U, B> keys{};
        std::array<S, B> values = units(std::make_index_sequence<B>{});
        S product = e();
        unsigned count = 0;
    };
    struct Node {
        U middle = 0;
        unsigned left = 0, right = 0;
        S product = e();
    };
    std::vector<Leaf> leaves;
    std::vector<Node> nodes{Node{}};
    unsigned root = 0, points = 0;
    Coord first, last;
    static U encode(Coord x) { return U(x) ^ sign; }
    static Coord decode(U x) { return Coord(x ^ sign); }
    Leaf& leaf(unsigned id) { return leaves[id & ~tag]; }
    const Leaf& leaf(unsigned id) const { return leaves[id & ~tag]; }
    const S& product(unsigned id) const { return id & tag ? leaf(id).product : nodes[id].product; }
    unsigned make_leaf() {
        leaves.emplace_back();
        return unsigned(leaves.size() - 1) | tag;
    }
    unsigned make_node(U middle, unsigned left, unsigned right) {
        S value = op(product(left), product(right));
        nodes.push_back({middle, left, right, std::move(value)});
        return unsigned(nodes.size() - 1);
    }
    template <unsigned N>
    [[gnu::always_inline]] static inline S fold(const S* values, unsigned count) {
        if constexpr (N == 1)
            return values[0];
        else {
            if (count <= N / 2) return fold<N / 2>(values, count);
            return op(fold<N / 2>(values, N / 2), fold<N / 2>(values + N / 2, count - N / 2));
        }
    }
    static void pull(Leaf& node) { node.product = fold<B>(node.values.data(), node.count); }

    template <bool Append>
    static void write(Leaf& node, U key, S value) {
        unsigned p = std::lower_bound(node.keys.begin(), node.keys.begin() + node.count, key) - node.keys.begin();
        if (p < node.count && node.keys[p] == key) {
            if constexpr (Append)
                node.values[p] = op(node.values[p], value);
            else
                node.values[p] = std::move(value);
        } else {
            std::move_backward(node.keys.begin() + p, node.keys.begin() + node.count, node.keys.begin() + node.count + 1);
            std::move_backward(node.values.begin() + p, node.values.begin() + node.count, node.values.begin() + node.count + 1);
            node.keys[p] = key;
            node.values[p] = std::move(value);
            ++node.count;
        }
        pull(node);
    }
    template <bool Append>
    void update(U key, S value) {
        if (!root) {
            root = make_leaf();
            write<Append>(leaf(root), key, std::move(value));
            ++points;
            return;
        }
        std::array<unsigned, 64> path;
        unsigned depth = 0, id = root;
        bool outside = false;
        while (!(id & tag)) {
            unsigned h = std::countr_zero(nodes[id].middle);
            if (((key ^ nodes[id].middle) >> h) > 1) {
                outside = true;
                break;
            }
            assert(depth < path.size());
            path[depth++] = id;
            id = key < nodes[id].middle ? nodes[id].left : nodes[id].right;
        }
        unsigned replacement = id;
        if (outside) {
            assert(points < limit);
            unsigned added = make_leaf();
            write<Append>(leaf(added), key, std::move(value));
            unsigned h = std::bit_width(U(key ^ nodes[id].middle)) - 1;
            U middle = ((key | nodes[id].middle) >> h) << h;
            replacement = key < middle ? make_node(middle, added, id) : make_node(middle, id, added);
            ++points;
        } else {
            auto& node = leaf(id);
            unsigned p = std::lower_bound(node.keys.begin(), node.keys.begin() + node.count, key) - node.keys.begin();
            bool exists = p < node.count && node.keys[p] == key;
            assert(exists || points < limit);
            if (exists || node.count < B)
                write<Append>(node, key, std::move(value));
            else {
                U low = std::min(key, node.keys[0]), high = std::max(key, node.keys[B - 1]);
                unsigned h = std::bit_width(U(low ^ high)) - 1;
                U middle = ((low | high) >> h) << h;
                unsigned cut = std::lower_bound(node.keys.begin(), node.keys.end(), middle) - node.keys.begin();
                unsigned added = make_leaf();
                auto &a = leaf(id), &b = leaf(added);
                b.count = B - cut;
                for (unsigned i = 0; i < b.count; ++i)
                    b.keys[i] = a.keys[cut + i], b.values[i] = std::move(a.values[cut + i]);
                a.count = cut;
                if (key < middle)
                    write<Append>(a, key, std::move(value));
                else
                    write<Append>(b, key, std::move(value));
                if (key < middle)
                    pull(b);
                else
                    pull(a);
                replacement = make_node(middle, id, added);
            }
            points += !exists;
        }
        if (depth) {
            auto parent = path[depth - 1];
            if (key < nodes[parent].middle)
                nodes[parent].left = replacement;
            else
                nodes[parent].right = replacement;
        } else
            root = replacement;
        while (depth) {
            unsigned k = path[--depth];
            nodes[k].product = op(product(nodes[k].left), product(nodes[k].right));
        }
    }
    S prefix(unsigned id, U r) const {
        S result = e();
        while (!(id & tag)) {
            if (r <= nodes[id].middle)
                id = nodes[id].left;
            else {
                result = op(result, product(nodes[id].left));
                id = nodes[id].right;
            }
        }
        const auto& node = leaf(id);
        for (unsigned i = 0; i < node.count && node.keys[i] < r; ++i)
            result = op(result, node.values[i]);
        return result;
    }
    S suffix(unsigned id, U l) const {
        S result = e();
        while (!(id & tag)) {
            if (nodes[id].middle <= l)
                id = nodes[id].right;
            else {
                result = op(product(nodes[id].right), result);
                id = nodes[id].left;
            }
        }
        const auto& node = leaf(id);
        for (unsigned i = node.count; i && l <= node.keys[i - 1];)
            result = op(node.values[--i], result);
        return result;
    }
    template <class Pred>
    bool search_right(unsigned id, U lower, U l, S& value, Pred& f, Coord& answer) const {
        if (l <= lower) {
            S next = op(value, product(id));
            if (f(next)) {
                value = std::move(next);
                return false;
            }
        }
        if (id & tag) {
            const auto& node = leaf(id);
            unsigned p = std::lower_bound(node.keys.begin(), node.keys.begin() + node.count, l) - node.keys.begin();
            for (; p < node.count; ++p) {
                S next = op(value, node.values[p]);
                if (!f(next)) {
                    answer = decode(node.keys[p]);
                    return true;
                }
                value = std::move(next);
            }
            return false;
        }
        U middle = nodes[id].middle;
        if (l < middle && search_right(nodes[id].left, lower, l, value, f, answer)) return true;
        return search_right(nodes[id].right, middle, l, value, f, answer);
    }
    template <class Pred>
    bool search_left(unsigned id, U upper, U r, S& value, Pred& f, Coord& answer) const {
        if (upper <= r) {
            S next = op(product(id), value);
            if (f(next)) {
                value = std::move(next);
                return false;
            }
        }
        if (id & tag) {
            const auto& node = leaf(id);
            unsigned p = std::lower_bound(node.keys.begin(), node.keys.begin() + node.count, r) - node.keys.begin();
            while (p) {
                S next = op(node.values[--p], value);
                if (!f(next)) {
                    answer = decode(node.keys[p]) + 1;
                    return true;
                }
                value = std::move(next);
            }
            return false;
        }
        U middle = nodes[id].middle;
        if (middle < r && search_left(nodes[id].right, upper, r, value, f, answer)) return true;
        return search_left(nodes[id].left, middle, r, value, f, answer);
    }

   public:
    using value_type = S;
    using coordinate_type = Coord;
    explicit DynamicSegtree(Coord n = 0) : DynamicSegtree(0, n) {}
    DynamicSegtree(Coord l, Coord r) : first(l), last(r) { assert(l <= r); }
    DynamicSegtree(const DynamicSegtree&) = default;
    DynamicSegtree& operator=(const DynamicSegtree&) = default;
    DynamicSegtree(DynamicSegtree&& other) : DynamicSegtree(other.first, other.last) { swap(other); }
    DynamicSegtree& operator=(DynamicSegtree&& other) {
        if (this != &other) {
            DynamicSegtree moved(std::move(other));
            swap(moved);
        }
        return *this;
    }
    void swap(DynamicSegtree& other) noexcept {
        leaves.swap(other.leaves);
        nodes.swap(other.nodes);
        std::swap(root, other.root);
        std::swap(points, other.points);
        std::swap(first, other.first);
        std::swap(last, other.last);
    }
    void reserve(std::size_t n) {
        assert(n <= limit);
        std::size_t capacity = (3 * n + 2 * B - 1) / (2 * B);
        leaves.reserve(capacity);
        nodes.reserve(capacity + 1);
    }
    std::size_t node_count() const { return points; }
    std::size_t memory_usage() const { return leaves.capacity() * sizeof(Leaf) + nodes.capacity() * sizeof(Node); }
    U size() const { return U(last) - U(first); }
    Coord left_bound() const { return first; }
    Coord right_bound() const { return last; }
    void clear() {
        leaves.clear();
        nodes.erase(nodes.begin() + 1, nodes.end());
        root = points = 0;
    }
    void set(Coord p, S value) {
        assert(first <= p && p < last);
        update<false>(encode(p), std::move(value));
    }
    void multiply(Coord p, S value) {
        assert(first <= p && p < last);
        update<true>(encode(p), std::move(value));
    }
    S get(Coord p) const {
        assert(first <= p && p < last);
        if (!root) return e();
        U key = encode(p);
        unsigned id = root;
        while (!(id & tag))
            id = key < nodes[id].middle ? nodes[id].left : nodes[id].right;
        const auto& node = leaf(id);
        unsigned i = std::lower_bound(node.keys.begin(), node.keys.begin() + node.count, key) - node.keys.begin();
        return i < node.count && node.keys[i] == key ? node.values[i] : e();
    }
    S prod(Coord l, Coord r) const {
        assert(first <= l && l <= r && r <= last);
        if (!root || l == r) return e();
        if (l == first && r == last) return product(root);
        U a = encode(l), b = encode(r);
        unsigned id = root;
        while (!(id & tag)) {
            U middle = nodes[id].middle;
            if (b <= middle)
                id = nodes[id].left;
            else if (middle <= a)
                id = nodes[id].right;
            else
                return op(suffix(nodes[id].left, a), prefix(nodes[id].right, b));
        }
        const auto& node = leaf(id);
        unsigned i = std::lower_bound(node.keys.begin(), node.keys.begin() + node.count, a) - node.keys.begin();
        S result = e();
        for (; i < node.count && node.keys[i] < b; ++i)
            result = op(result, node.values[i]);
        return result;
    }
    S all_prod() const { return root ? product(root) : e(); }
    template <class Pred>
    Coord max_right(Coord l, Pred f) const {
        assert(first <= l && l <= last && f(e()));
        Coord answer = last;
        if (root && l != last) {
            S value = e();
            search_right(root, encode(first), encode(l), value, f, answer);
        }
        return answer;
    }
    template <class Pred>
    Coord min_left(Coord r, Pred f) const {
        assert(first <= r && r <= last && f(e()));
        Coord answer = first;
        if (root && r != first) {
            S value = e();
            search_left(root, encode(last), encode(r), value, f, answer);
        }
        return answer;
    }
};
