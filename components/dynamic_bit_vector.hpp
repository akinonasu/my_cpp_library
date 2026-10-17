#ifndef DYNAMIC_BIT_VECTOR_HPP
#define DYNAMIC_BIT_VECTOR_HPP

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

template <int BlockBits = 1024, int Branch = 8>
class DynamicBitVector {
	static_assert(BlockBits >= 128 && BlockBits <= 65536 && std::has_single_bit(unsigned(BlockBits)));
	static_assert(Branch >= 4 && std::has_single_bit(unsigned(Branch)));
	using Word = std::uint64_t;
	static constexpr int words = BlockBits / 64;
	struct Leaf {
		std::array<Word, words> data{};
		std::array<std::uint16_t, words> prefix{};
		int length = 0, ones = 0;
	};
	struct Node {
		std::array<int, Branch> ids{};
		std::array<int, Branch + 1> sizes{}, sums{};
		int count = 0;
	};
	std::vector<Leaf> leaves;
	std::vector<Node> nodes;
	int root = 0, n = 0, total = 0, free_leaf = 0, free_node = 0;
	static Word mask(int length) { return length == 64 ? ~Word(0) : (Word(1) << length) - 1; }
	static Word extract(const Word *data, int count, int p) {
		int w = p >> 6, shift = p & 63;
		Word value = w < count ? data[w] >> shift : 0;
		if (shift && w + 1 < count) value |= data[w + 1] << (64 - shift);
		return value;
	}
	int length(int id) const { return id & 1 ? leaves[id >> 1].length : nodes[id >> 1].sizes[nodes[id >> 1].count]; }
	int ones(int id) const { return id & 1 ? leaves[id >> 1].ones : nodes[id >> 1].sums[nodes[id >> 1].count]; }
	int make_leaf(const Leaf &leaf) {
		if (free_leaf) {
			int id = free_leaf;
			free_leaf = leaves[id >> 1].length;
			leaves[id >> 1] = leaf;
			return id;
		}
		assert(leaves.size() < (std::size_t(1) << 30));
		leaves.push_back(leaf);
		return (int(leaves.size()) - 1) * 2 + 1;
	}
	void release_leaf(int id) {
		leaves[id >> 1].length = free_leaf;
		free_leaf = id;
	}
	void pull(int id) {
		auto &node = nodes[id >> 1];
		int size = 0, count = 0;
		for (int i = 0; i < node.count; ++i) {
			node.sizes[i + 1] = size += length(node.ids[i]);
			node.sums[i + 1] = count += ones(node.ids[i]);
		}
	}
	int make_node(const int *children, int count) {
		Node node;
		node.count = count;
		int size = 0, sum = 0;
		for (int i = 0; i < count; ++i) node.ids[i] = children[i], node.sizes[i + 1] = size += length(children[i]), node.sums[i + 1] = sum += ones(children[i]);
		if (free_node) {
			int id = free_node;
			free_node = nodes[id >> 1].count;
			nodes[id >> 1] = node;
			return id;
		}
		if (nodes.empty()) nodes.emplace_back();
		assert(nodes.size() < (std::size_t(1) << 30));
		nodes.push_back(node);
		return (int(nodes.size()) - 1) * 2;
	}
	void release_node(int id) {
		nodes[id >> 1].count = free_node;
		free_node = id;
	}
	static int child(const Node &node, int p) {
		if constexpr (Branch <= 8) {
			int i = 0;
			while (i + 1 < node.count && p >= node.sizes[i + 1]) ++i;
			return i;
		} else {
			int l = 0, r = node.count - 1;
			while (l < r) {
				int m = (l + r) / 2;
				if (p < node.sizes[m + 1]) r = m;
				else l = m + 1;
			}
			return l;
		}
	}
	static int rank_leaf(const Leaf &leaf, int p) {
		if (p == leaf.length) return leaf.ones;
		int w = p >> 6;
		return leaf.prefix[w] + std::popcount(leaf.data[w] & mask(p & 63));
	}
	int rank_node(int id, int p) const {
		int result = 0;
		while (!(id & 1)) {
			const auto &node = nodes[id >> 1];
			int j = child(node, p);
			p -= node.sizes[j], result += node.sums[j];
			id = node.ids[j];
		}
		return result + rank_leaf(leaves[id >> 1], p);
	}
	static void insert_leaf(Leaf &leaf, int p, bool value) {
		int w = p >> 6, shift = p & 63;
		if (!(leaf.length & 63)) leaf.prefix[leaf.length >> 6] = leaf.ones;
		for (int i = leaf.length >> 6; i > w; --i) {
			int carry = leaf.data[i - 1] >> 63;
			leaf.data[i] = (leaf.data[i] << 1) | carry;
			leaf.prefix[i] += int(value) - carry;
		}
		Word low = mask(shift), x = leaf.data[w];
		leaf.data[w] = (x & low) | ((x & ~low) << 1) | (Word(value) << shift);
		++leaf.length;
		leaf.ones += value;
	}
	static void erase_leaf(Leaf &leaf, int p) {
		int w = p >> 6, last = (leaf.length - 1) >> 6;
		Word low = mask(p & 63), x = leaf.data[w];
		int value = (x >> (p & 63)) & 1;
		leaf.ones -= value;
		leaf.data[w] = (x & low) | ((x >> 1) & ~low);
		if (w < last) leaf.data[w] |= leaf.data[w + 1] << 63;
		for (int i = w + 1; i <= last; ++i) {
			leaf.prefix[i] += int(leaf.data[i] & 1) - value;
			leaf.data[i] = (leaf.data[i] >> 1) | (i < last ? leaf.data[i + 1] << 63 : 0);
		}
		--leaf.length;
	}
	static void assign_leaf(Leaf &leaf, const Word *data, int count, int start, int length) {
		leaf = Leaf{};
		leaf.length = length;
		for (int i = 0; i < (length + 63) / 64; ++i) {
			Word value = extract(data, count, start + i * 64);
			if (length - i * 64 < 64) value &= mask(length - i * 64);
			leaf.data[i] = value;
			leaf.prefix[i] = leaf.ones;
			leaf.ones += std::popcount(value);
		}
	}
	bool rebalance(int parent, int j) {
		auto &node = nodes[parent >> 1];
		int l = j ? j - 1 : 0, r = l + 1;
		int left = node.ids[l], right = node.ids[r];
		bool merge;
		if (left & 1) {
			auto &a = leaves[left >> 1], &b = leaves[right >> 1];
			std::array<Word, 2 * words> packed{};
			std::copy(a.data.begin(), a.data.end(), packed.begin());
			int offset = a.length >> 6, shift = a.length & 63, size = a.length + b.length;
			for (int i = 0; i < (b.length + 63) / 64; ++i) {
				packed[offset + i] |= b.data[i] << shift;
				if (shift) packed[offset + i + 1] |= b.data[i] >> (64 - shift);
			}
			merge = size <= BlockBits;
			int first = merge ? size : size / 2;
			assign_leaf(a, packed.data(), int(packed.size()), 0, first);
			if (merge) release_leaf(right);
			else assign_leaf(b, packed.data(), int(packed.size()), first, size - first);
		} else {
			auto &a = nodes[left >> 1], &b = nodes[right >> 1];
			std::array<int, 2 * Branch> children;
			int size = a.count + b.count;
			for (int i = 0; i < a.count; ++i) children[i] = a.ids[i];
			for (int i = 0; i < b.count; ++i) children[a.count + i] = b.ids[i];
			merge = size <= Branch;
			a.count = merge ? size : size / 2;
			for (int i = 0; i < a.count; ++i) a.ids[i] = children[i];
			pull(left);
			if (merge) release_node(right);
			else {
				b.count = size - a.count;
				for (int i = 0; i < b.count; ++i) b.ids[i] = children[a.count + i];
				pull(right);
			}
		}
		if (merge) {
			for (int i = r; i + 1 < node.count; ++i) node.ids[i] = node.ids[i + 1];
			--node.count;
		}
		pull(parent);
		return merge;
	}
	static int select_word(Word value, int k) {
		int result = 0;
		for (int shift : {32, 16, 8}) {
			int count = std::popcount(value & mask(shift));
			if (k >= count) k -= count, value >>= shift, result += shift;
			else value &= mask(shift);
		}
		while (k--) value &= value - 1;
		return result + std::countr_zero(value);
	}

  public:
	DynamicBitVector() = default;
	explicit DynamicBitVector(const std::vector<bool> &values) { build(values); }
	DynamicBitVector(const DynamicBitVector &) = default;
	DynamicBitVector &operator=(const DynamicBitVector &) = default;
	DynamicBitVector(DynamicBitVector &&other) noexcept { swap(other); }
	DynamicBitVector &operator=(DynamicBitVector &&other) noexcept {
		if (this != &other) clear(), swap(other);
		return *this;
	}
	void swap(DynamicBitVector &other) noexcept {
		leaves.swap(other.leaves);
		nodes.swap(other.nodes);
		std::swap(root, other.root), std::swap(n, other.n), std::swap(total, other.total);
		std::swap(free_leaf, other.free_leaf), std::swap(free_node, other.free_node);
	}
	int size() const { return n; }
	bool empty() const { return !n; }
	int count() const { return total; }
	std::size_t memory_usage() const { return leaves.capacity() * sizeof(Leaf) + nodes.capacity() * sizeof(Node); }
	void clear() {
		leaves.clear(), nodes.clear();
		root = n = total = free_leaf = free_node = 0;
	}
	void reserve(int size) {
		assert(size >= 0);
		if (!size) return;
		std::size_t capacity = std::size_t(size) / (BlockBits / 2) + 2;
		leaves.reserve(capacity);
		if (size > BlockBits) nodes.reserve(capacity / (Branch / 2 - 1) + 4);
	}
	void build(int size, const std::vector<Word> &packed) {
		assert(size >= 0 && packed.size() >= (std::size_t(size) + 63) / 64);
		clear();
		if (!size) return;
		n = size;
		int target = BlockBits * 3 / 4, count = size <= BlockBits ? 1 : (size - 1) / target + 1;
		leaves.reserve(count);
		std::vector<int> current;
		current.reserve(count);
		for (int i = 0, start = 0; i < count; ++i) {
			int length = size / count + (i < size % count);
			Leaf leaf;
			assign_leaf(leaf, packed.data(), int(packed.size()), start, length);
			total += leaf.ones;
			current.push_back(make_leaf(leaf));
			start += length;
		}
		while (current.size() > 1) {
			int count = (int(current.size()) - 1) / Branch + 1;
			std::vector<int> next;
			next.reserve(count);
			for (int i = 0, start = 0; i < count; ++i) {
				int length = int(current.size()) / count + (i < int(current.size()) % count);
				next.push_back(make_node(current.data() + start, length));
				start += length;
			}
			current.swap(next);
		}
		root = current[0];
	}
	void build(const std::vector<bool> &values) {
		assert(values.size() <= std::size_t(std::numeric_limits<int>::max()));
		std::vector<Word> packed((values.size() + 63) / 64);
		for (std::size_t i = 0; i < values.size(); ++i) packed[i >> 6] |= Word(values[i]) << (i & 63);
		build(int(values.size()), packed);
	}
	int rank1(int r) const {
		assert(0 <= r && r <= n);
		if (!total) return 0;
		if (total == n) return r;
		if (!r) return 0;
		if (r == n) return total;
		return rank_node(root, r);
	}
	int rank0(int r) const { return r - rank1(r); }
	std::pair<int, int> rank1(int l, int r) const {
		assert(0 <= l && l <= r && r <= n);
		if (!total) return {0, 0};
		if (total == n) return {l, r};
		if (!l) return {0, rank1(r)};
		if (r == n) return {rank1(l), total};
		int id = root, before = 0;
		while (!(id & 1)) {
			const auto &node = nodes[id >> 1];
			int a = child(node, l), b;
			if constexpr (Branch <= 8) {
				b = a;
				while (b + 1 < node.count && r >= node.sizes[b + 1]) ++b;
			} else b = child(node, r);
			if (a != b) {
				int x = node.sums[a], y = node.sums[b];
				return {before + x + rank_node(node.ids[a], l - node.sizes[a]),
						before + y + rank_node(node.ids[b], r - node.sizes[b])};
			}
			l -= node.sizes[a], r -= node.sizes[a], before += node.sums[a];
			id = node.ids[a];
		}
		const auto &leaf = leaves[id >> 1];
		return {before + rank_leaf(leaf, l), before + rank_leaf(leaf, r)};
	}
	std::pair<int, bool> access_rank(int p) const {
		assert(0 <= p && p < n);
		if (!total) return {0, false};
		if (total == n) return {p, true};
		int id = root, before = 0;
		while (!(id & 1)) {
			const auto &node = nodes[id >> 1];
			int j = child(node, p);
			p -= node.sizes[j], before += node.sums[j];
			id = node.ids[j];
		}
		const auto &leaf = leaves[id >> 1];
		return {before + rank_leaf(leaf, p), (leaf.data[p >> 6] >> (p & 63)) & 1};
	}
	bool access(int p) const {
		assert(0 <= p && p < n);
		if (!total || total == n) return total != 0;
		int id = root;
		while (!(id & 1)) {
			const auto &node = nodes[id >> 1];
			int j = child(node, p);
			p -= node.sizes[j];
			id = node.ids[j];
		}
		return (leaves[id >> 1].data[p >> 6] >> (p & 63)) & 1;
	}
	bool operator[](int p) const { return access(p); }
	int insert(int p, bool value) {
		assert(0 <= p && p <= n && n < std::numeric_limits<int>::max());
		int known = !total ? 0 : total == n ? p
							 : p == n		? total
											: -1;
		if (!root) {
			Leaf leaf;
			insert_leaf(leaf, 0, value);
			root = make_leaf(leaf);
			++n, total += value;
			return 0;
		}
		std::array<int, 32> path, slots;
		int id = root, depth = 0, before = 0;
		while (!(id & 1)) {
			const auto &node = nodes[id >> 1];
			int j = child(node, p);
			path[depth] = id, slots[depth++] = j;
			p -= node.sizes[j], before += node.sums[j];
			id = node.ids[j];
		}
		before = known >= 0 ? known : before + rank_leaf(leaves[id >> 1], p);
		int split = 0;
		if (leaves[id >> 1].length == BlockBits) {
			Leaf right;
			auto &left = leaves[id >> 1];
			right.length = BlockBits / 2;
			left.length = BlockBits / 2;
			std::copy(left.data.begin() + words / 2, left.data.end(), right.data.begin());
			std::fill(left.data.begin() + words / 2, left.data.end(), 0);
			right.ones = left.ones - left.prefix[words / 2];
			left.ones -= right.ones;
			for (int i = 0; i < words / 2; ++i) right.prefix[i] = left.prefix[i + words / 2] - left.ones;
			split = make_leaf(right);
			if (p <= BlockBits / 2) insert_leaf(leaves[id >> 1], p, value);
			else insert_leaf(leaves[split >> 1], p - BlockBits / 2, value);
		} else insert_leaf(leaves[id >> 1], p, value);
		while (depth) {
			id = path[--depth];
			int j = slots[depth];
			auto &node = nodes[id >> 1];
			if (!split) {
				for (int i = j; i < node.count; ++i) ++node.sizes[i + 1], node.sums[i + 1] += value;
				continue;
			}
			std::array<int, Branch + 1> children;
			for (int i = 0; i <= j; ++i) children[i] = node.ids[i];
			children[j + 1] = split;
			for (int i = j + 1; i < node.count; ++i) children[i + 1] = node.ids[i];
			int count = node.count + 1;
			node.count = count <= Branch ? count : count / 2;
			for (int i = 0; i < node.count; ++i) node.ids[i] = children[i];
			pull(id);
			int first = node.count;
			split = count <= Branch ? 0 : make_node(children.data() + first, count - first);
		}
		if (split) {
			std::array<int, 2> children{root, split};
			root = make_node(children.data(), 2);
		}
		++n, total += value;
		return before;
	}
	std::pair<int, bool> erase_rank(int p) {
		assert(0 <= p && p < n);
		int known = !total ? 0 : total == n ? p
											: -1;
		std::array<int, 32> path, slots;
		int id = root, depth = 0, before = 0;
		while (!(id & 1)) {
			const auto &node = nodes[id >> 1];
			int j = child(node, p);
			path[depth] = id, slots[depth++] = j;
			p -= node.sizes[j], before += node.sums[j];
			id = node.ids[j];
		}
		auto &leaf = leaves[id >> 1];
		before = known >= 0 ? known : before + rank_leaf(leaf, p);
		bool value = (leaf.data[p >> 6] >> (p & 63)) & 1;
		erase_leaf(leaf, p);
		--n, total -= value;
		if (!n) {
			release_leaf(root);
			root = 0;
			return {before, value};
		}
		for (int d = 0; d < depth; ++d) {
			auto &node = nodes[path[d] >> 1];
			for (int i = slots[d]; i < node.count; ++i) --node.sizes[i + 1], node.sums[i + 1] -= value;
		}
		while (depth && (id & 1 ? leaves[id >> 1].length < BlockBits / 2 : nodes[id >> 1].count < Branch / 2)) {
			int parent = path[--depth];
			if (!rebalance(parent, slots[depth])) break;
			id = parent;
		}
		if (!(root & 1) && nodes[root >> 1].count == 1) {
			int old = root;
			root = nodes[old >> 1].ids[0];
			release_node(old);
		}
		return {before, value};
	}
	bool erase(int p) { return erase_rank(p).second; }
	void push_back(bool value) { insert(n, value); }
	bool pop_back() { return erase(n - 1); }
	std::pair<int, bool> set_rank(int p, bool value) {
		assert(0 <= p && p < n);
		if ((!total && !value) || (total == n && value)) return {value ? p : 0, value};
		std::array<int, 32> path, slots;
		int id = root, depth = 0, before = 0;
		while (!(id & 1)) {
			const auto &node = nodes[id >> 1];
			int j = child(node, p);
			path[depth] = id, slots[depth++] = j;
			p -= node.sizes[j], before += node.sums[j];
			id = node.ids[j];
		}
		auto &leaf = leaves[id >> 1];
		before += rank_leaf(leaf, p);
		bool old = (leaf.data[p >> 6] >> (p & 63)) & 1;
		if (old != value) {
			int delta = int(value) - int(old);
			leaf.data[p >> 6] ^= Word(1) << (p & 63);
			leaf.ones += delta, total += delta;
			for (int i = (p >> 6) + 1; i <= (leaf.length - 1) / 64; ++i) leaf.prefix[i] += delta;
			for (int d = 0; d < depth; ++d) {
				auto &node = nodes[path[d] >> 1];
				for (int i = slots[d]; i < node.count; ++i) node.sums[i + 1] += delta;
			}
		}
		return {before, old};
	}
	std::tuple<int, bool, int> move_rank(int from, int to, bool value) {
		assert(0 <= from && from < n && 0 <= to && to < n);
		if (from == to) {
			auto [rank, old] = set_rank(from, value);
			return {rank, old, rank};
		}
		std::array<int, 32> left_path, right_path, left_slots, right_slots;
		int a = root, b = root, p = from, q = to + (from < to), x = 0, y = 0, left_depth = 0, right_depth = 0;
		while (!(a & 1)) {
			const auto &node = nodes[a >> 1];
			int j = child(node, p);
			left_path[left_depth] = a, left_slots[left_depth++] = j;
			p -= node.sizes[j], x += node.sums[j];
			a = node.ids[j];
		}
		while (!(b & 1)) {
			const auto &node = nodes[b >> 1];
			int j = child(node, q);
			right_path[right_depth] = b, right_slots[right_depth++] = j;
			q -= node.sizes[j], y += node.sums[j];
			b = node.ids[j];
		}
		auto &left = leaves[a >> 1], &right = leaves[b >> 1];
		if (a != b && (left.length <= BlockBits / 2 || right.length == BlockBits)) {
			auto [rank, old] = erase_rank(from);
			int next = insert(to, value);
			return {rank, old, next};
		}
		x += rank_leaf(left, p), y += rank_leaf(right, q);
		bool old = (left.data[p >> 6] >> (p & 63)) & 1;
		y -= (from < to) && old;
		erase_leaf(left, p);
		insert_leaf(right, q - (a == b && from < to), value);
		total += int(value) - int(old);
		int common = 0;
		while (common < left_depth && common < right_depth && left_path[common] == right_path[common]) {
			auto &node = nodes[left_path[common] >> 1];
			int l = left_slots[common], r = right_slots[common];
			int delta = int(value) - int(old);
			for (int i = std::min(l, r); i < std::max(l, r); ++i) {
				node.sizes[i + 1] += l < r ? -1 : 1;
				node.sums[i + 1] += l < r ? -int(old) : int(value);
			}
			if (delta)
				for (int i = std::max(l, r); i < node.count; ++i) node.sums[i + 1] += delta;
			++common;
		}
		for (int d = common; d < left_depth; ++d) {
			auto &node = nodes[left_path[d] >> 1];
			for (int i = left_slots[d]; i < node.count; ++i) --node.sizes[i + 1], node.sums[i + 1] -= old;
		}
		for (int d = common; d < right_depth; ++d) {
			auto &node = nodes[right_path[d] >> 1];
			for (int i = right_slots[d]; i < node.count; ++i) ++node.sizes[i + 1], node.sums[i + 1] += value;
		}
		return {x, old, y};
	}

	void set(int p, bool value) { set_rank(p, value); }
	int select(bool value, int k) const {
		if (k < 0 || k >= (value ? total : n - total)) return -1;
		if (!total || total == n) return k;
		int id = root, before = 0;
		while (!(id & 1)) {
			const auto &node = nodes[id >> 1];
			int l = 0, r = node.count - 1;
			while (l < r) {
				int m = (l + r) / 2;
				int count = value ? node.sums[m + 1] : node.sizes[m + 1] - node.sums[m + 1];
				if (k < count) r = m;
				else l = m + 1;
			}
			if (l) {
				before += node.sizes[l];
				k -= value ? node.sums[l] : node.sizes[l] - node.sums[l];
			}
			id = node.ids[l];
		}
		const auto &leaf = leaves[id >> 1];
		for (int i = 0; i < (leaf.length + 63) / 64; ++i) {
			Word word = value ? leaf.data[i] : ~leaf.data[i];
			if (leaf.length - i * 64 < 64) word &= mask(leaf.length - i * 64);
			int count = std::popcount(word);
			if (k < count) return before + i * 64 + select_word(word, k);
			k -= count;
		}
		return -1;
	}
	int select0(int k) const { return select(false, k); }
	int select1(int k) const { return select(true, k); }
};

#endif
