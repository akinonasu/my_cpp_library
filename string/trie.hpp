#include <array>
#include <cassert>
#include <iostream>
#include <optional>
#include <print>
#include <string>
#include <string_view>
#include <vector>

/**
 * @brief 競技プログラミング用 Trie木 (完全版)
 * @tparam ALPHABET_SIZE 文字種数 (デフォルト: 英小文字 26)
 * @tparam BASE 基準となる最小文字 (デフォルト: 'a')
 */
template <int ALPHABET_SIZE = 26, char BASE = 'a'>
struct Trie {
    struct Node {
        std::array<int, ALPHABET_SIZE> next;
        int exist_count = 0;    // このノードで終了する単語数
        int subtree_count = 0;  // このノードを通過する単語の総数
        int word_id = -1;       // 単語に付与された固有ID（不要なら -1）

        Node() { next.fill(-1); }
    };

    std::vector<Node> nodes;

    Trie() {
        nodes.emplace_back();  // ルートノード (index: 0)
    }

    // ノード数の事前確保（メモリ再確保のオーバーヘッドを削減）
    void reserve(int max_nodes) { nodes.reserve(max_nodes); }

    // ツリーの初期化
    void clear() {
        nodes.clear();
        nodes.emplace_back();
    }

    // 単語の追加 (ID指定可能, count個追加) - O(|word|)
    void insert(std::string_view word, int id = -1, int count = 1) {
        int node_id = 0;
        nodes[node_id].subtree_count += count;

        for (char ch : word) {
            int c = ch - BASE;
            assert(0 <= c && c < ALPHABET_SIZE);
            if (nodes[node_id].next[c] == -1) {
                nodes[node_id].next[c] = static_cast<int>(nodes.size());
                nodes.emplace_back();
            }
            node_id = nodes[node_id].next[c];
            nodes[node_id].subtree_count += count;
        }
        nodes[node_id].exist_count += count;
        if (id != -1) {
            nodes[node_id].word_id = id;
        }
    }

    // 単語の削除 (存在する場合のみ count 個削除) - O(|word|)
    bool erase(std::string_view word, int count = 1) {
        if (this->count(word) < count) return false;

        int node_id = 0;
        nodes[node_id].subtree_count -= count;

        for (char ch : word) {
            int c = ch - BASE;
            node_id = nodes[node_id].next[c];
            nodes[node_id].subtree_count -= count;
        }
        nodes[node_id].exist_count -= count;
        if (nodes[node_id].exist_count == 0) {
            nodes[node_id].word_id = -1;
        }
        return true;
    }

    // 指定した接頭辞を持つノード番号を取得 (存在しない場合は -1) - O(|prefix|)
    int find(std::string_view prefix) const {
        int node_id = 0;
        for (char ch : prefix) {
            int c = ch - BASE;
            if (c < 0 || c >= ALPHABET_SIZE || nodes[node_id].next[c] == -1) {
                return -1;
            }
            node_id = nodes[node_id].next[c];
        }
        return node_id;
    }

    // 単語の完全一致数を取得 - O(|word|)
    int count(std::string_view word) const {
        int node_id = find(word);
        return (node_id == -1) ? 0 : nodes[node_id].exist_count;
    }

    // 単語が存在するか判定 - O(|word|)
    bool contains(std::string_view word) const { return count(word) > 0; }

    // 単語のIDを取得 (存在しない場合は -1) - O(|word|)
    int get_id(std::string_view word) const {
        int node_id = find(word);
        return (node_id == -1) ? -1 : nodes[node_id].word_id;
    }

    // 指定した接頭辞から始まる単語の総数を取得 - O(|prefix|)
    int count_prefix(std::string_view prefix) const {
        int node_id = find(prefix);
        return (node_id == -1) ? 0 : nodes[node_id].subtree_count;
    }

    // 辞書順で k 番目 (0-indexed) に小さい文字列を取得 - O(|word| * ALPHABET_SIZE)
    std::optional<std::string> kth_element(int k) const {
        if (k < 0 || k >= size()) return std::nullopt;

        int node_id = 0;
        std::string result = "";

        while (true) {
            if (k < nodes[node_id].exist_count) {
                return result;  // 現在のノードで終わる単語にヒット
            }
            k -= nodes[node_id].exist_count;

            for (int c = 0; c < ALPHABET_SIZE; ++c) {
                int next_id = nodes[node_id].next[c];
                if (next_id == -1) continue;

                if (k < nodes[next_id].subtree_count) {
                    result += static_cast<char>(BASE + c);
                    node_id = next_id;
                    break;
                }
                k -= nodes[next_id].subtree_count;
            }
        }
    }

    // --- 部分文字列検索用機能 (Suffix Trie) ---

    // 文字列の全接尾辞を登録する - O(|word|^2)
    void insert_suffixes(std::string_view word) {
        for (int i = 0; i < static_cast<int>(word.size()); ++i) {
            insert(word.substr(i));
        }
    }

    // 部分文字列として含まれているか判定 - O(|pattern|)
    bool contains_substring(std::string_view pattern) const { return find(pattern) != -1; }

    // --- 情報取得 ---

    // 挿入されている単語の総数
    int size() const { return nodes[0].subtree_count; }

    // 作成されたノード総数
    int node_count() const { return static_cast<int>(nodes.size()); }
};
