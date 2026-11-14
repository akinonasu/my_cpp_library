#pragma once
#include <algorithm>
#include <cassert>
#include <concepts>
#include <functional>
#include <iterator>
#include <memory>
#include <set>
#include <type_traits>
#include <utility>

// 評価関数1: 最近傍までの距離（要素数1のときは 0）
template <typename T, typename Res = long>
struct NearestNeighborDist {
    Res operator()(const T* L, const T& cur, const T* R) const {
        if (!L && !R) return Res{0};
        if (!L) return static_cast<Res>(*R - cur);
        if (!R) return static_cast<Res>(cur - *L);
        return static_cast<Res>(std::min(cur - *L, *R - cur));
    }
};

// 評価関数2: 右隣との差の2乗（隣接区間の関数の総和の例）
template <typename T, typename Res = long>
struct AdjacentDiffSq {
    Res operator()(const T* /*L*/, const T& cur, const T* R) const {
        if (!R) return Res{0};
        Res d = static_cast<Res>(*R - cur);
        return d * d;
    }
};

template <typename T, typename SumT = long, typename Eval = NearestNeighborDist<T, SumT>, typename Compare = std::less<T>, bool AllowDuplicate = true, bool TrackMinMax = false>
class NeighborContributionSet {
   private:
    using SetType = std::multiset<T, Compare>;
    using Iter = typename SetType::const_iterator;

    struct Empty {};

    Compare comp_{};
    Eval eval_{};
    SetType st_{comp_};
    SumT total_{};
    [[no_unique_address]] std::conditional_t<TrackMinMax, std::multiset<SumT>, Empty> contribs_{};

    bool equiv(const T& a, const T& b) const { return !comp_(a, b) && !comp_(b, a); }

    SumT calc(Iter it) const {
        const T* L = (it != st_.begin()) ? std::addressof(*std::prev(it)) : nullptr;
        Iter nxt = std::next(it);
        const T* R = (nxt != st_.end()) ? std::addressof(*nxt) : nullptr;
        return eval_(L, *it, R);
    }

    void add_contrib(Iter it) {
        SumT v = calc(it);
        total_ += v;
        if constexpr (TrackMinMax) {
            contribs_.insert(v);
        }
    }

    void remove_contrib(Iter it) {
        SumT v = calc(it);
        total_ -= v;
        if constexpr (TrackMinMax) {
            auto cit = contribs_.find(v);
            assert(cit != contribs_.end());
            contribs_.erase(cit);
        }
    }

   public:
    NeighborContributionSet() = default;

    explicit NeighborContributionSet(Eval eval, Compare comp = Compare{}) : comp_(std::move(comp)), eval_(std::move(eval)), st_(comp_) {}

    // 要素を追加する（AllowDuplicate == false で重複時は false を返す）
    bool insert(const T& x) {
        Iter it_next = st_.upper_bound(x);
        bool has_prev = (it_next != st_.begin());
        Iter it_prev = has_prev ? std::prev(it_next) : st_.end();

        if constexpr (!AllowDuplicate) {
            if (has_prev && equiv(*it_prev, x)) {
                return false;
            }
        }

        bool has_next = (it_next != st_.end());

        if (has_prev) remove_contrib(it_prev);
        if (has_next) remove_contrib(it_next);

        // upper_bound の直前に挿入することで it_prev と it_next の間に入ることが保証される
        Iter it = st_.insert(it_next, x);

        if (has_prev) add_contrib(it_prev);
        add_contrib(it);
        if (has_next) add_contrib(it_next);

        return true;
    }

    // 要素を1つ削除する（存在しない場合は false を返す）
    bool erase(const T& x) {
        Iter it = st_.find(x);
        if (it == st_.end()) return false;

        bool has_prev = (it != st_.begin());
        Iter it_prev = has_prev ? std::prev(it) : st_.end();
        Iter it_next = std::next(it);
        bool has_next = (it_next != st_.end());

        if (has_prev) remove_contrib(it_prev);
        remove_contrib(it);
        if (has_next) remove_contrib(it_next);

        st_.erase(it);

        if (has_prev) add_contrib(it_prev);
        if (has_next) add_contrib(it_next);

        return true;
    }

    void clear() {
        st_.clear();
        total_ = SumT{};
        if constexpr (TrackMinMax) {
            contribs_.clear();
        }
    }

    int size() const { return static_cast<int>(st_.size()); }
    bool empty() const { return st_.empty(); }
    int count(const T& x) const { return static_cast<int>(st_.count(x)); }
    bool contains(const T& x) const { return st_.find(x) != st_.end(); }

    // 寄与の総和を返す
    SumT total() const { return total_; }

    // 各要素の寄与の最小値・最大値を返す（TrackMinMax == true のときのみ有効）
    std::pair<bool, SumT> min_contrib() const
        requires TrackMinMax
    {
        if (contribs_.empty()) return {false, SumT{}};
        return {true, *contribs_.begin()};
    }

    std::pair<bool, SumT> max_contrib() const
        requires TrackMinMax
    {
        if (contribs_.empty()) return {false, SumT{}};
        return {true, *contribs_.rbegin()};
    }

    // 現在の集合内で x に最も近い要素を 1 つ返す
    std::pair<bool, T> nearest_element(const T& x) const
        requires std::same_as<Compare, std::less<T>> && requires(T a, T b) { a - b; }
    {
        if (st_.empty()) return {false, T{}};
        Iter it = st_.lower_bound(x);
        if (it == st_.begin()) return {true, *it};
        if (it == st_.end()) return {true, *std::prev(it)};
        Iter pit = std::prev(it);
        if (x - *pit <= *it - x) return {true, *pit};
        return {true, *it};
    }

    const SetType& data() const { return st_; }
};
