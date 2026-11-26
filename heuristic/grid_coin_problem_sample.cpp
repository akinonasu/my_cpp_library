#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

#include "heuristic.hpp"

constexpr int H = 300;
constexpr int W = 300;
constexpr long TURNS = 400;
constexpr long BEAM_WIDTH = 300;
constexpr int DR[4] = {-1, 1, 0, 0};  // 0:上, 1:下, 2:左, 3:右 (dir ^ 1 が逆方向)
constexpr int DC[4] = {0, 0, -1, 1};

struct GridProblem {
    std::array<long, H * W> coin{};
    std::array<std::uint64_t, H * W> z_coin{};
    std::array<std::uint64_t, H * W> z_pos{};

    GridProblem() {
        Random rng(999);
        for (int i = 0; i < H * W; ++i) {
            coin[i] = rng.next_long(1, 101);
            z_coin[i] = rng.next_u64();
            z_pos[i] = rng.next_u64();
        }
    }
};

// パターン A: 標準ビームサーチ (2段階展開・コピー型)
void solve_with_beam_search(const GridProblem& prob) {
    struct State {
        int r = 0;
        int c = 0;
        long score = 0;
        std::uint64_t hash = 0;
        std::array<std::uint8_t, H * W> visited{};
    };

    Timer timer;
    State init;
    init.visited[0] = 1;
    init.score = prob.coin[0];
    init.hash = prob.z_pos[0] ^ prob.z_coin[0];

    BeamSearch<State, int, long> bs(init, init.score);
    bs.reserve_cands(BEAM_WIDTH * 4);

    for (long t = 0; t < TURNS; ++t) {
        for (long i = 0; i < bs.size(); ++i) {
            const State& s = bs.state(i);
            for (int dir = 0; dir < 4; ++dir) {
                int nr = s.r + DR[dir];
                int nc = s.c + DC[dir];
                if (nr < 0 || nr >= H || nc < 0 || nc >= W) {
                    continue;
                }
                int nidx = nr * W + nc;
                int oidx = s.r * W + s.c;
                long nscore = s.score + (s.visited[nidx] ? 0 : prob.coin[nidx]);
                std::uint64_t nhash = s.hash ^ prob.z_pos[oidx] ^ prob.z_pos[nidx] ^ (s.visited[nidx] ? 0ULL : prob.z_coin[nidx]);
                bs.push(nscore, i, dir, nhash);
            }
        }
        // 生き残った上位 BEAM_WIDTH 個にのみ状態を複製して遷移を適用
        bs.advance(BEAM_WIDTH, [&](State& s, int dir) {
            int oidx = s.r * W + s.c;
            s.r += DR[dir];
            s.c += DC[dir];
            int nidx = s.r * W + s.c;
            s.hash ^= prob.z_pos[oidx] ^ prob.z_pos[nidx];
            if (!s.visited[nidx]) {
                s.visited[nidx] = 1;
                s.score += prob.coin[nidx];
                s.hash ^= prob.z_coin[nidx];
            }
        });
    }

    bs.dump();
    std::vector<int> actions = bs.restore_actions();
    std::cout << "[BeamSearch]      score=" << bs.best_score() << ", steps=" << std::ssize(actions) << ", time=" << timer.elapsed_us()
              << " us\n";
}

// パターン B: ロールバック型ビームサーチ (状態を 1 つだけ持ち、コピーしない)
void solve_with_euler_beam_search(const GridProblem& prob) {
    struct SingleState {
        int r = 0;
        int c = 0;
        long score = 0;
        std::uint64_t hash = 0;
        std::array<int, H * W> visit_cnt{};
    };

    Timer timer;
    SingleState state;
    state.visit_cnt[0] = 1;
    state.score = prob.coin[0];
    state.hash = prob.z_pos[0] ^ prob.z_coin[0];

    auto apply = [&](int dir) {
        int oidx = state.r * W + state.c;
        state.r += DR[dir];
        state.c += DC[dir];
        int nidx = state.r * W + state.c;
        state.hash ^= prob.z_pos[oidx] ^ prob.z_pos[nidx];
        if (state.visit_cnt[nidx]++ == 0) {
            state.score += prob.coin[nidx];
            state.hash ^= prob.z_coin[nidx];
        }
    };

    auto undo = [&](int dir) {
        int oidx = state.r * W + state.c;
        if (--state.visit_cnt[oidx] == 0) {
            state.score -= prob.coin[oidx];
            state.hash ^= prob.z_coin[oidx];
        }
        int rev = dir ^ 1;
        state.r += DR[rev];
        state.c += DC[rev];
        int nidx = state.r * W + state.c;
        state.hash ^= prob.z_pos[oidx] ^ prob.z_pos[nidx];
    };

    EulerBeamSearch<int, long> bs;
    bs.reserve_cands(BEAM_WIDTH * 4);

    bs.run(TURNS, BEAM_WIDTH, apply, undo, [&]() {
        int oidx = state.r * W + state.c;
        for (int dir = 0; dir < 4; ++dir) {
            int nr = state.r + DR[dir];
            int nc = state.c + DC[dir];
            if (nr < 0 || nr >= H || nc < 0 || nc >= W) {
                continue;
            }
            int nidx = nr * W + nc;
            long nscore = state.score + (state.visit_cnt[nidx] == 0 ? prob.coin[nidx] : 0);
            std::uint64_t nhash =
                state.hash ^ prob.z_pos[oidx] ^ prob.z_pos[nidx] ^ (state.visit_cnt[nidx] == 0 ? prob.z_coin[nidx] : 0ULL);
            bs.push(nscore, dir, nhash);
        }
    });

    bs.dump();
    std::vector<int> actions = bs.best_actions();
    std::cout << "[EulerBeamSearch] score=" << bs.best_score() << ", steps=" << std::ssize(actions) << ", time=" << timer.elapsed_us()
              << " us\n";
}

int main() {
    std::cout << "----\n";
    GridProblem prob;
    solve_with_beam_search(prob);
    solve_with_euler_beam_search(prob);
}