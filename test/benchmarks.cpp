#include <catch2/catch_all.hpp>
#include <functional>
#include <random>
#include <sstream>
#include <nanobench.h>
#include <queue>
#include <algorithm>

#include "interval/interval.hpp"
#include "interval/box.hpp"
#include "pappus.hpp"

namespace fp = pappus::fp;

const auto inf = fp::inf;
const auto max = std::numeric_limits<double>::max();
const auto min = std::numeric_limits<double>::min();

using I = pappus::interval<double>;
const auto F = I::infinite();
const auto Z = I::zero();
const auto E = I::empty();

// ---------------------------------------------------------------------------
// interval splitting / optimizer benchmark
// ---------------------------------------------------------------------------

TEST_CASE("splitting", "[.]")
{
    I a(-1, 1), b(-1, 1), c(-1, 1);
    using X = std::vector<I>;
    auto f = [](auto&& x) { return x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * x[2] * (1 - x[2]); };

    X x { a, b, c };

    auto split_naive = [](auto&& x, auto&& f, int n = 1000) {
        int total = n * n;
        std::vector<int> counters(x.size(), 0ul);
        auto z = E;
        std::vector<I> box(x.size());
        while (total > 0) {
            for (auto i = counters.rbegin(); i < counters.rend(); ++i) {
                auto d = std::distance(i, counters.rend()) - 1;
                box[d] = x[d].segment(*i, n);
                z |= f(box);
                if (i > counters.rbegin()) {
                    auto j = i - 1;
                    if (*j == n) { *j = 0; (*i)++; }
                } else {
                    (*i)++;
                }
            }
            --total;
        }
        return z;
    };

    auto w = pappus::optimize_bounds(f, x, 1e-4, 1000);
    ankerl::nanobench::doNotOptimizeAway(w);

    SECTION("grid split performance")
    {
        ankerl::nanobench::Bench bench;
        bench.performanceCounters(true);
        I z;
        std::ostringstream ss;
        for (size_t i = 100; i <= 1000; i += 100) {
            z = split_naive(x, f, i);
            ss.str(""); ss.clear(); ss << z;
            bench.run("splitting " + std::to_string(i) + " = " + ss.str(), [&]() {
                ankerl::nanobench::doNotOptimizeAway(z = split_naive(x, f, static_cast<int>(i)));
            });
        }
    }
}

// ---------------------------------------------------------------------------
// IA vs AA comparison
// ---------------------------------------------------------------------------

TEST_CASE("affine vs interval arithmetic", "[.]")
{
    using af = pappus::affine_form<double>;

    ankerl::nanobench::Bench bench;
    bench.performanceCounters(true).warmup(5).epochIterations(50'000);

    // f(x) = x^2 - x on [-1, 1]
    // Classic dependency problem: IA naively overestimates, AA tracks correlations.
    // True range: [-0.25, 2]

    SECTION("IA: x^2 - x")
    {
        I x(-1, 1), result;
        bench.run("IA:         x^2 - x", [&]() {
            ankerl::nanobench::doNotOptimizeAway(result = x * x - x);
        });
    }

    SECTION("AA: x^2 - x")
    {
        I result;
        bench.run("AA:         x^2 - x", [&]() {
            pappus::affine_context c;
            af x(c, I(-1, 1));
            ankerl::nanobench::doNotOptimizeAway(result = (x * x - x).to_interval());
        });
    }

    for (int depth = 1; depth <= 3; ++depth) {
        SECTION("AA+bisect(" + std::to_string(depth) + "): x^2 - x")
        {
            I result;
            bench.run("AA+bisect(" + std::to_string(depth) + "): x^2 - x", [&]() {
                pappus::affine_context c;
                af x(c, I(-1, 1));
                auto f = [](af v) { return v * v - v; };
                ankerl::nanobench::doNotOptimizeAway(
                    result = pappus::evaluate_bisected(f, x, depth));
            });
        }
    }
}

// ---------------------------------------------------------------------------
// condense benchmark
// ---------------------------------------------------------------------------

TEST_CASE("condense", "[.]")
{
    using af = pappus::affine_form<double>;

    ankerl::nanobench::Bench bench;
    bench.performanceCounters(true).warmup(5).epochIterations(10'000);

    // Build a form with many terms by repeated multiplication
    auto make_wide = [](int n_ops) {
        pappus::affine_context c;
        af x(c, pappus::interval<double>(-1, 1));
        af y(c, pappus::interval<double>(0, 1));
        for (int i = 0; i < n_ops; ++i)
            x = x * y + x;
        return std::make_pair(std::move(c), std::move(x));
    };

    SECTION("AA with many terms, no condense")
    {
        I result;
        bench.run("AA wide (no condense)", [&]() {
            auto [c, x] = make_wide(8);
            ankerl::nanobench::doNotOptimizeAway(result = x.to_interval());
        });
    }

    SECTION("AA with many terms, condense to 8")
    {
        I result;
        bench.run("AA wide (condense 8)", [&]() {
            auto [c, x] = make_wide(8);
            ankerl::nanobench::doNotOptimizeAway(result = x.condense(8).to_interval());
        });
    }
}
