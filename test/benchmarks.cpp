#include "interval/interval.hpp"
#include <doctest/doctest.h>
#include <functional>
#include <iomanip>
#include <random>
#include <sstream>
#define ANKERL_NANOBENCH_IMPLEMENT
#include "nanobench.h"
#include <queue>
#include <algorithm>

#include "interval/box.hpp"

namespace dt = doctest;
namespace fp = pappus::fp;

const auto inf = fp::inf;
const auto max = std::numeric_limits<double>::max();
const auto min = std::numeric_limits<double>::min();

using I = pappus::interval;
const auto F = I::infinite();
const auto Z = I::zero();
const auto E = I::empty();

TEST_CASE("splitting")
{
    I a(-1, 1), b(-1, 1), c(-1, 1);
    using X = std::vector<I>;
    auto f = [](auto&& x) { return x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * x[2] * (1 - x[2]); };

    X x { a, b, c };
    auto y = f(x);
    std::cout << "IA y = " << y << "\n";

    auto split_naive = [](auto&& x, auto&& f, int n = 1000) {
        int total = n * n;
        std::vector<int> counters(x.size(), 0ul);

        // we are going to start with an empty interval
        auto z = E;

        // we are going to reuse a box (sequence of intervals)
        std::vector<I> box(x.size());

        size_t count = 0;
        // while our cartesian product is not done
        while (total > 0) {
            // update the counters for each sequence in the product
            for (auto i = counters.rbegin(); i < counters.rend(); ++i) {
                auto d = std::distance(i, counters.rend()) - 1;
                box[d] = x[d].segment(*i, n);

                z |= f(box);

                // check the counter to the right for wrap-around
                if (i > counters.rbegin()) {
                    auto j = i - 1;

                    // wrapped-around?
                    if (*j == n) {
                        *j = 0;
                        (*i)++;
                    }
                } else {
                    // this is the rightmost counter, just increment it
                    (*i)++;
                }
            }
            --total;
        }
        return z;
    };

    auto split_random = [](auto&& x, auto&& f, int n = 10, int s = 1000) {
        std::vector<I> box(x.size());
        std::vector<size_t> limits(x.size());
        std::vector<int> counters(x.size(), 0ul);
        std::default_random_engine rd(std::random_device {}());
        std::uniform_int_distribution<size_t> dist(0, n);

        auto y = F;
        while (s-- > 0) {
            std::generate(limits.begin(), limits.end(), [&]() { return std::pow(2ul, dist(rd)); });
            auto total = std::reduce(limits.begin(), limits.end(), 1ul, std::multiplies {});
            std::fill(counters.begin(), counters.end(), 0ul);

            auto z = E;
            while (total > 0) {
                for (auto i = counters.rbegin(); i < counters.rend(); ++i) {
                    auto d = std::distance(i, counters.rend()) - 1;
                    box[d] = x[d].segment(*i, limits[d]);

                    z |= f(box);

                    // check the counter to the right for wrap-around
                    if (i > counters.rbegin()) {
                        auto j = i - 1;

                        // wrapped-around?
                        if (*j == limits[std::distance(j, counters.rend()) - 1]) {
                            *j = 0;
                            (*i)++;
                        }
                    } else {
                        // this is the rightmost counter, just increment it
                        (*i)++;
                    }
                }
                --total;
            }
            y &= z;
        }
        return y;
    };

    auto w = pappus::optimize_bounds(f, x, 1e-4, 1000);
    std::cout << std::setprecision(10) << "w = " << w << "\n";

    //auto z = split_naive(x, f, 1000);
    //std::cout << "z = " << z << "\n";

    //ankerl::nanobench::Bench bench;
    //bench.performanceCounters(true);

    //SUBCASE("grid split performance")
    //{
    //    I z = split_naive(x, f, 10);
    //    std::ostringstream ss;
    //    ss << z;
    //    bench.run("splitting " + std::to_string(10) + " = " + ss.str(), [&]() {
    //        ankerl::nanobench::doNotOptimizeAway(z = split_naive(x, f, 10));
    //    });
    //    for (size_t i = 100; i <= 1000; i += 100) {
    //        z = split_naive(x, f, i);
    //        ss.str("");
    //        ss.clear();
    //        ss << z;
    //        bench.run("splitting " + std::to_string(i) + " = " + ss.str(), [&]() {
    //            ankerl::nanobench::doNotOptimizeAway(z = split_naive(x, f, i));
    //        });
    //    }
    //}

    //SUBCASE("random split performance")
    //{
    //    I w;
    //    auto s = 10;
    //    bench.run("sampling", [&]() {
    //        ankerl::nanobench::doNotOptimizeAway(w = split_random(x, f, s, 100));
    //    });
    //    std::cout << w << "\n";
    //};
}
