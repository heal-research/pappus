#ifndef PAPPUS_INTERVAL_BOX_HPP
#define PAPPUS_INTERVAL_BOX_HPP

#include <algorithm>
#include <concepts>
#include <limits>
#include <numeric>
#include <queue>
#include <vector>

#include "interval.hpp"

namespace pappus {

template<std::floating_point T = double>
using box = std::vector<interval<T>>;

template<std::floating_point T = double, typename F>
interval<T> optimize_bounds(F&& f, box<T> const& x, bool m = false, T w = T(1e-5), size_t n = 1000)
{
    static_assert(std::is_invocable_r_v<interval<T>, F, box<T> const&>);

    using Tuple = std::tuple<T, T, box<T>>;

    auto get_bound  = [&](interval<T> iv) { return m ? -iv.sup() : iv.inf(); };
    auto get_volume = [](auto const& bx) {
        return std::transform_reduce(cbegin(bx), cend(bx), T(1), std::multiplies<T>{},
                                      [](auto iv) { return iv.diameter(); });
    };
    auto widest_splittable_dimension = [&](auto const& bx) -> std::optional<std::size_t> {
        std::optional<std::size_t> best;
        T best_diameter = w;
        for (std::size_t i = 0; i < bx.size(); ++i) {
            auto diameter = bx[i].diameter();
            if (diameter > best_diameter) {
                best = i;
                best_diameter = diameter;
            }
        }
        return best;
    };
    auto compare = [](auto const& a, auto const& b) {
        return std::tie(std::get<0>(a), std::get<1>(a)) > std::tie(std::get<0>(b), std::get<1>(b));
    };

    std::priority_queue<Tuple, std::vector<Tuple>, decltype(compare)> q(compare);
    q.push({ get_bound(f(x)), get_volume(x), x });

    T best_bound = fp::inf_v<T>;

    while (!q.empty() && n-- > 0) {
        auto [bound, volume, bx] = q.top(); q.pop();
        auto split_dimension = widest_splittable_dimension(bx);
        if (!split_dimension) {
            best_bound = std::min(best_bound, bound);
            continue;
        }

        auto [left, right] = bx[*split_dimension].split();
        auto push_child = [&](interval<T> child) {
            auto next = bx;
            next[*split_dimension] = child;
            q.push({ get_bound(f(next)), get_volume(next), std::move(next) });
        };

        push_child(left);
        push_child(right);
    }
    if (!q.empty())
        best_bound = std::min(best_bound, std::get<0>(q.top()));

    return m ? interval<T>(-fp::inf_v<T>, -best_bound)
             : interval<T>(best_bound, fp::inf_v<T>);
}

template<std::floating_point T = double, typename F>
interval<T> optimize_bounds(F&& f, box<T> const& x, T w = T(1e-5), size_t n = 1000)
{
    return optimize_bounds<T>(f, x, false, w, n) & optimize_bounds<T>(f, x, true, w, n);
}

} // namespace pappus
#endif
