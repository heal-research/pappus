#ifndef PAPPUS_INTERVAL_BOX_HPP
#define PAPPUS_INTERVAL_BOX_HPP

#include <algorithm>
#include <numeric>
#include <vector>
#include <queue>

#include "interval.hpp"

namespace pappus {

using box = std::vector<interval>;

template<typename F>
interval optimize_bounds(F&& f, box const& x, bool m = false, double w = 1e-5, size_t n = 1000)
{
    static_assert(std::is_invocable_r_v<interval, F, box const&>);

    using T = std::tuple<double, double, box>;

    auto get_bound = [&](interval iv) { return m ? -iv.sup() : iv.inf(); };

    auto get_volume = [](auto&& box) {
        return std::transform_reduce(cbegin(box), cend(box), -1, std::multiplies{}, [](auto iv) { return iv.diameter(); });
    };

    auto compare = [](auto&& a, auto&& b) {
        return std::tie(std::get<0>(a), std::get<1>(a)) > std::tie(std::get<0>(b), std::get<1>(b));
    };

    // a priority queue to keep track of the bounds
    std::priority_queue<T, std::vector<T>, decltype(compare)> q(compare); 

    //auto best_bound = get_bound(f(x)); // best bound so far
    q.push({ get_bound(f(x)), get_volume(x), x });

    auto best_bound = fp::inf;
    box new_box(x.size());
    std::vector<box> splits(x.size());

    while(!q.empty() && n-- > 0) {
        auto [bound, volume, box] = q.top(); q.pop();
        // can we split the current hyperbox into smaller regions
        if (std::none_of(cbegin(box), cend(box), [&](auto iv) { return iv.diameter() > w; })) {
            best_bound = std::min(best_bound, bound);
            continue;
        }

        for (size_t i = 0; i < box.size(); ++i) {
            auto iv = box[i];
            // do the splitting and put the results in the splits vector
            splits[i].clear();
            if (iv.diameter() > w) {
                auto [left, right] = iv.split();
                splits[i].push_back(left);
                splits[i].push_back(right);
            } else {
                splits[i].push_back(iv);
            }
        }

        auto best_ = std::numeric_limits<double>::max();
        // do a cartesian product in-place to avoid temporaries using this lambda
        auto perm = [&](auto i, auto&& perm) {
            if (i == x.size()) {
                auto tmp = new_box; // make a copy
                q.push({ get_bound(f(tmp)), get_volume(tmp), tmp });
                return;
            }
            for (auto iv : splits[i]) {
                new_box[i] = iv; 
                perm(i+1, perm);
            }
        };
        perm(0, perm);
    }
    if (!q.empty()) {
        best_bound = std::min(best_bound, std::get<0>(q.top()));
    }
    // we could return a double here (the best bound) but an interval is cheap enough 
    return m ? interval(-fp::inf, -best_bound) : interval(best_bound, fp::inf);
}

template<typename F>
interval optimize_bounds(F&& f, box const& x, double w = 1e-5, size_t n = 1000)
{
    return optimize_bounds(f, x, false, w, n) & optimize_bounds(f, x, true, w, n);
}

} // namespace
#endif

