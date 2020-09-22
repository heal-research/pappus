#ifndef PAPPUS_FPUTIL_HPP
#define PAPPUS_FPUTIL_HPP

#include <cfenv>
#include <cmath>
#include <iostream>
#include <sstream>

#define EXPECT(cond)                                                                                    \
    if (!(cond)) {                                                                                      \
        std::cerr << "precondition " << #cond << " failed at " << __FILE__ << ": " << __LINE__ << "\n"; \
        std::terminate();                                                                               \
    }

#define ENSURE(cond)                                                                                     \
    if (!(cond)) {                                                                                       \
        std::cerr << "postcondition " << #cond << " failed at " << __FILE__ << ": " << __LINE__ << "\n"; \
        std::terminate();                                                                                \
    }

namespace pappus {
namespace fp {
    const auto pi      = std::acos(-1);
    const auto half_pi = pi / 2;
    const auto two_pi  = 2 * pi;
    const auto tau     = two_pi;

    constexpr auto inf     = std::numeric_limits<double>::infinity();
    constexpr auto nan     = std::numeric_limits<double>::quiet_NaN();

    template <int ROUND_MODE>
    double from_string(std::string const& s)
    {
#if defined(DIRECTED_ROUNDING)
        static_assert(ROUND_MODE == FE_UPWARD || ROUND_MODE == FE_DOWNWARD);
        auto rnd = std::fegetround();
        std::fesetround(ROUND_MODE);
#endif
        std::istringstream is(s);
        double v;
        is.precision(std::numeric_limits<double>::max_digits10);
        is >> v;
#if defined(DIRECTED_ROUNDING)
        std::fesetround(rnd);
#endif
        return v;
    }
}
}

#endif
