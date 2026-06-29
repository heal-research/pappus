#ifndef PAPPUS_FPUTIL_HPP
#define PAPPUS_FPUTIL_HPP

#include <cfenv>
#include <concepts>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace pappus::fp {

// double constants kept for backward compat
inline constexpr double pi      = 3.141592653589793238462643383279502884L;
inline constexpr double half_pi = pi / 2.0;
inline constexpr double two_pi  = 2.0 * pi;
inline constexpr double tau     = two_pi;
inline constexpr double inf     = std::numeric_limits<double>::infinity();
inline constexpr double nan     = std::numeric_limits<double>::quiet_NaN();

// type-parametric versions for use in generic code
template<std::floating_point T> inline constexpr T pi_v      = T(3.141592653589793238462643383279502884L);
template<std::floating_point T> inline constexpr T half_pi_v = pi_v<T> / T(2);
template<std::floating_point T> inline constexpr T two_pi_v  = T(2) * pi_v<T>;
template<std::floating_point T> inline constexpr T tau_v     = two_pi_v<T>;
template<std::floating_point T> inline constexpr T inf_v     = std::numeric_limits<T>::infinity();
template<std::floating_point T> inline constexpr T nan_v     = std::numeric_limits<T>::quiet_NaN();

template<int ROUND_MODE, std::floating_point T = double>
T from_string(std::string const& s)
{
    // String-parsed interval construction is a cold path (run once at setup),
    // never inside an evaluation loop, so the fesetround thread-safety concern
    // that motivates the eve-based ropd/ropu path in math.hpp doesn't apply.
    // istringstream uses strtod which honours fesetround, giving 1-ULP outward
    // expansion (matching the rounding direction requested) instead of the
    // 2-ULP-wide interval that parse-then-nextafter would produce.
    static_assert(ROUND_MODE == FE_UPWARD || ROUND_MODE == FE_DOWNWARD);
    auto rnd = std::fegetround();
    std::fesetround(ROUND_MODE);
    std::istringstream is(s);
    T v = std::numeric_limits<T>::quiet_NaN();
    is.precision(std::numeric_limits<T>::max_digits10);
    is >> v;
    std::fesetround(rnd);
    if (!is || !is.eof())
        throw std::invalid_argument("pappus::fp::from_string: invalid floating-point literal");
    return v;
}

} // namespace pappus::fp

#endif
