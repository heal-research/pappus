#ifndef PAPPUS_FPOPS_HPP
#define PAPPUS_FPOPS_HPP

// Directed-rounding primitives for sound interval arithmetic.
//
// Arithmetic ops use eve's `lower`/`upper` decorators: correctly-rounded
// per-instruction rounding via the hardware MXCSR/FPCR save+restore. No
// process-global state, safe to call from multiple threads concurrently.
//
// Transcendentals use eve's polynomial approximations (round-to-nearest)
// followed by a 1-ULP outward expansion via eve::prev / eve::next. This
// yields sound interval bounds with a small (~1-4 ULP) overestimation
// versus a correctly-rounded implementation.
//
// When DIRECTED_ROUNDING is not defined, ropd/ropu degenerate to plain
// round-to-nearest computation: fast but not sound. Useful for performance
// testing or when soundness is guaranteed by other means.

#include "util.hpp"

#include <cmath>
#include <concepts>
#include <limits>
#include <type_traits>

#include <eve/eve.hpp>
#include <eve/module/core.hpp>
#include <eve/module/math.hpp>

namespace pappus::fp {

namespace detail {

// 1-ULP outward expansion for sound bounds.
// Preserve infinities/NaNs verbatim. Some functions (sin/tan/asin/atan/sinh/
// tanh/log at specific exact points) should also preserve an exact zero result
// instead of inflating it to a symmetric subnormal interval around zero.
template<typename T>
T outward_lo(T x)
{
    if constexpr (eve::value<T>) {
        auto nan_mask = eve::is_nan(x);
        auto pos_inf  = x == eve::inf(eve::as<T>{});
        auto neg_inf  = x == -eve::inf(eve::as<T>{});
        auto y = eve::prev(x);
        y = eve::if_else(pos_inf, eve::valmax(eve::as<T>{}), y);
        y = eve::if_else(neg_inf, x, y);
        return eve::if_else(nan_mask, x, y);
    } else {
        if (std::isnan(x)) return x;
        if (x == std::numeric_limits<T>::infinity()) return std::numeric_limits<T>::max();
        if (x == -std::numeric_limits<T>::infinity()) return x;
        return std::nextafter(x, -std::numeric_limits<T>::infinity());
    }
}

template<typename T>
T outward_hi(T x)
{
    if constexpr (eve::value<T>) {
        auto nan_mask = eve::is_nan(x);
        auto pos_inf  = x == eve::inf(eve::as<T>{});
        auto neg_inf  = x == -eve::inf(eve::as<T>{});
        auto y = eve::next(x);
        y = eve::if_else(pos_inf, x, y);
        y = eve::if_else(neg_inf, -eve::valmax(eve::as<T>{}), y);
        return eve::if_else(nan_mask, x, y);
    } else {
        if (std::isnan(x)) return x;
        if (x == std::numeric_limits<T>::infinity()) return x;
        if (x == -std::numeric_limits<T>::infinity()) return -std::numeric_limits<T>::max();
        return std::nextafter(x, +std::numeric_limits<T>::infinity());
    }
}

template<typename T>
T outward_lo_keep_zero(T x)
{
    if constexpr (eve::value<T>) {
        auto special = eve::is_not_finite(x) || eve::is_eqz(x);
        return eve::if_else(special, x, eve::prev(x));
    } else {
        return (!std::isfinite(x) || x == T(0)) ? x : std::nextafter(x, -std::numeric_limits<T>::infinity());
    }
}

template<typename T>
T outward_hi_keep_zero(T x)
{
    if constexpr (eve::value<T>) {
        auto special = eve::is_not_finite(x) || eve::is_eqz(x);
        return eve::if_else(special, x, eve::next(x));
    } else {
        return (!std::isfinite(x) || x == T(0)) ? x : std::nextafter(x, +std::numeric_limits<T>::infinity());
    }
}

} // namespace detail

// Operation tag types — used only for dispatch in ropd/ropu below.
struct op_add  {};
struct op_sub  {};
struct op_mul  {};
struct op_div  {};
struct op_abs  {};
struct op_exp  {};
struct op_log  {};
struct op_sqrt {};
struct op_sin  {};
struct op_cos  {};
struct op_tan  {};
struct op_asin {};
struct op_acos {};
struct op_atan {};
struct op_sinh {};
struct op_cosh {};
struct op_tanh {};
struct op_fmod {};
struct op_cbrt {};
struct op_log1p {};
struct op_floor {};
struct op_ceil  {};

// ---------------------------------------------------------------------------
// ropd: round-op toward -inf. ropu: round-op toward +inf.
// ---------------------------------------------------------------------------
template<typename OP, typename T, typename... Args>
T ropd(T first, Args... rest)
{
#if defined(DIRECTED_ROUNDING)
    if constexpr (std::is_same_v<OP, op_add>)        return eve::add[eve::lower](first, rest...);
    else if constexpr (std::is_same_v<OP, op_sub>)   return eve::sub[eve::lower](first, rest...);
    else if constexpr (std::is_same_v<OP, op_mul>)   return eve::mul[eve::lower](first, rest...);
    else if constexpr (std::is_same_v<OP, op_div>)   return eve::div[eve::lower](first, rest...);
    else if constexpr (std::is_same_v<OP, op_sqrt>)  return eve::sqrt[eve::lower](first);
    else if constexpr (std::is_same_v<OP, op_abs>)   return eve::abs(first);
    // Transcendentals: eve approximation + 1-ULP outward expansion.
    else if constexpr (std::is_same_v<OP, op_exp>)   return detail::outward_lo(eve::exp(first));
    else if constexpr (std::is_same_v<OP, op_log>)   return detail::outward_lo_keep_zero(eve::log(first));
    else if constexpr (std::is_same_v<OP, op_sin>)   return detail::outward_lo_keep_zero(eve::sin(first));
    else if constexpr (std::is_same_v<OP, op_cos>)   return detail::outward_lo(eve::cos(first));
    else if constexpr (std::is_same_v<OP, op_tan>)   return detail::outward_lo_keep_zero(eve::tan(first));
    else if constexpr (std::is_same_v<OP, op_asin>)  return detail::outward_lo_keep_zero(eve::asin(first));
    else if constexpr (std::is_same_v<OP, op_acos>)  return detail::outward_lo_keep_zero(eve::acos(first));
    else if constexpr (std::is_same_v<OP, op_atan>)  return detail::outward_lo_keep_zero(eve::atan(first));
    else if constexpr (std::is_same_v<OP, op_sinh>)  return detail::outward_lo_keep_zero(eve::sinh(first));
    else if constexpr (std::is_same_v<OP, op_cosh>)  return detail::outward_lo(eve::cosh(first));
    else if constexpr (std::is_same_v<OP, op_tanh>)  return detail::outward_lo_keep_zero(eve::tanh(first));
    else if constexpr (std::is_same_v<OP, op_cbrt>)  return detail::outward_lo(eve::cbrt(first));
    else if constexpr (std::is_same_v<OP, op_log1p>) return detail::outward_lo_keep_zero(eve::log1p(first));
    // floor/ceil are exact operations — no rounding adjustment needed.
    else if constexpr (std::is_same_v<OP, op_floor>) return eve::floor(first);
    else if constexpr (std::is_same_v<OP, op_ceil>)  return eve::ceil(first);
    else static_assert(sizeof(OP) == 0, "unsupported op tag for ropd");
#else
    if constexpr (std::is_same_v<OP, op_add>)        return eve::add(first, rest...);
    else if constexpr (std::is_same_v<OP, op_sub>)   return eve::sub(first, rest...);
    else if constexpr (std::is_same_v<OP, op_mul>)   return eve::mul(first, rest...);
    else if constexpr (std::is_same_v<OP, op_div>)   return eve::div(first, rest...);
    else if constexpr (std::is_same_v<OP, op_sqrt>)  return eve::sqrt(first);
    else if constexpr (std::is_same_v<OP, op_abs>)   return eve::abs(first);
    else if constexpr (std::is_same_v<OP, op_exp>)   return eve::exp(first);
    else if constexpr (std::is_same_v<OP, op_log>)   return eve::log(first);
    else if constexpr (std::is_same_v<OP, op_sin>)   return eve::sin(first);
    else if constexpr (std::is_same_v<OP, op_cos>)   return eve::cos(first);
    else if constexpr (std::is_same_v<OP, op_tan>)   return eve::tan(first);
    else if constexpr (std::is_same_v<OP, op_asin>)  return eve::asin(first);
    else if constexpr (std::is_same_v<OP, op_acos>)  return eve::acos(first);
    else if constexpr (std::is_same_v<OP, op_atan>)  return eve::atan(first);
    else if constexpr (std::is_same_v<OP, op_sinh>)  return eve::sinh(first);
    else if constexpr (std::is_same_v<OP, op_cosh>)  return eve::cosh(first);
    else if constexpr (std::is_same_v<OP, op_tanh>)  return eve::tanh(first);
    else if constexpr (std::is_same_v<OP, op_cbrt>)  return eve::cbrt(first);
    else if constexpr (std::is_same_v<OP, op_log1p>) return eve::log1p(first);
    else if constexpr (std::is_same_v<OP, op_floor>) return eve::floor(first);
    else if constexpr (std::is_same_v<OP, op_ceil>)  return eve::ceil(first);
    else static_assert(sizeof(OP) == 0, "unsupported op tag for ropd");
#endif
}

template<typename OP, typename T, typename... Args>
T ropu(T first, Args... rest)
{
#if defined(DIRECTED_ROUNDING)
    if constexpr (std::is_same_v<OP, op_add>)        return eve::add[eve::upper](first, rest...);
    else if constexpr (std::is_same_v<OP, op_sub>)   return eve::sub[eve::upper](first, rest...);
    else if constexpr (std::is_same_v<OP, op_mul>)   return eve::mul[eve::upper](first, rest...);
    else if constexpr (std::is_same_v<OP, op_div>)   return eve::div[eve::upper](first, rest...);
    else if constexpr (std::is_same_v<OP, op_sqrt>)  return eve::sqrt[eve::upper](first);
    else if constexpr (std::is_same_v<OP, op_abs>)   return eve::abs(first);
    else if constexpr (std::is_same_v<OP, op_exp>)   return detail::outward_hi(eve::exp(first));
    else if constexpr (std::is_same_v<OP, op_log>)   return detail::outward_hi_keep_zero(eve::log(first));
    else if constexpr (std::is_same_v<OP, op_sin>)   return detail::outward_hi_keep_zero(eve::sin(first));
    else if constexpr (std::is_same_v<OP, op_cos>)   return detail::outward_hi(eve::cos(first));
    else if constexpr (std::is_same_v<OP, op_tan>)   return detail::outward_hi_keep_zero(eve::tan(first));
    else if constexpr (std::is_same_v<OP, op_asin>)  return detail::outward_hi_keep_zero(eve::asin(first));
    else if constexpr (std::is_same_v<OP, op_acos>)  return detail::outward_hi_keep_zero(eve::acos(first));
    else if constexpr (std::is_same_v<OP, op_atan>)  return detail::outward_hi_keep_zero(eve::atan(first));
    else if constexpr (std::is_same_v<OP, op_sinh>)  return detail::outward_hi_keep_zero(eve::sinh(first));
    else if constexpr (std::is_same_v<OP, op_cosh>)  return detail::outward_hi(eve::cosh(first));
    else if constexpr (std::is_same_v<OP, op_tanh>)  return detail::outward_hi_keep_zero(eve::tanh(first));
    else if constexpr (std::is_same_v<OP, op_cbrt>)  return detail::outward_hi(eve::cbrt(first));
    else if constexpr (std::is_same_v<OP, op_log1p>) return detail::outward_hi_keep_zero(eve::log1p(first));
    else if constexpr (std::is_same_v<OP, op_floor>) return eve::floor(first);
    else if constexpr (std::is_same_v<OP, op_ceil>)  return eve::ceil(first);
    else static_assert(sizeof(OP) == 0, "unsupported op tag for ropu");
#else
    return ropd<OP>(first, rest...);
#endif
}

namespace trig {
    template<std::floating_point T>
    int get_quadrant(T a)
    {
        auto x = std::fmod(a, two_pi_v<T>);
        if (x == T(0)) return 0;
        if (x < T(0)) x += two_pi_v<T>;
        return static_cast<int>(x / half_pi_v<T>);
    }
}

} // namespace pappus::fp

#endif
