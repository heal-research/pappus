#ifndef PAPPUS_FPOPS_HPP
#define PAPPUS_FPOPS_HPP

// Directed-rounding primitives for sound interval arithmetic.
//
// Exact operations (+ - * sqrt) use eve's `lower`/`upper` decorators (no
// rounding-mode register involved: AVX-512 encodes direction as an
// instruction operand; elsewhere eve uses an error-free transform + ULP
// nudge). Division is the exception: eve::div[lower] is unsound (see
// detail::div_rd below) and is computed as -div[upper](-a, b) instead.
//
// Transcendentals use eve's round-to-nearest kernels, which are NOT
// correctly rounded: the measured worst error is ~3.6 ULP (cbrt), not
// 0.5-1 ULP (see detail::eve_ulp_margin). The result is therefore widened by a
// fixed number of ULPs of the *evaluation* type:
//   - double: eve::<fn>(x) widened by eve_ulp_margin double ULPs;
//   - float : eve::<fn>(double(x)) (eve's double kernels, error well below
//             one float ULP after the narrowing conversion) widened by
//             eve_ulp_margin double ULPs, then converted to float rounding
//             toward -inf / +inf.
// The result is a sound but slightly wider enclosure than a correctly
// rounded implementation would produce: ~2-3 ULP(float), ~16 ULP(double).
//
// When DIRECTED_ROUNDING is not defined, ropd/ropu degenerate to plain
// round-to-nearest computation: fast but not sound. Useful for performance
// testing or when soundness is guaranteed by other means.

#include "util.hpp"

#include <cmath>
#include <cstdint>
#include <concepts>
#include <limits>
#include <type_traits>

#include <eve/eve.hpp>
#include <eve/module/core.hpp>
#include <eve/module/math.hpp>

namespace pappus::fp {

namespace detail {

// Number of ULPs (of the evaluation type) added on each side of an eve
// transcendental. MPFR sweeps found a maximum 3.6-ULP error in the exact
// result's grid (cbrt). `prev`/`next` step in the kernel result's grid,
// which can be half as large at a binade boundary, so 16 steps preserve an
// 8-exact-ULP safety margin.
inline constexpr int eve_ulp_margin = 16;

// Number of ULPs added around a libm (std::exp, std::cbrt, ...) value by
// widen_lo / widen_hi. glibc documents up to 4 ULP (cbrt, double) and was
// measured at 3.3 ULP (float: below 1); a correctly rounded libm needs 1.
inline constexpr int libm_ulp_margin = 8;

// Move y by k ULPs toward -inf / +inf. NaN and infinities are preserved;
// callers clamp an unsafe infinity (see transc).
template<typename T>
T step_dn(T y, int k) { return eve::if_else(eve::is_not_finite(y), y, eve::prev(y, k)); }

template<typename T>
T step_up(T y, int k) { return eve::if_else(eve::is_not_finite(y), y, eve::next(y, k)); }

// ---------------------------------------------------------------------------
// Exact operations: eve's [lower]/[upper] plus two range-edge repairs.
//
// eve implements directed + - * / sqrt as round-to-nearest followed by an
// error-free transform (two_sum / two_prod via fma) that decides whether to
// step one ULP. Two situations break that, measured against MPFR:
//   * overflow: a finite-operand result beyond the largest finite value
//     comes back as +inf for the lower bound (-inf for the upper bound),
//     which is not a bound of the (finite) exact result;
//   * underflow: when the exact product/quotient/root is below ~2^(emin+p)
//     the error term itself underflows, so the transform sees "no error" and
//     returns the round-to-nearest value (0 for a tiny product rounding to
//     0), which can lie on the wrong side of the exact result.
// Both are repaired here: an unsafe-side infinity from finite operands is
// replaced by the largest finite value, and below tiny_threshold the result
// is moved one ULP outward (the round-to-nearest value and eve's result are
// within 1 ULP of the exact one, so one extra step is sufficient).
// ---------------------------------------------------------------------------

// 2^(emin + digits + 2): float 2^-100, double 2^-967.
template<typename T>
T tiny_threshold()
{
    using E = eve::element_type_t<T>;
    return T(std::ldexp(std::numeric_limits<E>::min(), std::numeric_limits<E>::digits + 2));
}

// +inf as a lower bound / -inf as an upper bound of a result computed from
// finite operands -> largest finite value.
template<bool DOWN, typename T, typename M>
T clamp_overflow(T r, M finite_operands)
{
    if constexpr (DOWN) return eve::if_else(finite_operands && r == eve::inf(eve::as(r)), eve::valmax(eve::as(r)), r);
    else                return eve::if_else(finite_operands && r == -eve::inf(eve::as(r)), eve::valmin(eve::as(r)), r);
}

template<bool DOWN, typename T, typename M>
T nudge_if(T r, M mask)
{
    if constexpr (DOWN) return eve::if_else(mask, eve::prev(r), r);
    else                return eve::if_else(mask, eve::next(r), r);
}

template<bool DOWN, typename T>
T add_rd(T a, T b)
{
    T const r = DOWN ? eve::add[eve::lower](a, b) : eve::add[eve::upper](a, b);
    return clamp_overflow<DOWN>(r, eve::is_finite(a) && eve::is_finite(b));
}

template<bool DOWN, typename T>
T sub_rd(T a, T b)
{
    T const r = DOWN ? eve::sub[eve::lower](a, b) : eve::sub[eve::upper](a, b);
    return clamp_overflow<DOWN>(r, eve::is_finite(a) && eve::is_finite(b));
}

template<bool DOWN, typename T>
T mul_rd(T a, T b)
{
    T r = DOWN ? eve::mul[eve::lower](a, b) : eve::mul[eve::upper](a, b);
    r = nudge_if<DOWN>(r, eve::abs(r) < tiny_threshold<T>() && eve::is_nez(a) && eve::is_nez(b));
    return clamp_overflow<DOWN>(r, eve::is_finite(a) && eve::is_finite(b));
}

// Quotient.
//
// eve::div[eve::lower] is broken independently of the range edges: in
// eve/module/core/regular/impl/div.hpp the lower branch tests `is_ltz(e)`
// where it must test `is_gtz(e)` (the upper branch is correct), so ~45% of
// inexact quotients come out one ULP too high, i.e. not a lower bound. Until
// that is fixed upstream the lower bound is RD(a / b) = -RU(-a / b), exact
// because negation is exact; replace the DOWN branch with
// eve::div[eve::lower](a, b) once eve is fixed. (Underflow: the error term
// of d * b underflows when |a| < tiny_threshold, see above.)
template<bool DOWN, typename T>
T div_rd(T a, T b)
{
    if constexpr (DOWN) {
        return -div_rd<false>(-a, b);
    } else {
        T r = eve::div[eve::upper](a, b);
        r = nudge_if<false>(r, eve::abs(a) < tiny_threshold<T>() && eve::is_nez(a) && eve::is_nez(b) && eve::is_finite(b));
        return clamp_overflow<false>(r, eve::is_finite(a) && eve::is_finite(b) && eve::is_nez(b));
    }
}

template<bool DOWN, typename T>
T sqrt_rd(T x)
{
    T r = DOWN ? eve::sqrt[eve::lower](x) : eve::sqrt[eve::upper](x);
    return nudge_if<DOWN>(r, x < tiny_threshold<T>() && x > T(0));
}

// Operation tags handled by transc (defined after the tag structs below).
template<typename OP> inline constexpr bool is_transcendental_v = false;

// Evaluation type of a transcendental on T: float lanes are evaluated in
// double lanes of the same cardinality, everything else in T itself.
template<typename T>
using eval_t = std::conditional_t<
    std::is_same_v<eve::element_type_t<T>, float>,
    std::conditional_t<eve::scalar_value<T>, double, eve::wide<double, eve::cardinal_t<T>>>,
    T>;

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

namespace detail {

template<> inline constexpr bool is_transcendental_v<op_exp>   = true;
template<> inline constexpr bool is_transcendental_v<op_log>   = true;
template<> inline constexpr bool is_transcendental_v<op_sin>   = true;
template<> inline constexpr bool is_transcendental_v<op_cos>   = true;
template<> inline constexpr bool is_transcendental_v<op_tan>   = true;
template<> inline constexpr bool is_transcendental_v<op_asin>  = true;
template<> inline constexpr bool is_transcendental_v<op_acos>  = true;
template<> inline constexpr bool is_transcendental_v<op_atan>  = true;
template<> inline constexpr bool is_transcendental_v<op_sinh>  = true;
template<> inline constexpr bool is_transcendental_v<op_cosh>  = true;
template<> inline constexpr bool is_transcendental_v<op_tanh>  = true;
template<> inline constexpr bool is_transcendental_v<op_cbrt>  = true;
template<> inline constexpr bool is_transcendental_v<op_log1p> = true;

// The eve kernel behind a transcendental tag.
template<typename OP, typename C>
C eval_kernel(C v)
{
    if constexpr (std::is_same_v<OP, op_exp>)        return eve::exp(v);
    else if constexpr (std::is_same_v<OP, op_log>)   return eve::log(v);
    else if constexpr (std::is_same_v<OP, op_sin>)   return eve::sin(v);
    else if constexpr (std::is_same_v<OP, op_cos>)   return eve::cos(v);
    else if constexpr (std::is_same_v<OP, op_tan>)   return eve::tan(v);
    else if constexpr (std::is_same_v<OP, op_asin>)  return eve::asin(v);
    else if constexpr (std::is_same_v<OP, op_acos>)  return eve::acos(v);
    else if constexpr (std::is_same_v<OP, op_atan>)  return eve::atan(v);
    else if constexpr (std::is_same_v<OP, op_sinh>)  return eve::sinh(v);
    else if constexpr (std::is_same_v<OP, op_cosh>)  return eve::cosh(v);
    else if constexpr (std::is_same_v<OP, op_tanh>)  return eve::tanh(v);
    else if constexpr (std::is_same_v<OP, op_cbrt>)  return eve::cbrt(v);
    else if constexpr (std::is_same_v<OP, op_log1p>) return eve::log1p(v);
    else static_assert(sizeof(OP) == 0, "not a transcendental op tag");
}

// Preserve an exact zero only when the mathematical function is zero at the
// exact input. Kernel output alone is insufficient: a kernel may flush a
// nonzero subnormal result to zero, which must still be widened outward.
template<typename OP, typename C>
auto exact_zero_input_mask(C x)
{
    if constexpr (std::is_same_v<OP, op_sin> || std::is_same_v<OP, op_tan> ||
                  std::is_same_v<OP, op_asin> || std::is_same_v<OP, op_atan> ||
                  std::is_same_v<OP, op_sinh> || std::is_same_v<OP, op_tanh> ||
                  std::is_same_v<OP, op_cbrt> || std::is_same_v<OP, op_log1p>)
        return eve::is_eqz(x);
    else if constexpr (std::is_same_v<OP, op_log> || std::is_same_v<OP, op_acos>)
        return x == C(1);
    else
        return eve::false_(eve::as<C> {});
}

// Directed-rounded transcendental: DOWN -> lower bound, else upper bound.
//   1. float is evaluated in double (eve's double kernels), everything else
//      in its own type;
//   2. the result is moved eve_ulp_margin ULPs of the evaluation type outward;
//   3. an infinity on the unsafe side (+inf lower bound / -inf upper bound)
//      produced by overflow of a finite argument is replaced by a finite
//      bound. Infinite arguments preserve their mathematical infinite limit.
//   4. for float, narrow to float rounding toward -inf / +inf (the double
//      value is exact in float or lies strictly between two floats).
template<bool DOWN, typename OP, typename T>
T transc(T x)
{
    using E = eve::element_type_t<T>;
    using C = eval_t<T>;

    C xc;
    if constexpr (std::is_same_v<E, float>) xc = eve::convert(x, eve::as<double>{});
    else                                    xc = x;

    C const yc = eval_kernel<OP>(xc);
    C b = DOWN ? step_dn(yc, eve_ulp_margin) : step_up(yc, eve_ulp_margin);

    // An unsafe-side infinity from a finite argument is not a bound of a
    // finite exact result. Keep infinities generated by infinite arguments.
    auto unsafe_inf = DOWN ? (b == eve::inf(eve::as(b))) : (b == -eve::inf(eve::as(b)));
    unsafe_inf = unsafe_inf && eve::is_finite(xc);
    if constexpr (std::is_same_v<OP, op_exp>)
        b = eve::if_else(unsafe_inf, eve::if_else(xc >= C(710), eve::valmax(eve::as(b)), C(std::ldexp(1.0, 1023))), b);
    else
        b = eve::if_else(unsafe_inf, DOWN ? eve::valmax(eve::as(b)) : eve::valmin(eve::as(b)), b);

    auto const exact_zero = exact_zero_input_mask<OP>(xc) && eve::is_eqz(yc);
    b = eve::if_else(exact_zero, yc, b);
    // eve::exp(double) flushes subnormal results to 0 (x <= -708.4), so the
    // kernel value says nothing about them. Any such result is below
    // 2^-1022 * (1 + 1e-13), hence 2^-1021 is a sound upper bound.
    if constexpr (std::is_same_v<OP, op_exp> && !DOWN) {
        constexpr auto tiny = std::numeric_limits<eve::element_type_t<C>>::min();
        b = eve::if_else(yc < tiny, eve::max(b, C(tiny * 2)), b);
    }

    if constexpr (std::is_same_v<E, float>) {
        T r = eve::convert(b, eve::as<float>{});               // round to nearest
        C const back = eve::convert(r, eve::as<double>{});
        // Nearest may have rounded the wrong way (and overflow rounds to
        // +-inf): step one float in the requested direction. prev(+inf) and
        // next(-inf) are the largest finite floats.
        if constexpr (DOWN) return eve::if_else(back > b, eve::prev(r), r);
        else                return eve::if_else(back < b, eve::next(r), r);
    } else {
        return b;
    }
}

} // namespace detail

// ---------------------------------------------------------------------------
// ropd: round-op toward -inf. ropu: round-op toward +inf.
// ---------------------------------------------------------------------------
template<typename OP, typename T, typename... Args>
T ropd(T first, Args... rest)
{
#if defined(DIRECTED_ROUNDING)
    if constexpr (std::is_same_v<OP, op_add>)        return detail::add_rd<true>(first, rest...);
    else if constexpr (std::is_same_v<OP, op_sub>)   return detail::sub_rd<true>(first, rest...);
    else if constexpr (std::is_same_v<OP, op_mul>)   return detail::mul_rd<true>(first, rest...);
    else if constexpr (std::is_same_v<OP, op_div>)   return detail::div_rd<true>(first, rest...);
    else if constexpr (std::is_same_v<OP, op_sqrt>)  return detail::sqrt_rd<true>(first);
    else if constexpr (std::is_same_v<OP, op_abs>)   return eve::abs(first);
    // Transcendentals: eve kernel + eve_ulp_margin outward expansion (see detail::transc).
    else if constexpr (detail::is_transcendental_v<OP>) return detail::transc<true, OP>(first);
    // floor/ceil are exact operations — they produce exact integers, so
    // directed rounding is a no-op. Both ropd and ropu resolve to the same
    // plain eve::floor/eve::ceil. Soundness comes from mathematical exactness,
    // not from rounding mode. The ropd/ropu wrapper is kept for API uniformity.
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
    if constexpr (std::is_same_v<OP, op_add>)        return detail::add_rd<false>(first, rest...);
    else if constexpr (std::is_same_v<OP, op_sub>)   return detail::sub_rd<false>(first, rest...);
    else if constexpr (std::is_same_v<OP, op_mul>)   return detail::mul_rd<false>(first, rest...);
    else if constexpr (std::is_same_v<OP, op_div>)   return detail::div_rd<false>(first, rest...);
    else if constexpr (std::is_same_v<OP, op_sqrt>)  return detail::sqrt_rd<false>(first);
    else if constexpr (std::is_same_v<OP, op_abs>)   return eve::abs(first);
    else if constexpr (detail::is_transcendental_v<OP>) return detail::transc<false, OP>(first);
    else if constexpr (std::is_same_v<OP, op_floor>) return eve::floor(first);
    else if constexpr (std::is_same_v<OP, op_ceil>)  return eve::ceil(first);
    else static_assert(sizeof(OP) == 0, "unsupported op tag for ropu");
#else
    return ropd<OP>(first, rest...);
#endif
}

// Outward-widening primitives for callers (e.g. affine_form's unary-op
// approximations) that need to guard a plain round-to-nearest *libm*
// value (std::exp, std::cbrt, ...) against its own error without going
// through a tagged ropd/ropu operation. libm functions are not correctly
// rounded either: widen by detail::libm_ulp_margin (8) rather than 1 ULP.
// +inf (lo) / -inf (hi) is clamped to the largest finite value, NaN and the
// safe-side infinity are preserved.
template<typename T>
T widen_lo(T x)
{
    auto y = detail::step_dn(x, detail::libm_ulp_margin);
    return eve::if_else(y == eve::inf(eve::as(y)), eve::valmax(eve::as(y)), y);
}

template<typename T>
T widen_hi(T x)
{
    auto y = detail::step_up(x, detail::libm_ulp_margin);
    return eve::if_else(y == -eve::inf(eve::as(y)), eve::valmin(eve::as(y)), y);
}

namespace trig {
    template<eve::floating_value T>
    auto get_quadrant(T a)
    {
        if constexpr (eve::value<T>) {
            auto x = eve::fmod(a, two_pi_v<T>);
            x = eve::if_else(x < T(0), x + two_pi_v<T>, x);
            // x + 2*pi rounds to exactly 2*pi for tiny negative a (|a| < ulp(2*pi)/2):
            // that is quadrant 3 (just below 2*pi), not a fifth quadrant; quadrant 4
            // would make [a, b] look like it wraps through a critical point.
            return eve::min(eve::floor(x / half_pi_v<T>), T(3));
        } else {
            auto x = std::fmod(a, two_pi_v<T>);
            if (x == T(0)) return 0;
            if (x < T(0)) x += two_pi_v<T>;
            return std::min(static_cast<int>(x / half_pi_v<T>), 3);
        }
    }

    // True iff the closed interval [a, b] (a <= b) contains a tangent
    // asymptote, i.e. a point of the form (k + 1/2)*pi for some integer k.
    // Touching a boundary counts as containing it (tan is undefined
    // there) -- the standard sound-interval-arithmetic convention of
    // treating a domain-boundary singularity as unsafe to exclude.
    //
    // Derivation: k0 is the smallest integer with (k0 + 1/2)*pi >= a; a
    // pole lies in [a, b] iff that nearest-from-above pole is <= b. This
    // replaces what were two independent, inconsistent, and (for the
    // quadrant-comparison one) demonstrably buggy tests -- interval<T>::tan()
    // used to compare get_quadrant(a) against get_quadrant(b) directly
    // (wrong for e.g. [pi/2, pi], which it treated as pole-free) while
    // affine_form<T>::tan()/tan_domain_ok used a separate floor-based test.
    // All three now share this one derivation.
    //
    // Precision caveat: pi_v<T> is itself a rounded constant, so for |a|,
    // |b| many periods from the origin, k0*pi drifts from the true
    // multiple of pi the same way std::fmod/std::floor-based reduction
    // does elsewhere in this file. Not a regression introduced by this
    // predicate -- an accepted, pre-existing limitation of doing periodic-
    // domain reduction without extended-precision pi.
    template<eve::floating_value T>
    auto has_tan_pole(T a, T b)
    {
        if constexpr (eve::value<T>) {
            auto const pi = pi_v<T>, half_pi = half_pi_v<T>;
            auto k0 = eve::ceil(a / pi - T(0.5));
            auto pole = k0 * pi + half_pi;
            auto finite_mask = eve::is_finite(a) && eve::is_finite(b) && eve::is_finite(pole);
            return eve::if_else(finite_mask, pole <= b, eve::true_(eve::as<T>{}));
        } else {
            // Non-finite bounds, or a finite `a`/`b` so large that k0*pi
            // overflows: cannot safely rule out a pole in either case (an
            // overflowed/non-finite `pole` compared via `<=` can spuriously
            // read as "false", i.e. "no pole", when the domain is actually
            // far too wide/extreme to say either way). Default to
            // conservative true -- a false "contains a pole" only costs
            // tightness (the caller rejects/treats as unbounded); a false
            // "no pole" would be a real soundness break.
            if (!std::isfinite(a) || !std::isfinite(b)) { return true; }
            auto const pi = pi_v<T>, half_pi = half_pi_v<T>;
            auto k0 = std::ceil(a / pi - T(0.5));
            auto pole = k0 * pi + half_pi;
            if (!std::isfinite(pole)) { return true; }
            return pole <= b;
        }
    }
}

} // namespace pappus::fp

#endif
