#ifndef PAPPUS_OPS_DOMAIN_HPP
#define PAPPUS_OPS_DOMAIN_HPP

#include <cmath>
#include <concepts>
#include <optional>

#include "pappus/affine/affine.hpp"
#include "pappus/fp/util.hpp"
#include "pappus/interval/interval.hpp"

// Domain predicates and non-throwing wrappers for affine_form operations.
//
// Affine arithmetic requires finite Chebyshev bounds, so functions with
// restricted domains throw when the input interval violates them.  The
// interval<T> counterparts instead return empty() or a clamped result.
// These utilities let callers check safety before evaluation, or opt into
// a nullopt-on-failure style without exceptions.

namespace pappus {

// ---------------------------------------------------------------------------
// Domain predicates on interval<T>
// Return true iff calling the named affine_form member on a form whose
// bounding interval equals `iv` will not throw a domain error.
// ---------------------------------------------------------------------------

template<std::floating_point T>
bool log_domain_ok(interval<T> const& iv) noexcept
{
    return !iv.is_empty() && iv.inf() > T(0);
}

template<std::floating_point T>
bool log1p_domain_ok(interval<T> const& iv) noexcept
{
    return !iv.is_empty() && iv.inf() > T(-1);
}

template<std::floating_point T>
bool sqrt_domain_ok(interval<T> const& iv) noexcept
{
    return !iv.is_empty() && iv.inf() >= T(0);
}

template<std::floating_point T>
bool isqrt_domain_ok(interval<T> const& iv) noexcept
{
    return !iv.is_empty() && iv.inf() > T(0);
}

template<std::floating_point T>
bool inv_domain_ok(interval<T> const& iv) noexcept
{
    return !iv.is_empty() && (iv.inf() > T(0) || iv.sup() < T(0));
}

template<std::floating_point T>
bool tan_domain_ok(interval<T> const& iv) noexcept
{
    if (iv.is_empty()) return false;
    auto hp = fp::half_pi_v<T>, pi = fp::pi_v<T>;
    return static_cast<long>(std::floor((iv.inf() + hp) / pi))
        == static_cast<long>(std::floor((iv.sup() + hp) / pi));
}

template<std::floating_point T>
bool asin_domain_ok(interval<T> const& iv) noexcept
{
    return !iv.is_empty() && iv.inf() >= T(-1) && iv.sup() <= T(1);
}

template<std::floating_point T>
bool acos_domain_ok(interval<T> const& iv) noexcept
{
    return !iv.is_empty() && iv.inf() >= T(-1) && iv.sup() <= T(1);
}

template<std::floating_point T>
bool pow_domain_ok(interval<T> const& iv, T exponent) noexcept
{
    if (iv.is_empty()) return false;
    bool is_int = std::isfinite(exponent) && std::trunc(exponent) == exponent;
    if (!is_int && iv.inf() < T(0)) return false;
    if (exponent < T(0) && iv.inf() <= T(0)) return false;
    return true;
}

// ---------------------------------------------------------------------------
// affine_form overloads — delegate to to_interval()
// ---------------------------------------------------------------------------

template<std::floating_point T> bool log_domain_ok  (affine_form<T> const& x) noexcept { return log_domain_ok  (x.to_interval()); }
template<std::floating_point T> bool log1p_domain_ok (affine_form<T> const& x) noexcept { return log1p_domain_ok(x.to_interval()); }
template<std::floating_point T> bool sqrt_domain_ok  (affine_form<T> const& x) noexcept { return sqrt_domain_ok (x.to_interval()); }
template<std::floating_point T> bool isqrt_domain_ok (affine_form<T> const& x) noexcept { return isqrt_domain_ok(x.to_interval()); }
template<std::floating_point T> bool inv_domain_ok   (affine_form<T> const& x) noexcept { return inv_domain_ok  (x.to_interval()); }
template<std::floating_point T> bool tan_domain_ok   (affine_form<T> const& x) noexcept { return tan_domain_ok  (x.to_interval()); }
template<std::floating_point T> bool asin_domain_ok  (affine_form<T> const& x) noexcept { return asin_domain_ok (x.to_interval()); }
template<std::floating_point T> bool acos_domain_ok  (affine_form<T> const& x) noexcept { return acos_domain_ok (x.to_interval()); }
template<std::floating_point T> bool pow_domain_ok   (affine_form<T> const& x, T e) noexcept { return pow_domain_ok(x.to_interval(), e); }

// ---------------------------------------------------------------------------
// Non-throwing wrappers — return nullopt if the domain check fails
// ---------------------------------------------------------------------------

template<std::floating_point T>
std::optional<affine_form<T>> try_log(affine_form<T> const& x)
{
    if (!log_domain_ok(x)) return std::nullopt;
    return x.log();
}

template<std::floating_point T>
std::optional<affine_form<T>> try_log1p(affine_form<T> const& x)
{
    if (!log1p_domain_ok(x)) return std::nullopt;
    return x.log1p();
}

template<std::floating_point T>
std::optional<affine_form<T>> try_sqrt(affine_form<T> const& x)
{
    if (!sqrt_domain_ok(x)) return std::nullopt;
    return x.sqrt();
}

template<std::floating_point T>
std::optional<affine_form<T>> try_isqrt(affine_form<T> const& x)
{
    if (!isqrt_domain_ok(x)) return std::nullopt;
    return x.isqrt();
}

template<std::floating_point T>
std::optional<affine_form<T>> try_inv(affine_form<T> const& x)
{
    if (!inv_domain_ok(x)) return std::nullopt;
    return x.inv();
}

template<std::floating_point T>
std::optional<affine_form<T>> try_tan(affine_form<T> const& x)
{
    if (!tan_domain_ok(x)) return std::nullopt;
    return x.tan();
}

template<std::floating_point T>
std::optional<affine_form<T>> try_asin(affine_form<T> const& x)
{
    if (!asin_domain_ok(x)) return std::nullopt;
    return x.asin();
}

template<std::floating_point T>
std::optional<affine_form<T>> try_acos(affine_form<T> const& x)
{
    if (!acos_domain_ok(x)) return std::nullopt;
    return x.acos();
}

template<std::floating_point T>
std::optional<affine_form<T>> try_pow(affine_form<T> const& x, T exponent)
{
    if (!pow_domain_ok(x, exponent)) return std::nullopt;
    return x.pow(exponent);
}

// ---------------------------------------------------------------------------
// safe_* wrappers — mirror interval<T> behaviour
//
// Unlike try_*, these clamp the input to the valid subdomain before computing,
// matching what interval<T> does implicitly.  The affine correlation of the
// original form is lost on clamping — the clamped form is a fresh noise
// symbol.  Return nullopt only when interval<T> would return empty(), i.e.
// when the interval is entirely outside the valid domain.
//
// inv and tan are excluded: when the domain is violated the true result is
// unbounded, which affine arithmetic cannot represent.
// ---------------------------------------------------------------------------

namespace {
    template<std::floating_point T>
    constexpr T pos_min() noexcept { return std::numeric_limits<T>::min(); }
}

// When clamping, the lower bound must survive the affine form's center-radius
// arithmetic: min() = (lo+hi)/2 - (hi-lo)/2.  If lo << hi the subtraction
// cancels and min() rounds to 0.  We use a relative epsilon (4 ULPs of hi)
// that is guaranteed to survive, then verify with the predicate after construction.

template<std::floating_point T>
std::optional<affine_form<T>> safe_log(affine_form<T> const& x)
{
    auto iv = x.to_interval();
    if (iv.is_empty() || iv.sup() <= T(0)) return std::nullopt;
    if (iv.inf() > T(0)) return x.log();
    auto lo = std::max(pos_min<T>(),
        std::ldexp(iv.sup(), -(std::numeric_limits<T>::digits - 2)));
    affine_form<T> clamped(x.context(), interval<T>(lo, iv.sup()));
    if (!log_domain_ok(clamped)) return std::nullopt;
    return clamped.log();
}

template<std::floating_point T>
std::optional<affine_form<T>> safe_log1p(affine_form<T> const& x)
{
    auto iv = x.to_interval();
    if (iv.is_empty() || iv.sup() <= T(-1)) return std::nullopt;
    if (iv.inf() > T(-1)) return x.log1p();
    auto span = iv.sup() - T(-1);
    auto lo = T(-1) + std::max(pos_min<T>(),
        std::ldexp(span, -(std::numeric_limits<T>::digits - 2)));
    affine_form<T> clamped(x.context(), interval<T>(lo, iv.sup()));
    if (!log1p_domain_ok(clamped)) return std::nullopt;
    return clamped.log1p();
}

template<std::floating_point T>
std::optional<affine_form<T>> safe_sqrt(affine_form<T> const& x)
{
    auto iv = x.to_interval();
    if (iv.is_empty() || iv.sup() < T(0)) return std::nullopt;
    if (iv.inf() >= T(0)) return x.sqrt();
    return affine_form<T>(x.context(), interval<T>(T(0), iv.sup())).sqrt();
}

template<std::floating_point T>
std::optional<affine_form<T>> safe_isqrt(affine_form<T> const& x)
{
    auto iv = x.to_interval();
    if (iv.is_empty() || iv.sup() <= T(0)) return std::nullopt;
    if (iv.inf() > T(0)) return x.isqrt();
    return affine_form<T>(x.context(), interval<T>(pos_min<T>(), iv.sup())).isqrt();
}

template<std::floating_point T>
std::optional<affine_form<T>> safe_asin(affine_form<T> const& x)
{
    auto iv = x.to_interval();
    if (iv.is_empty()) return std::nullopt;
    auto lo = std::max(iv.inf(), T(-1)), hi = std::min(iv.sup(), T(1));
    if (lo > hi) return std::nullopt;
    if (lo == iv.inf() && hi == iv.sup()) return x.asin();
    return affine_form<T>(x.context(), interval<T>(lo, hi)).asin();
}

template<std::floating_point T>
std::optional<affine_form<T>> safe_acos(affine_form<T> const& x)
{
    auto iv = x.to_interval();
    if (iv.is_empty()) return std::nullopt;
    auto lo = std::max(iv.inf(), T(-1)), hi = std::min(iv.sup(), T(1));
    if (lo > hi) return std::nullopt;
    if (lo == iv.inf() && hi == iv.sup()) return x.acos();
    return affine_form<T>(x.context(), interval<T>(lo, hi)).acos();
}

} // namespace pappus

#endif
