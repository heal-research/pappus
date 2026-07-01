#ifndef PAPPUS_OPS_FUNCTIONS_HPP
#define PAPPUS_OPS_FUNCTIONS_HPP

#include <concepts>
#include <utility>

#include "value.hpp"

namespace pappus::ops {

template<std::floating_point T>
inline affine_form<T> finalize(affine_context<T> const& context, affine_form<T> value)
{
    if (context.max_terms != 0 && value.size() > context.max_terms)
        value.condense(context.max_terms);
    return value;
}

template<std::floating_point T>
inline interval<T> add(interval<T> lhs, interval<T> const& rhs)
{
    lhs += rhs;
    return lhs;
}

template<std::floating_point T>
inline affine_form<T> add(affine_form<T> lhs, affine_form<T> const& rhs)
{
    lhs += rhs;
    return lhs;
}

template<std::floating_point T>
inline affine_form<T> add(affine_context<T> const& context, affine_form<T> lhs, affine_form<T> const& rhs)
{
    return finalize(context, add(std::move(lhs), rhs));
}

template<std::floating_point T>
inline interval<T> sub(interval<T> lhs, interval<T> const& rhs)
{
    lhs -= rhs;
    return lhs;
}

template<std::floating_point T>
inline affine_form<T> sub(affine_form<T> lhs, affine_form<T> const& rhs)
{
    lhs -= rhs;
    return lhs;
}

template<std::floating_point T>
inline affine_form<T> sub(affine_context<T> const& context, affine_form<T> lhs, affine_form<T> const& rhs)
{
    return finalize(context, sub(std::move(lhs), rhs));
}

template<std::floating_point T>
inline interval<T> mul(interval<T> lhs, interval<T> const& rhs)
{
    lhs *= rhs;
    return lhs;
}

template<std::floating_point T>
inline affine_form<T> mul(affine_form<T> lhs, affine_form<T> const& rhs)
{
    lhs *= rhs;
    return lhs;
}

template<std::floating_point T>
inline affine_form<T> mul(affine_context<T> const& context, affine_form<T> lhs, affine_form<T> const& rhs)
{
    return finalize(context, mul(std::move(lhs), rhs));
}

template<std::floating_point T>
inline interval<T> div(interval<T> lhs, interval<T> const& rhs)
{
    lhs /= rhs;
    return lhs;
}

template<std::floating_point T>
inline affine_form<T> div(affine_form<T> lhs, affine_form<T> const& rhs)
{
    lhs /= rhs;
    return lhs;
}

template<std::floating_point T>
inline affine_form<T> div(affine_context<T> const& context, affine_form<T> lhs, affine_form<T> const& rhs)
{
    return finalize(context, div(std::move(lhs), rhs));
}

template<std::floating_point T>
inline interval<T> neg(interval<T> value)
{
    return -value;
}

template<std::floating_point T>
inline affine_form<T> neg(affine_form<T> value)
{
    return -value;
}

template<std::floating_point T>
inline interval<T> inv(interval<T> value)
{
    return value.inv();
}

template<std::floating_point T>
inline affine_form<T> inv(affine_form<T> value)
{
    return value.inv();
}

template<std::floating_point T>
inline affine_form<T> inv(affine_context<T> const& context, affine_form<T> value)
{
    return finalize(context, value.inv());
}

template<std::floating_point T>
inline interval<T> square(interval<T> value)
{
    return value.square();
}

template<std::floating_point T>
inline affine_form<T> square(affine_form<T> value)
{
    return value * value;
}

template<std::floating_point T>
inline affine_form<T> square(affine_context<T> const& context, affine_form<T> value)
{
    return finalize(context, value * value);
}

#define PAPPUS_DEFINE_UNARY_OP(name) \
template<std::floating_point T> \
inline interval<T> name(interval<T> value) \
{ \
    return value.name(); \
} \
template<std::floating_point T> \
inline affine_form<T> name(affine_form<T> value) \
{ \
    return value.name(); \
} \
template<std::floating_point T> \
inline affine_form<T> name(affine_context<T> const& context, affine_form<T> value) \
{ \
    return finalize(context, value.name()); \
}

PAPPUS_DEFINE_UNARY_OP(sqrt)
PAPPUS_DEFINE_UNARY_OP(exp)
PAPPUS_DEFINE_UNARY_OP(log)
PAPPUS_DEFINE_UNARY_OP(sin)
PAPPUS_DEFINE_UNARY_OP(cos)
PAPPUS_DEFINE_UNARY_OP(tan)
PAPPUS_DEFINE_UNARY_OP(asin)
PAPPUS_DEFINE_UNARY_OP(acos)
PAPPUS_DEFINE_UNARY_OP(atan)
PAPPUS_DEFINE_UNARY_OP(sinh)
PAPPUS_DEFINE_UNARY_OP(cosh)
PAPPUS_DEFINE_UNARY_OP(tanh)
PAPPUS_DEFINE_UNARY_OP(abs)
PAPPUS_DEFINE_UNARY_OP(cbrt)
PAPPUS_DEFINE_UNARY_OP(log1p)
PAPPUS_DEFINE_UNARY_OP(floor)
PAPPUS_DEFINE_UNARY_OP(ceil)

#undef PAPPUS_DEFINE_UNARY_OP

// min/max: element-wise interval min/max (not hull/intersection).
// min([a1,b1], [a2,b2]) = [min(a1,a2), min(b1,b2)]
// max([a1,b1], [a2,b2]) = [max(a1,a2), max(b1,b2)]
template<std::floating_point T>
inline interval<T> min(interval<T> const& lhs, interval<T> const& rhs)
{
    return interval<T>(std::fmin(lhs.inf(), rhs.inf()), std::fmin(lhs.sup(), rhs.sup()));
}

template<std::floating_point T>
inline interval<T> max(interval<T> const& lhs, interval<T> const& rhs)
{
    return interval<T>(std::fmax(lhs.inf(), rhs.inf()), std::fmax(lhs.sup(), rhs.sup()));
}

// Affine min/max: the result's interval enclosure is the element-wise min/max
// of the input intervals. The affine approximation uses the secant with a
// delta covering the deviation. For domains where one input dominates, the
// result is simply that input.
template<std::floating_point T>
inline affine_form<T> min(affine_form<T> const& lhs, affine_form<T> const& rhs)
{
    auto a1 = lhs.min(), b1 = lhs.max();
    auto a2 = rhs.min(), b2 = rhs.max();
    // If ranges don't overlap, the result is the smaller form.
    if (b1 <= a2) return lhs;
    if (b2 <= a1) return rhs;
    // Overlapping: return the hull of the min enclosure as a constant
    // approximation (conservative but sound).
    auto lo = std::fmin(a1, a2);
    auto hi = std::fmin(b1, b2);
    return affine_form<T>(lhs.context(), interval<T>(lo, hi));
}

template<std::floating_point T>
inline affine_form<T> min(affine_context<T> const& context, affine_form<T> const& lhs, affine_form<T> const& rhs)
{
    return finalize(context, min(lhs, rhs));
}

template<std::floating_point T>
inline affine_form<T> max(affine_form<T> const& lhs, affine_form<T> const& rhs)
{
    auto a1 = lhs.min(), b1 = lhs.max();
    auto a2 = rhs.min(), b2 = rhs.max();
    if (b1 <= a2) return rhs;
    if (b2 <= a1) return lhs;
    auto lo = std::fmax(a1, a2);
    auto hi = std::fmax(b1, b2);
    return affine_form<T>(lhs.context(), interval<T>(lo, hi));
}

template<std::floating_point T>
inline affine_form<T> max(affine_context<T> const& context, affine_form<T> const& lhs, affine_form<T> const& rhs)
{
    return finalize(context, max(lhs, rhs));
}

template<std::floating_point T>
inline interval<T> pow(interval<T> base, int exponent)
{
    return base.pow(exponent);
}

template<std::floating_point T>
inline interval<T> pow(interval<T> base, T exponent)
{
    return base.pow(exponent);
}

template<std::floating_point T>
inline interval<T> pow(interval<T> base, interval<T> const& exponent)
{
    return base.pow(exponent);
}

template<std::floating_point T>
inline interval<T> pow(T base, interval<T> const& exponent)
{
    return interval<T>(base).pow(exponent);
}

template<std::floating_point T>
inline affine_form<T> pow(affine_form<T> base, int exponent)
{
    return base.pow(exponent);
}

template<std::floating_point T>
inline affine_form<T> pow(affine_context<T> const& context, affine_form<T> base, int exponent)
{
    return finalize(context, base.pow(exponent));
}

template<std::floating_point T>
inline affine_form<T> pow(affine_form<T> base, T exponent)
{
    return base.pow(exponent);
}

template<std::floating_point T>
inline affine_form<T> pow(affine_context<T> const& context, affine_form<T> base, T exponent)
{
    return finalize(context, base.pow(exponent));
}

template<std::floating_point T>
inline affine_form<T> pow(affine_form<T> base, affine_form<T> const& exponent)
{
    return base.pow(exponent);
}

template<std::floating_point T>
inline affine_form<T> pow(affine_context<T> const& context, affine_form<T> base, affine_form<T> const& exponent)
{
    return finalize(context, base.pow(exponent));
}

template<std::floating_point T>
inline affine_form<T> pow(T base, affine_form<T> const& exponent)
{
    return affine_form<T>::pow(base, exponent);
}

template<std::floating_point T>
inline affine_form<T> pow(affine_context<T> const& context, T base, affine_form<T> const& exponent)
{
    return finalize(context, affine_form<T>::pow(base, exponent));
}

// ---------------------------------------------------------------------------
// Composite ops: aq, sqrtabs, logabs
// ---------------------------------------------------------------------------
// Each reduces N finalize() calls (one per intermediate ops::xxx call) to 1.
// With max_terms=0 (default) finalize is a no-op, so the gain is mainly
// avoiding the intermediate affine_form heap allocation in the affine case.

// aq(x, y) = x / sqrt(1 + y^2)
// Affine: y is taken by value so *=y is a self-multiply (safe: operator*=
// builds a temporary, then swaps).  The +=1 shift is a center-only update
// (no new noise term).  One finalize at the end vs three in the naive path.
template<std::floating_point T>
inline interval<T> aq(interval<T> x, interval<T> const& y)
{
    return x / (y.square() + interval<T>(T(1))).sqrt();
}

template<std::floating_point T>
inline affine_form<T> aq(affine_form<T> x, affine_form<T> y)
{
    y *= y;          // y^2 in-place
    y += T(1);       // 1 + y^2: center shift, no new noise term
    // 1+y_input^2 ≥ 1 always, but the affine form can overapproximate below 0
    // for wide inputs.  When that happens, sqrt() would throw.  Fall back to an
    // interval-arithmetic denominator: sqrt([1, y.max()]) ⊇ sqrt(1+y_input^2).
    if (y.min() <= T(0)) {
        auto denom = affine_form<T>(y.context(), interval<T>(T(1), y.max()).sqrt());
        return x / denom;
    }
    return x / y.sqrt();
}

template<std::floating_point T>
inline affine_form<T> aq(affine_context<T> const& context, affine_form<T> x, affine_form<T> y)
{
    return finalize(context, aq(std::move(x), std::move(y)));
}

// sqrtabs(x) = sqrt(|x|)
// Reduces 2 finalize calls (abs then sqrt) to 1.
template<std::floating_point T>
inline interval<T> sqrtabs(interval<T> value)
{
    return value.abs().sqrt();
}

template<std::floating_point T>
inline affine_form<T> sqrtabs(affine_form<T> value)
{
    auto lo = value.min(), hi = value.max();
    if (lo >= T(0)) return value.sqrt();
    if (hi <= T(0)) return (-value).sqrt();
    // Zero-crossing: abs() overapproximates with negative min(), causing sqrt to throw.
    // Use a fresh interval [0, max(|lo|, |hi|)] for the absolute value enclosure.
    auto bound = std::fmax(-lo, hi);
    return affine_form<T>(value.context(), interval<T>(T(0), bound)).sqrt();
}

template<std::floating_point T>
inline affine_form<T> sqrtabs(affine_context<T> const& context, affine_form<T> value)
{
    return finalize(context, sqrtabs(std::move(value)));
}

// logabs(x) = log(|x|)
// Reduces 2 finalize calls (abs then log) to 1.
template<std::floating_point T>
inline interval<T> logabs(interval<T> value)
{
    return value.abs().log();
}

template<std::floating_point T>
inline affine_form<T> logabs(affine_form<T> value)
{
    return value.abs().log();
}

template<std::floating_point T>
inline affine_form<T> logabs(affine_context<T> const& context, affine_form<T> value)
{
    return finalize(context, value.abs().log());
}

} // namespace pappus::ops

#endif
