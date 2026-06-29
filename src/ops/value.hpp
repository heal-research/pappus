#ifndef PAPPUS_OPS_VALUE_HPP
#define PAPPUS_OPS_VALUE_HPP

#include <concepts>

#include "context.hpp"
#include "affine/affine.hpp"
#include "interval/interval.hpp"

namespace pappus::ops {

template<std::floating_point T>
using interval_value = interval<T>;

template<std::floating_point T>
using affine_value = affine_form<T>;

template<std::floating_point T>
inline interval<T> constant(T value)
{
    return interval<T>(value);
}

template<std::floating_point T>
inline interval<T> variable(T lower, T upper)
{
    return interval<T>(lower, upper);
}

template<std::floating_point T>
inline affine_form<T> constant(affine_context<T>& context, T value)
{
    return affine_form<T>(context.state, value);
}

template<std::floating_point T>
inline affine_form<T> variable(affine_context<T>& context, T lower, T upper)
{
    return affine_form<T>(context.state, interval<T>(lower, upper));
}

} // namespace pappus::ops

#endif
