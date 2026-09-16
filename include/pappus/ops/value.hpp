#ifndef PAPPUS_OPS_VALUE_HPP
#define PAPPUS_OPS_VALUE_HPP

#include "context.hpp"
#include "pappus/affine/affine.hpp"
#include "pappus/interval/interval.hpp"

namespace pappus::ops {

template<eve::floating_value T>
using interval_value = interval<T>;

template<eve::floating_value T>
using affine_value = affine_form<T>;

template<eve::floating_value T>
inline interval<T> constant(T value)
{
    return interval<T>(value);
}

template<eve::floating_value T>
inline interval<T> variable(T lower, T upper)
{
    return interval<T>(lower, upper);
}

template<eve::floating_value T>
inline affine_form<T> constant(affine_context<T>& context, T value)
{
    return affine_form<T>(context.state, value);
}

template<eve::floating_value T>
inline affine_form<T> variable(affine_context<T>& context, T lower, T upper)
{
    return affine_form<T>(context.state, interval<T>(lower, upper));
}

} // namespace pappus::ops

#endif
