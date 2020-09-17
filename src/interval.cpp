#include "interval.hpp"

namespace pappus {

interval interval::operator-() const
{
    return interval(-upper(), -lower());
}

interval interval::operator+(interval const& other) const
{
    auto [a1, a2] = bounds();
    auto [b1, b2] = other.bounds();
    return interval(ropd<op_add>(a1, b1), ropu<op_add>(a2, b2));
}

interval interval::operator-(interval const& other) const
{
    auto [a1, a2] = bounds();
    auto [b1, b2] = other.bounds();
    return interval(ropd<op_sub>(a1, b2), ropu<op_sub>(a2, b1));
}

// implementation based on the algorithm presented in:
// Interval Arithmetic from Principles to Implementation
// Timothy J. Hickey, Qun Ju, and Maarten H. van Emden
// J. ACM 48:5, p. 1038--1068, sep. 2001
interval interval::operator*(interval const& other) const
{
    auto [a1, a2] = bounds();
    auto [b1, b2] = other.bounds();

    const auto D = ropd<op_mul>;
    const auto U = ropu<op_mul>;

    if (is_zero() || other.is_zero())
        return interval(0.0, 0.0);

    if (a2 <= 0) {
        if (b2 <= 0)
            return interval(D(a2, b2), U(a1, b1));

        if (other.contains_strict(0.0))
            return interval(D(a1, b2), U(a1, b1));

        if (b1 >= 0)
            return interval(D(a1, b2), U(a2, b1));
    } else if (contains_strict(0.0)) {
        if (b2 <= 0)
            return interval(D(a2, b1), U(a1, b1));

        if (other.contains_strict(0.0))
            return interval(std::min(D(a1, b2), D(a2, b1)), std::min(U(a1, b1), U(a2, b2)));

        if (b1 >= 0)
            return interval(D(a1, b2), D(a2, b2));
    } else if (a1 >= 0) {
        if (b2 <= 0)
            return interval(D(a2, b1), U(a1, b2));

        if (other.contains_strict(0.0))
            return interval(D(a2, b1), U(a2, b2));

        if (b1 >= 0)
            return interval(D(a1, b1), U(a2, b2));
    }

    return interval(NAN, NAN);
}

interval interval::operator/(interval const& other) const
{
    auto [a1, a2] = bounds();
    auto [b1, b2] = other.bounds();

    const auto D = ropd<op_div>;
    const auto U = ropu<op_div>;

    if (isinfinite() || other.isinfinite())
        return interval(-INFINITY, INFINITY);

    if (other.is_zero()) {
            return contains(0.0) ? interval(-INFINITY, INFINITY)
                                 : interval(NAN, NAN);
    }

    if (other.contains(0.0)) {
        if (contains(0.0))
            return interval(-INFINITY, INFINITY);

        if (a2 < 0) {
            if (other.is_zero())
                return interval(NAN, NAN); // empty interval

            if (b2 == 0)
                return interval(D(a2, b1), INFINITY);

            if (b1 == 0)
                return interval(-INFINITY, D(a2, b2));

            // in the case b1 < 0 < b2 we could return [-inf, inf]
            // or two distinct intervals: [-inf, U(a2, b2)] and [D(a2, b1), inf]
            // we choose the simple variant for now
            return interval(-INFINITY, INFINITY);
        }

    } else {

    }

    return interval(NAN, NAN);
}

interval& interval::operator+=(interval const& other)
{
    auto tmp = *this + other;
    std::swap(*this, tmp);
    return *this;
}

interval& interval::operator-=(interval const& other)
{
    auto tmp = *this - other;
    std::swap(*this, tmp);
    return *this;
}
}
