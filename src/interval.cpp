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
    return interval(fp::ropd<fp::op_add>(a1, b1), fp::ropu<fp::op_add>(a2, b2));
}

interval interval::operator-(interval const& other) const
{
    auto [a1, a2] = bounds();
    auto [b1, b2] = other.bounds();
    return interval(fp::ropd<fp::op_sub>(a1, b2), fp::ropu<fp::op_sub>(a2, b1));
}

// implementation based on the algorithm presented in:
// Interval Arithmetic from Principles to Implementation
// Timothy J. Hickey, Qun Ju, and Maarten H. van Emden
// J. ACM 48:5, p. 1038--1068, sep. 2001
interval interval::operator*(interval const& other) const
{
    auto [a1, a2] = bounds();
    auto [b1, b2] = other.bounds();

    const auto D = fp::ropd<fp::op_mul>;
    const auto U = fp::ropu<fp::op_mul>;

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

    return interval(fp::nan, fp::nan);
}

interval interval::operator/(interval const& other) const
{
    if (isempty())
        return interval::emptyset();

    EXPECT(!(std::isnan(lower()) || std::isnan(upper())));

    auto [a1, a2] = bounds();
    auto [b1, b2] = other.bounds();

    const auto D = fp::ropd<fp::op_div>;
    const auto U = fp::ropu<fp::op_div>;

    if (isinfinite() || other.isinfinite())
        return interval::unbounded();

    if (other.is_zero()) {
            return contains(0.0) ? interval::unbounded()
                                 : interval::emptyset();
    }

    if (other.contains(0.0)) {
        if (contains(0.0))
            return interval::unbounded();

        if (a2 < 0) {
            if (other.is_zero())
                return interval::emptyset(); // empty interval

            if (b2 == 0)
                return interval(D(a2, b1), fp::inf);

            if (b1 == 0)
                return interval(-fp::inf, U(a2, b2));
        }
        return interval::unbounded();
    } else {
        if (is_zero())
            return interval::emptyset();

        if (a2 <= 0) {
            if (b2 < 0)
                return interval(D(a2, b1), U(a1, b2));

            if (b1 > 0)
                return interval(D(a1, b1), U(a2, b2));
        }

        if (contains_strict(0.0)) {
            if (b2 < 0)
                return interval(D(a2, b2), U(a1, b2));

            if (b1 > 0)
                return interval(D(a1, b1), U(a2, b1));
        }

        if (0 <= a1) {
            if (b2 < 0)
                return interval(D(a2, b2), U(a1, b1));

            if (b1 > 0)
                return interval(D(a1, b2), U(a2, b1));
        }
    }

    return interval(fp::nan, fp::nan);
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
