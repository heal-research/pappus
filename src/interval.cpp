#include "interval.hpp"

namespace pappus {

using namespace pappus::fp;

interval interval::operator+() const
{
    return *this;
}

interval interval::operator+(double v) const
{
    return interval(ropd<op_add>(lower(), v), ropu<op_add>(upper(), v));
}

interval interval::operator+(interval const& other) const
{
    auto [a, b] = bounds();
    auto [c, d] = other.bounds();
    return interval(ropd<op_add>(a, c), ropu<op_add>(b, d));
}

interval interval::operator-() const
{
    return interval(-upper(), -lower());
}

interval interval::operator-(double v) const
{
    return interval(ropd<op_sub>(lower(), v), ropu<op_sub>(upper(), v));
}

interval interval::operator-(interval const& other) const
{
    auto [a, b] = bounds();
    auto [c, d] = other.bounds();
    return interval(ropd<op_sub>(a, d), ropu<op_sub>(b, c));
}

interval interval::operator*(double v) const
{
    if (is_empty() || std::isnan(v))
        return interval::empty();

    if (v == 0)
        return interval::zero();

    return v < 0
        ? interval(ropd<op_mul>(upper(), v), ropu<op_mul>(lower(), v))
        : interval(ropd<op_mul>(lower(), v), ropu<op_mul>(upper(), v));
}

interval interval::operator*(interval const& other) const
{
    if (is_empty() || other.is_empty())
        return interval::empty();

    if (is_zero() || other.is_zero())
        return interval::zero();

    if (is_infinite() || other.is_infinite())
        return interval::infinite();

    auto [a, b] = bounds();
    auto [c, d] = other.bounds();

    const auto D = ropd<op_mul>;
    const auto U = ropu<op_mul>;

    if (b <= 0) { // N
        if (d <= 0) // N
            return interval(D(b, d), U(a, c));

        if (other.contains_strict(0.0)) // M
            return interval(D(a, d), U(a, c));

        if (0 <= c) // P
            return interval(D(a, d), U(b, c));
    } else if (contains_strict(0.0)) { // M
        if (d <= 0) // N
            return interval(D(b, c), U(a, c));

        if (other.contains_strict(0.0)) // M
            return interval(std::min(D(a, d), D(b, c)), std::max(U(a, c), U(b, d)));

        if (0 <= c) // P
            return interval(D(a, d), D(b, d));
    } else if (0 <= a) { // P
        if (d <= 0) // N
            return interval(D(b, c), U(a, d));

        if (other.contains_strict(0.0)) // M
            return interval(D(b, c), U(b, d));

        if (0 <= c) // P
            return interval(D(a, c), U(b, d));
    }
    return interval::empty();
}

interval interval::operator/(double v) const
{
    if (is_empty() || v == 0.0 || std::isnan(v))
        return interval::empty();

    return v < 0
        ? interval(ropd<op_div>(upper(), v), ropu<op_div>(lower(), v))
        : interval(ropd<op_div>(lower(), v), ropu<op_div>(upper(), v));
}

interval interval::operator/(interval const& other) const
{
    if (is_empty() || other.is_empty() || other.is_zero())
        return interval::empty();

    if (is_zero())
        return interval::zero();

    if (other.contains_strict(0.0))
        return interval::infinite();

    auto [a, b] = bounds();
    auto [c, d] = other.bounds();

    const auto D = ropd<op_div>;
    const auto U = ropu<op_div>;

    if (0 <= c) { // P
        if (0 < a) // P1
            return interval(D(a, d), U(b, c));

        if (0 == a) // P0
            return interval(0.0, U(b, c));

        // a < 0
        if (0 < b) // M
            return interval(D(a, c), U(b, c));

        if (0 == b) // N0
            return interval(D(a, c), 0.0);

        if (b < 0) // N1
            return interval(D(a, c), U(b, d));
    } else if (d <= 0) { // N
        if (0 < a) // P1
            return interval(D(b, d), U(a, c));

        if (0 == a) // P0
            return interval(D(b, d), 0.0);

        // a < 0
        if (0 < b) // M
            return interval(D(b, d), U(a, d));

        if (0 == b) // N0
            return interval(0.0, U(a, d));

        if (b < 0) // N1
            return interval(D(b, c), U(a, d));
    }

    return interval::empty();
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
