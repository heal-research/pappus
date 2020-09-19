#include "interval.hpp"
#include "fpops.hpp"

#include <cmath>

namespace pappus {

using namespace pappus::fp;

double interval::mid() const
{
        // lower
        auto lo = fp::ropd<fp::op_mul>(lower_, 0.5);
        // upper
        auto up = fp::ropu<fp::op_mul>(upper_, 0.5);
        return lo + up;
}

double interval::radius() const
{
        double m = mid();
        auto [a, b] = bounds();
        return std::fmax(fp::ropd<fp::op_sub>(m, a), fp::ropu<fp::op_sub>(b, m));
}

double interval::diameter() const
{
        auto [a, b] = bounds();
        return std::fmax(fp::ropd<fp::op_sub>(b, a), fp::ropu<fp::op_sub>(b, a));
}

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

    const auto D = [](auto x, auto y) { return ropd<op_mul>(x, y); };
    const auto U = [](auto x, auto y) { return ropu<op_mul>(x, y); };

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
            return interval(std::fmin(D(a, d), D(b, c)), std::fmax(U(a, c), U(b, d)));

        if (0 <= c) // P
            return interval(D(a, d), U(b, d));
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

    const auto D = [](auto x, auto y) { return ropd<op_div>(x, y); };
    const auto U = [](auto x, auto y) { return ropu<op_div>(x, y); };

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

interval interval::inv() const
{
    if (is_empty() || is_zero())
        return interval::empty();

    if (contains_strict(0.0))
        return interval::infinite();

    return interval(ropd<op_div>(1, upper()), ropu<op_div>(1, lower()));
}

interval interval::exp() const
{
    // ensure that the result is strictly positive
    return interval(std::fmax(0.0, ropd<op_exp>(lower())), ropu<op_exp>(upper()));
}

interval interval::log() const
{
    // cannot take log of a negative interval
    if (upper() < 0)
        return interval::empty();

    auto i = interval(std::fmax(0.0, lower()), upper());
    return interval(ropd<op_log>(lower()), ropu<op_log>(upper()));
}

interval interval::sin() const
{
    if (is_empty() || is_zero())
        return *this;

    if (is_infinite() || radius() > two_pi)
        return interval(-1, 1);

    auto [a, b] = bounds();
    int x = trig::get_quadrant(a);
    int y = trig::get_quadrant(b);

    const auto D = [](auto x) { return ropd<op_sin>(x); };
    const auto U = [](auto x) { return ropu<op_sin>(x); };

    if (x == y) { // same quadrant
        if (diameter() > fp::pi) // wrap around
            return interval(-1.0, 1.0);

        if (x == 0 || x == 3) // Q1,Q4: increasing
            return interval(D(a), U(b));

        if (x == 1 || x == 2) // Q2,Q3: decreasing
            return interval(D(b), U(a));
    } else {
        if (x == 3 && y == 0)
            return interval(D(a), U(b));

       if (x == 1 && y == 2)
            return interval(D(b), U(a));

       if ((x == 0 || x == 3) && (y == 1 || y == 2))
           return interval(std::fmin(D(a), D(b)), 1.0);

       if ((x ==1 || x == 2) && (y == 0 || y == 3))
           return interval(-1.0, std::fmax(U(a), U(b)));
    }

    return interval(-1, 1);
}

interval interval::cos() const
{
    if (is_empty() || is_zero())
        return *this;

    if (is_infinite() || radius() > two_pi)
        return interval(-1, 1);

    auto [a, b] = bounds();
    auto x = trig::get_quadrant(a);
    auto y = trig::get_quadrant(b);

    const auto D = [](auto x) { return ropd<op_cos>(x); };
    const auto U = [](auto x) { return ropu<op_cos>(x); };

    if (x == y) { // same quadrant
        if (diameter() > fp::pi) // wrap around
            return interval(-1.0, 1.0);

        if (x == 0 || x == 1)
            return interval(D(b), U(a));

        if (x == 2 || x == 3)
            return interval(D(a), U(b));
    } else {
        if (x == 2 && y == 3)
            return interval(D(a), U(b));

        if (x == 0 && y == 1)
            return interval(D(b), U(a));

        if ((x == 2 || x == 3) && (y == 0 || y == 1))
            return interval(std::fmin(D(a), D(b)), 1.0);

        if ((x == 0 || x == 1) && (y == 2 || y == 3))
            return interval(-1.0, std::fmax(U(a), U(b)));
    }

    return interval(-1, 1);
}

interval interval::square() const
{
    if (is_zero() || is_empty())
        return *this;

    auto [a, b] = bounds();
    const auto D = [](auto x, auto y) { return ropd<op_mul>(x, y); };
    const auto U = [](auto x, auto y) { return ropu<op_mul>(x, y); };

    if (0 <= a)
        return interval(D(a,a), U(b,b));

    if (b <= 0)
        return interval(D(b,b), U(a,a));

    EXPECT(contains_strict(0.0));
    return interval(0.0, std::fmax(U(a,a), U(b,b)));
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

interval& interval::operator*=(interval const& other)
{
    auto tmp = *this * other;
    std::swap(*this, tmp);
    return *this;
}

interval& interval::operator/=(interval const& other)
{
    auto tmp = *this / other;
    std::swap(*this, tmp);
    return *this;
}
}
