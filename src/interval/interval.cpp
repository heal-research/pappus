#include "interval.hpp"
#include "fp/math.hpp"

#include <cmath>

namespace pappus {

using namespace pappus::fp;
using num = std::numeric_limits<double>;

double interval::mid() const
{
    if (is_empty())
        return num::quiet_NaN();

    if (is_infinite())
        return 0.0;

    if (std::isinf(inf()))
        return num::lowest();

    if (std::isinf(sup()))
        return num::max();

    return ropd<op_mul>(inf(), 0.5) + ropu<op_mul>(sup(), 0.5);
}

double interval::radius() const
{
    if (is_empty())
        return num::quiet_NaN();

    double m = mid();
    return std::fmax(ropd<op_sub>(m, inf()), ropu<op_sub>(sup(), m));
}

double interval::diameter() const
{
    if (is_empty())
        return num::quiet_NaN();

    return ropu<op_sub>(sup(), inf());
}

double interval::mig() const
{
    if (is_empty())
        return fp::nan;

    if (contains(0))
        return 0;

    return std::min(ropd<op_abs>(inf()), ropd<op_abs>(sup()));
}

double interval::mag() const
{
    if (is_empty())
        return fp::nan;

    return std::max(ropd<op_abs>(inf()), ropd<op_abs>(sup()));
}

// divide this interval into n segments and return the ith segment
interval interval::segment(std::size_t i, std::size_t n) const
{
    EXPECT(i < n);
    auto h = diameter() / n;
    auto a = ropd<op_add>(inf(), i * h);
    auto b = ropu<op_add>(inf(), (i + 1) * h);
    return interval(a, b);
}

interval interval::operator+() const
{
    return *this;
}

interval interval::operator+(double v) const
{
    return interval(ropd<op_add>(inf(), v), ropu<op_add>(sup(), v));
}

interval interval::operator+(interval const other) const
{
    auto [a, b] = bounds();
    auto [c, d] = other.bounds();
    return interval(ropd<op_add>(a, c), ropu<op_add>(b, d));
}

interval interval::operator-() const
{
    return interval(-sup(), -inf());
}

interval interval::operator-(double v) const
{
    return interval(ropd<op_sub>(inf(), v), ropu<op_sub>(sup(), v));
}

interval interval::operator-(interval const other) const
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
        ? interval(ropd<op_mul>(sup(), v), ropu<op_mul>(inf(), v))
        : interval(ropd<op_mul>(inf(), v), ropu<op_mul>(sup(), v));
}

interval interval::operator*(interval const other) const
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
        ? interval(ropd<op_div>(sup(), v), ropu<op_div>(inf(), v))
        : interval(ropd<op_div>(inf(), v), ropu<op_div>(sup(), v));
}

interval interval::operator/(interval const other) const
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

    return interval(ropd<op_div>(1, sup()), ropu<op_div>(1, inf()));
}

interval interval::exp() const
{
    // ensure that the result is strictly positive
    return interval(std::fmax(0.0, ropd<op_exp>(inf())), ropu<op_exp>(sup()));
}

interval interval::log() const
{
    // cannot take log of a negative interval
    if (sup() < 0)
        return interval::empty();

    auto i = interval(std::fmax(0.0, inf()), sup());
    return interval(ropd<op_log>(inf()), ropu<op_log>(sup()));
}

interval interval::sin() const
{
    if (is_empty() || is_zero())
        return *this;

    if (is_infinite() || diameter() > two_pi)
        return interval(-1, 1);

    auto [a, b] = bounds();
    int x = trig::get_quadrant(a);
    int y = trig::get_quadrant(b);

    const auto D = [](auto x) { return ropd<op_sin>(x); };
    const auto U = [](auto x) { return ropu<op_sin>(x); };

    if (x == y) { // same quadrant
        if (diameter() > pi) // wrap around
            return interval(-1, 1);

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

        if ((x == 1 || x == 2) && (y == 0 || y == 3))
            return interval(-1.0, std::fmax(U(a), U(b)));
    }

    return interval(-1, 1);
}

interval interval::cos() const
{
    if (is_empty())
        return *this;

    if (is_infinite() || diameter() > two_pi)
        return interval(-1, 1);

    auto [a, b] = bounds();
    auto x = trig::get_quadrant(a);
    auto y = trig::get_quadrant(b);

    const auto D = [](auto x) { return ropd<op_cos>(x); };
    const auto U = [](auto x) { return ropu<op_cos>(x); };

    if (x == y) { // same quadrant
        if (diameter() > pi) // wrap around
            return interval(-1, 1);

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

interval interval::tan() const
{
    if (is_empty())
        return *this;

    if (diameter() > pi)
        return interval::infinite();

    auto [a, b] = bounds();
    auto x = trig::get_quadrant(a);
    auto y = trig::get_quadrant(b);

    auto m = std::fmod(x, 2);
    auto n = std::fmod(y, 2);

    if (m <= n && x < y) {
        return interval::infinite();
    }

    return interval(ropd<op_tan>(a), ropu<op_tan>(b));
}

interval interval::asin() const
{
    auto t = *this & interval(-1, 1);
    if (t.is_empty())
        return interval::empty();
    return interval(ropd<op_asin>(t.inf()), ropu<op_asin>(t.sup()));
}

interval interval::acos() const
{
    auto t = *this & interval(-1, 1);
    if (t.is_empty())
        return interval::empty();
    return interval(ropd<op_acos>(t.sup()), ropu<op_acos>(t.inf()));
}

interval interval::atan() const
{
    if (is_empty())
        return interval::empty();
    return interval(ropd<op_atan>(inf()), ropu<op_atan>(sup()));
}

interval interval::sinh() const
{
    if (is_empty())
        return interval::empty();

    return interval(ropd<op_sinh>(inf()), ropu<op_sinh>(sup()));
}

interval interval::cosh() const
{
    if (is_empty())
        return interval::empty();

    return interval(ropd<op_cosh>(mig()), ropu<op_cosh>(mag()));
}

interval interval::tanh() const
{
    if (is_empty())
        return interval::empty();

    return interval(ropd<op_tanh>(inf()), ropu<op_tanh>(sup()));
}

interval interval::square() const
{
    if (is_zero() || is_empty())
        return *this;

    auto [a, b] = bounds();
    const auto D = [](auto x, auto y) { return ropd<op_mul>(x, y); };
    const auto U = [](auto x, auto y) { return ropu<op_mul>(x, y); };

    if (0 <= a)
        return interval(D(a, a), U(b, b));

    if (b <= 0)
        return interval(D(b, b), U(a, a));

    return interval(0.0, std::fmax(U(a, a), U(b, b)));
}

interval interval::pow(int p) const
{
    if (is_empty())
        return interval::empty();

    if (p == 0)
        return interval(1.0);

    if (p == 1)
        return *this;

    if (p < 0 && is_zero())
        return interval::empty();

    auto [a, b] = bounds();

    if (p % 2 == 0) { // even power
        if (p > 0) {
            if (a >= 0)
                return interval(ropd<op_pow>(a, p), ropu<op_pow>(b, p));

            if (b <= 0)
                return interval(ropd<op_pow>(b, p), ropu<op_pow>(a, p));

            return interval(ropd<op_pow>(mig(), p), ropu<op_pow>(mag(), p));
        } else {
            if (a >= 0)
                return interval(ropd<op_pow>(b, p), ropu<op_pow>(a, p));

            if (b <= 0)
                return interval(ropd<op_pow>(a, p), ropu<op_pow>(b, p));

            return interval(ropd<op_pow>(mag(), p), ropu<op_pow>(mig(), p));
        }
    } else { // odd power
        if (is_infinite())
            return interval::infinite();

        if (p > 0) {
            if (a == 0)
                return interval(0, ropu<op_pow>(b, p));

            if (b == 0)
                return interval(ropd<op_pow>(a, p), 0);

            return interval(ropd<op_pow>(a, p), ropu<op_pow>(b, p));
        } else {
            if (a == 0)
                return interval(ropd<op_pow>(b, p), fp::inf);

            if (b == 0)
                return interval(-fp::inf, ropu<op_pow>(a, p));

            if (contains(0.0))
                return interval::infinite();

            return interval(ropd<op_pow>(b, p), ropu<op_pow>(a, p));
        }
    }
}

interval& interval::operator+=(interval const other)
{
    auto tmp = *this + other;
    std::swap(*this, tmp);
    return *this;
}

interval& interval::operator-=(interval const other)
{
    auto tmp = *this - other;
    std::swap(*this, tmp);
    return *this;
}

interval& interval::operator*=(interval const other)
{
    auto tmp = *this * other;
    std::swap(*this, tmp);
    return *this;
}

interval& interval::operator/=(interval const other)
{
    auto tmp = *this / other;
    std::swap(*this, tmp);
    return *this;
}
}
