#ifndef PAPPUS_INTERVAL_HPP
#define PAPPUS_INTERVAL_HPP

#include <cassert>
#include <cmath>
#include <limits>
#include <ostream>
#include <type_traits>
#include <utility>

#include "fp/util.hpp"
#include "util/contracts.hpp"
#include "iterator.hpp"

namespace pappus {

enum interval_class { M, Z, P, P0, P1, N, N0, N1, NONE };

class interval {
public:
    interval()
        : interval(0, 0)
    {
    }

    explicit interval(double v) : interval(v, v)
    {
    }

    explicit interval(double l, double h)
        : interval(check_bounds(l, h))
    {
    }

    explicit interval(std::string const& s)
        : interval(fp::from_string<FE_DOWNWARD>(s), fp::from_string<FE_UPWARD>(s))
    {
    }

    explicit interval(std::string const& l, std::string const& u)
        : interval(fp::from_string<FE_DOWNWARD>(l), fp::from_string<FE_UPWARD>(u))
    {
    }

    interval& operator=(const interval& other)
    {
        if (this != &other) {
            inf_ = other.inf_;
            sup_ = other.sup_;
        }
        return *this;
    }

    double inf() const { return inf_; }
    void set_lower(double value) { inf_ = value; }

    double sup() const { return sup_; }
    void set_upper(double value) { sup_ = value; }

    void set_bounds(double lower, double upper)
    {
        auto [lo, up] = check_bounds(lower, upper);
        inf_ = lo;
        sup_ = up;
    };

    void set_bounds(std::pair<double, double> bounds)
    {
        set_bounds(bounds.first, bounds.second);
    }

    std::pair<double, double> bounds() const
    {
        return std::pair<double, double>(inf(), sup());
    }

    double mid() const;
    double radius() const;
    double diameter() const;
    double mig() const;
    double mag() const;

    bool contains(double v) const
    {
        return inf() <= v && v <= sup();
    }

    bool contains_strict(double v) const
    {
        return inf() < v && v < sup();
    }

    bool contains(interval const other) const
    {
        return inf() <= other.inf() && other.sup() <= sup();
    }

    bool contains_strict(interval const other) const
    {
        return inf() < other.inf() && other.sup() < sup();
    }

    // properties (loosely inspired from GAOL)
    bool is_finite() const
    {
        return std::isfinite(inf()) && std::isfinite(sup());
    }

    bool is_infinite() const
    {
        return std::isinf(inf()) && std::isinf(sup());
    }

    bool is_empty() const
    {
        return !(inf() <= sup()); // negation to handle NaNs
    }

    bool is_zero() const
    {
        return inf() == 0 && sup() == 0;
    }

    bool is_symmetric() const
    {
        return !is_empty() && (-inf() == sup());
    }

    bool is_positive() const
    {
        return is_empty() || (inf() >= 0.0);
    }

    bool is_strictly_positive() const
    {
        return is_empty() || (inf() > 0.0);
    }

    std::pair<interval, interval> split() const
    {
        auto left = interval { inf(), mid() };
        auto right = interval { std::nextafter(mid(), sup()), sup() };
        return std::make_pair(left, right);
    }

    subdivision split(size_t n) const
    {
        return subdivision(*this, n);
    }

    // intersection
    interval operator&(interval const other) const
    {
        if (is_empty() || other.is_empty())
            return interval::empty();

        if (sup() < other.inf() || other.sup() < inf())
            return interval::empty();

        return interval(std::fmax(inf(), other.inf()), std::fmin(sup(), other.sup()));
    }

    interval& operator&=(interval const other)
    {
        auto tmp = *this & other;
        std::swap(*this, tmp);
        return *this;
    }

    // technically this is the hull (not union) 
    interval operator|(interval const other) const
    {
        if (is_empty()) 
            return other;

        if (other.is_empty()) 
            return *this;

        return interval(std::fmin(inf(), other.inf()), std::fmax(sup(), other.sup()));
    }

    interval& operator|=(interval const other)
    {
        auto tmp = *this | other;
        std::swap(*this, tmp);
        return *this;
    }

    // this assumes set semantics where two intervals
    // are equal if they have the same bounds
    bool operator==(interval const other) const
    {
        return (is_empty() && other.is_empty()) ||
            (inf() == other.inf() && sup() == other.sup());
    }

    bool operator!=(interval const other) const
    {
        return !(*this == other);
    }

    bool operator<(interval const other) const
    {
        if (is_empty())
            return !other.is_empty();

        return sup() < other.inf();
    }

    bool operator<=(interval const other) const
    {
        if (is_empty())
            return true;

        return sup() <= other.inf();
    }

    // arithmetic operators
    interval operator+(interval const other) const;
    interval operator+() const;
    interval operator-(interval const other) const;
    interval operator-() const;
    interval operator*(interval const other) const;
    interval operator/(interval const other) const;
    interval inv() const;
    interval exp() const;
    interval log() const;
    interval sin() const;
    interval cos() const;
    interval tan() const;
    interval asin() const;
    interval acos() const;
    interval atan() const;
    interval square() const;
    interval pow(interval const other) const;
    interval pow(int) const;
    interval pow(double) const;
    interval sinh() const;
    interval cosh() const;
    interval tanh() const;

    interval& operator+=(interval const);
    interval& operator-=(interval const);
    interval& operator*=(interval const);
    interval& operator/=(interval const);

    interval operator+(double) const;
    interval operator-(double) const;
    interval operator*(double) const;
    interval operator/(double) const;

    // friends
    friend interval operator+(double v, interval const i) { return i + v; }
    friend interval operator-(double v, interval const i) { return -i + v; }
    friend interval operator*(double v, interval const i) { return i * v; }
    friend interval operator/(double v, interval const i) { return i.inv() * v; }
    friend interval pow(double v, interval const i) { return interval(+0.0, -0.0); };

    template <interval_class C>
    bool is() const
    {
        static_assert(C >= interval_class::M && C <= interval_class::N1, "Unknown interval class.");
        auto [a, b] = bounds();
        switch (C) {
        case M:
            return a < 0 && b > 0;
        case Z:
            return a == 0 && b == 0;
        case P:
            return a >= 0 && b > 0;
        case P0:
            return a = 0 && b > 0;
        case P1:
            return a > 0 && b > 0;
        case N:
            return a < 0 && b < 0;
        case N0:
            return a < 0 && b == 0;
        case N1:
            return a < 0 && b < 0;
        default:
            return false;
        }
    }

    // printing
    friend std::ostream& operator<<(std::ostream& s, interval const interval)
    {
        s << "[" << interval.inf() << ", " << interval.sup() << "]";
        return s;
    }

    // constants
    static interval empty()    { return interval(+fp::nan, -fp::nan); }
    static interval zero()     { return interval(+0.0, -0.0); }
    static interval infinite() { return interval(-fp::inf, +fp::inf); }

    static interval force_interval(double a, double b)
    {
        if (a > b) std::swap(a, b);
        return interval(a, b);
    }

    interval segment(std::size_t i, std::size_t n) const;

private:
    double inf_;
    double sup_;

    static std::pair<double, double> check_bounds(double lo, double hi)
    {
        // signed zero convention - see T. Hickey, Q. Ju, M. H. van Emden
        // Interval Arithmetic - From Principles to Implementation - Section 5.2
        // -0 used for right endpoints while +0 is used for left endpoints
        EXPECT(!(lo > hi));
        if (lo == 0) lo = +0.0;
        if (hi == 0) hi = -0.0;
        return std::pair<double, double>(lo, hi);
    }

    explicit interval(std::pair<double, double> interval)
        : inf_(interval.first)
        , sup_(interval.second)
    {
    }
};


} // namespace
#endif
