#ifndef PAPPUS_INTERVAL_HPP
#define PAPPUS_INTERVAL_HPP

#include <cassert>
#include <cmath>
#include <limits>
#include <ostream>
#include <type_traits>
#include <utility>

#include "fputil.hpp"

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
            lower_ = other.lower_;
            upper_ = other.upper_;
        }
        return *this;
    }

    double lower() const { return lower_; }
    void set_lower(double value) { lower_ = value; }

    double upper() const { return upper_; }
    void set_upper(double value) { upper_ = value; }

    void set_bounds(double lower, double upper)
    {
        auto [lo, up] = check_bounds(lower, upper);
        lower_ = lo;
        upper_ = up;
    };

    void set_bounds(std::pair<double, double> bounds)
    {
        set_bounds(bounds.first, bounds.second);
    }

    std::pair<double, double> bounds() const
    {
        return std::pair<double, double>(lower(), upper());
    }

    double mid() const;

    double radius() const;

    double diameter() const;

    bool contains(double v) const
    {
        return lower() <= v && v <= upper();
    }

    bool contains_strict(double v) const
    {
        return lower() < v && v < upper();
    }

    bool contains(interval const& other) const
    {
        return lower() <= other.lower() && other.upper() <= upper();
    }

    bool contains_strict(interval const& other) const
    {
        return lower() < other.lower() && other.upper() < upper();
    }

    // properties (loosely inspired from GAOL)
    bool is_finite() const
    {
        return std::isfinite(lower()) && std::isfinite(upper());
    }

    bool is_infinite() const
    {
        return std::isinf(lower()) && std::isinf(upper());
    }

    bool is_empty() const
    {
        return !(lower() <= upper()); // negation to handle NaNs
    }

    bool is_zero() const
    {
        return lower() == 0 && upper() == 0;
    }

    bool is_symmetric() const
    {
        return !is_empty() && (-lower() == upper());
    }

    bool is_positive() const
    {
        return is_empty() || (lower() >= 0.0);
    }

    bool is_strictly_positive() const
    {
        return is_empty() || (lower() > 0.0);
    }

    std::pair<interval, interval> split() const
    {
        auto left = interval { lower(), mid() };
        auto right = interval { std::nextafter(mid(), upper()), upper() };
        return std::make_pair(left, right);
    }

    // intersection
    interval operator&(interval const& other) const
    {
        if (is_empty() || other.is_empty())
            return interval::empty();

        if (upper() < other.lower() || other.upper() < lower())
            return interval::empty();

        return interval(std::fmax(lower(), other.lower()), std::fmin(upper(), other.upper()));
    }

    interval& operator&=(interval const& other)
    {
        auto tmp = *this & other;
        std::swap(*this, tmp);
        return *this;
    }

    // technically this is the hull (not union) 
    interval operator|(interval const& other) const
    {
        if (is_empty()) 
            return other;

        if (other.is_empty()) 
            return *this;

        return interval(std::fmin(lower(), other.lower()), std::fmax(upper(), other.upper()));
    }

    interval& operator|=(interval const& other)
    {
        auto tmp = *this | other;
        std::swap(*this, tmp);
        return *this;
    }

    // this assumes set semantics where two intervals
    // are equal if they have the same bounds
    bool operator==(interval const& other) const
    {
        return (is_empty() && other.is_empty()) ||
            (lower() == other.lower() && upper() == other.upper());
    }

    bool operator!=(interval const& other) const
    {
        return !(*this == other);
    }

    bool operator<(interval const& other) const
    {
        if (is_empty())
            return !other.is_empty();

        return upper() < other.lower();
    }

    bool operator<=(interval const& other) const
    {
        if (is_empty())
            return true;

        return upper() <= other.lower();
    }

    // arithmetic operators
    interval operator+(interval const& other) const;
    interval operator+() const;
    interval operator-(interval const& other) const;
    interval operator-() const;
    interval operator*(interval const& other) const;
    interval operator/(interval const& other) const;
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
    interval pow(interval const& other) const;

    interval& operator+=(interval const& other);
    interval& operator-=(interval const& other);
    interval& operator*=(interval const& other);
    interval& operator/=(interval const& other);

    interval operator+(double v) const;
    interval operator-(double v) const;
    interval operator*(double v) const;
    interval operator/(double v) const;
    interval pow(double v) const;

    // friends
    friend interval operator+(double v, interval const& i) { return i + v; }
    friend interval operator-(double v, interval const& i) { return -i + v; }
    friend interval operator*(double v, interval const& i) { return i * v; }
    friend interval operator/(double v, interval const& i) { return i.inv() * v; }
    friend interval pow(double v, interval const& i) { return interval(+0.0, -0.0); };

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
    friend std::ostream& operator<<(std::ostream& s, interval const& interval)
    {
        s << "[" << interval.lower() << ", " << interval.upper() << "]";
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

private:
    double lower_;
    double upper_;

    static std::pair<double, double> check_bounds(double lo, double hi)
    {
        // signed zero convention - see T. Hickey, Q. Ju, M. H. van Emden
        // Interval Arithmetic - From Principles to Implementation - Section 5.2
        // -0 used for right endpoints while +0 is used for left endpoints
        if (lo == 0) lo = +0.0;
        if (hi == 0) hi = -0.0;
        EXPECT(!(lo > hi));
        return std::pair<double, double>(lo, hi);
    }

    explicit interval(std::pair<double, double> interval)
        : lower_(interval.first)
        , upper_(interval.second)
    {
    }
};

} // namespace
#endif
