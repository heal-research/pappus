#ifndef PAPPUS_INTERVAL_HPP
#define PAPPUS_INTERVAL_HPP

#include <cassert>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
#include <ostream>
#include <sstream>
#include <type_traits>
#include <utility>

#define EXPECT(cond) \
    if(!(cond)) \
    { \
        std::cout << "precondition " << #cond << " failed at " << __FILE__ << ": " << __LINE__ << "\n"; \
        std::terminate(); \
    } 

#define ENSURE(cond) \
    if(!(cond)) \
    { \
        std::cout << "postcondition " << #cond << " failed at " << __FILE__ << ": " << __LINE__ << "\n"; \
        std::terminate(); \
    }

#ifndef M_PI
#define M_PI 3.14159265358979323846 // pi
#endif

namespace pappus {

namespace {
    template <int ROUND_MODE>
    double from_string(std::string const& s)
    {
        static_assert(ROUND_MODE == FE_UPWARD || ROUND_MODE == FE_DOWNWARD);
        auto rounding_mode = std::fegetround();
        std::fesetround(ROUND_MODE);
        std::istringstream is(s);
        double v;
        is >> v;
        std::fesetround(rounding_mode);
        return v;
    }

    // rounded op
    template <typename OP, int ROUND_MODE>
    double rop(double a, double b)
    {
        static_assert(ROUND_MODE == FE_UPWARD || ROUND_MODE == FE_DOWNWARD);
        static_assert(std::is_invocable_r<double, OP, double, double>::value);
        auto rnd = std::fegetround();
        std::fesetround(ROUND_MODE);
        auto c = OP()(a, b);
        std::fesetround(rnd);
        // necessary to abide to the arithmetic rules underlying IA, details in:
        // Complete Interval Arithmetic and its Implementation on the Computer
        // Ulrich W. Kulisch
        if (std::isnan(c) || c == -0.0)
            return 0.0;

        return c;
    }

    // rounded op, downwards
    template <typename OP>
    const auto ropd = rop<OP, FE_DOWNWARD>;

    // rounded op, upwards
    template <typename OP>
    const auto ropu = rop<OP, FE_UPWARD>;

    using op_add = std::plus<double>;
    using op_sub = std::minus<double>;
    using op_mul = std::multiplies<double>;
    using op_div = std::divides<double>;
}

enum interval_class { M, Z, P, P0, P1, N, N0, N1, NONE };

class interval {
public:
    interval()
        : interval(0, 0)
    {
    }
    explicit interval(double l, double h)
        : interval(check_bounds(l, h))
    {
    }

    explicit interval(std::string const& s)
        : interval(from_string<FE_DOWNWARD>(s), from_string<FE_UPWARD>(s))
    {
    }

    explicit interval(std::string const& l, std::string const& u)
        : interval(from_string<FE_DOWNWARD>(l), from_string<FE_UPWARD>(u))
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
        return std::pair<double, double>(lower_, upper_);
    }

    double mid() const
    {
        auto rounding_mode = std::fegetround();

        // lower
        auto lo = ropd<op_mul>(lower_, 0.5);
        // upper
        auto up = ropu<op_mul>(upper_, 0.5);
        return lo + up;
    }

    double radius() const
    {
        double m = mid();
        auto [a, b] = bounds();
        return std::fmax(ropd<op_sub>(m, a), ropu<op_sub>(b, m));
    }

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
    bool isfinite() const
    {
        return std::isfinite(lower()) && std::isfinite(upper());
    }

    bool isinfinite() const
    {
        return std::isinf(lower()) && std::isinf(upper());
    }

    bool isempty() const
    {
        return !(lower() <= upper()); // negation to handle NaNs
    }

    bool is_zero() const
    {
        return is<interval_class::Z>();
    }

    bool issymmetric() const
    {
        return !isempty() && (-lower() == upper());
    }

    bool is_positive() const
    {
        return isempty() || (lower() >= 0.0);
    }

    bool is_strictly_positive() const
    {
        return isempty() || (lower() > 0.0);
    }

    // intersection
    interval operator&(interval const& other) const
    {
        auto lo = std::max(lower(), other.lower());
        auto hi = std::min(upper(), other.upper());
        return interval(lo, hi);
    }

    // union
    interval operator|(interval const& other) const
    {
        auto lo = std::min(lower(), other.lower());
        auto hi = std::max(upper(), other.upper());
        return interval(lo, hi);
    }

    // printing
    friend std::ostream& operator<<(std::ostream& s, interval const& interval)
    {
        s << "[" << interval.lower() << ", " << interval.upper() << "]";
        return s;
    }

    bool operator==(interval const& other) const
    {
        return (isempty() && other.isempty()) || (lower() == other.lower() && upper() == other.upper());
    }

    std::pair<interval, interval> split() const
    {
        auto left = interval { lower(), mid() };
        auto right = interval { std::nextafter(mid(), upper()), upper() };
        return std::make_pair(left, right);
    }

    // arithmetic operators
    interval operator+(interval const& other) const;
    interval operator-(interval const& other) const;
    interval operator-() const;
    interval operator*(interval const& other) const;
    interval operator/(interval const& other) const;

    interval& operator+=(interval const& other);
    interval& operator-=(interval const& other);
    interval& operator*=(interval const& other);
    interval& operator/=(interval const& other);

    interval operator+(double v) const;
    interval operator-(double v) const;
    interval operator*(double v) const;
    interval operator/(double v) const;

    interval inv() const;

    template <interval_class C>
    constexpr bool is() const
    {
        static_assert(C >= interval_class::M && C <= interval_class::N1);
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
        return false;
    }

private:
    double lower_;
    double upper_;

    static std::pair<double, double> check_bounds(double lo, double hi)
    {
        EXPECT(!(lo > hi));
        return std::pair<double, double>(lo, hi);
    }

    explicit interval(std::pair<double, double> interval)
        : lower_(interval.first)
        , upper_(interval.second)
    {
    }
};

}

#endif
