#ifndef PAPPUS_INTERVAL_HPP
#define PAPPUS_INTERVAL_HPP

#include <cassert>
#include <cmath>
#include <cfenv>
#include <iostream>
#include <ostream>
#include <utility>

namespace pappus {
class affine_interval {
public:
    affine_interval()
        : affine_interval(0, 0)
    {
    }
    explicit affine_interval(double l, double h)
        : affine_interval(check_bounds(l, h))
    {
    }

    affine_interval& operator=(const affine_interval& other)
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
        auto checked = check_bounds(lower, upper);
        lower_ = lower;
        upper_ = upper;
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
        double l, u;

        auto rounding_mode = std::fegetround();

        // lower
        std::fesetround(FE_DOWNWARD);
        l = lower_ * 0.5;

        // upper
        std::fesetround(FE_UPWARD);
        u = upper_ * 0.5;

        // restore the original rounding mode
        std::fesetround(rounding_mode);

        return l + u;
    }

    double radius() const
    {
        double m = mid();

        auto rounding_mode = std::fegetround();

        std::fesetround(FE_DOWNWARD);
        auto l = m - lower_;

        std::fesetround(FE_UPWARD);
        auto u = upper_ - m;

        std::fesetround(rounding_mode);

        return std::fmax(l, u);
    }

    affine_interval minimal_periodic(const affine_interval& interval) 
    {
        auto [a, b] = interval.bounds();

        auto t1 = std::floor(a / (2 * M_PI));
        auto t2 = b - a;

        a = a - (t1 * 2 * M_PI);
        b = a + t2;

        return affine_interval(a, b);
    }

    // printing
    friend std::ostream& operator<<(std::ostream& s,
        const affine_interval& interval)
    {
        s << "[" << interval.lower() << ", " << interval.upper() << "]";
        return s;
    }

    bool operator==(affine_interval const& other) const 
    {
        auto eps = std::numeric_limits<double>::epsilon();
        auto l = std::abs(lower_ - other.lower_);
        auto u = std::abs(upper_ - other.upper_);

        return l < eps && u < eps;
    }

private:
    double lower_;
    double upper_;

    static std::pair<double, double> check_bounds(double lo, double hi)
    {
        assert(lo <= hi);
        return std::make_pair(lo, hi);
    }

    explicit affine_interval(std::pair<double, double> interval)
        : lower_(interval.first)
        , upper_(interval.second)
    {
    }
};

}

#endif

