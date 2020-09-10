#ifndef PAPPUS_INTERVAL_HPP
#define PAPPUS_INTERVAL_HPP

#ifndef M_PI
#define M_PI 3.14159265358979323846 // pi
#endif

#include <cassert>
#include <cfenv>
#include <cmath>
#include <iostream>
#include <limits>
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
        std::fesetround(FE_DOWNWARD);
        auto lo = lower_ * 0.5;

        // upper
        std::fesetround(FE_UPWARD);
        auto up = upper_ * 0.5;

        // restore the original rounding mode
        std::fesetround(rounding_mode);

        return lo + up;
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
        auto l = std::fabs(lower_ - other.lower_);
        auto u = std::fabs(upper_ - other.upper_);

        return l < eps && u < eps;
    }

private:
    double lower_;
    double upper_;

    static std::pair<double, double> check_bounds(double lo, double hi)
    {
        assert(!std::isnan(lo));
        assert(!std::isnan(hi));
        assert(lo <= hi);
        return std::pair<double, double>(lo, hi);
    }

    explicit affine_interval(std::pair<double, double> interval)
        : lower_(interval.first)
        , upper_(interval.second)
    {
    }
};

}

#endif
