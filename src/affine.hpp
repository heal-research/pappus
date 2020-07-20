#ifndef PAPPUS_AF_HPP
#define PAPPUS_AF_HPP

#include "interval.hpp"

#include <algorithm>
#include <cmath>
#include <execution>
#include <vector>

#include <Eigen/Eigen>

namespace pappus {
template <typename T>
using array = Eigen::Array<T, Eigen::Dynamic, 1>;

class affine_form {
public:
    explicit affine_form(const affine_interval& interval)
        : center_(interval.mid())
        , radius_(0)
    {
    }

    explicit affine_form(double v)
        : center_(v)
        , radius_(0)
    {
    }

    explicit affine_form(double v, array<double> const& deviations, array<size_t> const& indices, size_t len)
        : center_(v)
        , radius_(deviations.abs().sum())
        , deviations_(deviations)
        , indices_(indices)
        , length_(len)
    {
    }

    // not explicit because we want to allow implicit copy
    affine_form(const affine_form& other)
        : center_(other.center_)
        , radius_(other.radius_)
        , deviations_(other.deviations_)
        , indices_(other.indices_)
        , length_(other.length_)
    {
    }

    double center() const { return center_; }
    double radius() const { return radius_; }

    double max() const { return center() + radius(); }
    double min() const { return center() - radius(); }

    double abs_max() const
    {
        return std::max(std::abs(min()), std::abs(max()));
    }

    double abs_min() const
    {
        auto min_ = min();
        auto max_ = max();

        auto abs_min_ = std::abs(min_);
        auto abs_max_ = std::abs(max_);

        return std::signbit(min_) == std::signbit(max_)
            ? std::min(abs_min_, abs_max_)
            : 0;
    }

    size_t size() const { return indices_.size(); }

    // comparison operators
    bool operator<(const affine_form& other)
    {
        return max() < other.max();
    }

    bool operator<=(const affine_form& other)
    {
        return max() <= other.max();
    }

    bool operator>(const affine_form& other)
    {
        return !(*this <= other);
    }

    bool operator>=(const affine_form& other)
    {
        return !(*this < other);
    }

    bool operator==(const affine_form& other)
    {
        return false; // to be implemented
    }

    // interval representation
    affine_interval to_inverval() const
    {
        return affine_interval(min(), max());
    }

    // arithmetic operations
    affine_form& operator+=(double v)
    {
        center_ += v;
        return *this;
    }

    affine_form& operator-=(double v)
    {
        center_ -= v;
        return *this;
    }

    affine_form& operator*=(double v)
    {
        deviations_ *= v;
        center_ *= v;
        radius_ *= std::abs(v);
        return *this;
    }

    affine_form& operator/=(double v)
    {
        return operator*=(1.0 / v);
    }

    affine_form operator+(const affine_form& other) const
    {
        if (length_ == 0 && other.length_ == 0) {
            return affine_form(center_ + other.center_);
        } 

        if (length_ == 0) {
            affine_form f(other);
            f += center_;
            return f;
        } 

        if (other.length_ == 0) {
            affine_form f(*this);
            f += other.center_;
            return f;
        }

        array<size_t> idx(length_ + other.length_);
        array<double> dev(length_ + other.length_);

        size_t i = 0, j = 0, k = 0;

        // do a set union of the indices
        // the indices_ arrays are assumed to be sorted
        // we fill the deviations array at the same time
        while(i < indices_.size() && j < other.indices_.size()) {
            auto a = indices_(i);
            auto b = other.indices_(j);
            
            auto d = a < b ? deviations_(i) : other.deviations_(j);
            auto v = std::min(a, b);

            if (idx.size() == 0 || idx(idx.size()-1) != v)
            {
                idx(k) = v;
                dev(k) = d; 

                ++k;
            }

            i += a <= b;
            j += b <= a;
        }

        return affine_form(center_ + other.center_, dev, idx, k);
    }

    affine_form& operator+=(const affine_form& other) 
    {
        auto tmp = *this + other;
        std::swap(tmp, *this);
        return *this;
    }

private:
    double center_; // central v
    double radius_; // affine radius

    array<double> deviations_; // vector of partial deviations
    array<size_t> indices_; // vector of indices

    size_t length_; // the actual number of indices still in use by this affine form

    void update_radius()
    {
        radius_ = deviations_.abs().sum();
    }
};
} // namespace pp

#endif

