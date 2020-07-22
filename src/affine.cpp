#include "affine.hpp"

#include <algorithm>
#include <cmath>
#include <execution>
#include <functional>

namespace pappus {
bool affine_form::operator==(affine_form const& other) const
{
    if (length_ != other.length_) {
        return false;
    }

    auto not_equal = [](double x, double y) -> bool {
        auto a = std::fabs(x);
        auto b = std::fabs(y);
        auto c = std::fabs(x - y);
        if (!(a < 1 && b < 1)) {
            c /= (a + b);
        }
        return c > std::numeric_limits<double>::epsilon();
    };

    // no equivalence if the central value is not equal
    if (not_equal(center_, other.center_)) {
        return false;
    }

    for (size_t i = 0; i < length_; ++i) {
        if (indices_[i] != other.indices_[i]) {
            return false;
        }
        if (not_equal(deviations_[i], other.deviations_[i])) {
            return false;
        }
    }

    return true;
}

// assignment
affine_form& affine_form::operator=(double v)
{
    center_ = v;
    radius_ = 0;
    length_ = 0;
    return *this;
}

affine_form& affine_form::operator=(affine_form const& other)
{
    if (this != &other)
    {
        affine_form tmp(other);
        std::swap(tmp, *this);
    }
    return *this;
}

// affine_form arithmetic
affine_form& affine_form::operator+=(double v)
{
    center_ += v;
    return *this;
}

affine_form& affine_form::operator-=(double v)
{
    center_ -= v;
    return *this;
}

affine_form& affine_form::operator*=(double v)
{
    eigen_array(deviations_) *= v;
    center_ *= v;
    radius_ *= std::fabs(v);
    return *this;
}

affine_form& affine_form::operator/=(double v)
{
    return operator*=(1.0 / v);
}

affine_form affine_form::operator+(affine_form const& other) const
{
    if (length_ == 0 && other.length_ == 0) {
        return affine_form(context_, center_ + other.center_);
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

    std::vector<size_t> idx;
    idx.reserve(length_ + other.length_);

    std::set_union(indices_.begin(), indices_.end(), 
                   other.indices_.begin(), other.indices_.end(),
                   std::back_inserter(idx));

    std::vector<double> dev(idx.size());

    // fill the deviations array
    // reset indices
    size_t i = 0, j = 0;
    for(size_t k = 0; k < dev.size(); ++k) {
        if (k == length_ || indices_[i] != idx[k]) {
            dev[k] = other.deviations_[j++];
            continue;
        }
        if (k == other.length_ || other.indices_[j] != idx[k]) {
            dev[k] = deviations_[i++];
            continue;
        }
        dev[k] = deviations_[i] + other.deviations_[j];
        ++i;
        ++j;
    }

    return affine_form(context_, center_ + other.center_, dev, idx, idx.size());
}

// the following operators (+=, -, -=) reuse operator+
// this might be less efficient but avoids code duplication
// and makes the code much simpler and less bug-prone
affine_form& affine_form::operator+=(affine_form const& other)
{
    auto tmp = *this + other;
    std::swap(tmp, *this);
    return *this;
}

affine_form affine_form::operator-(affine_form const& other) const
{
    auto tmp(other);
    tmp *= -1.0;
    return *this + tmp;
}

affine_form& affine_form::operator-=(affine_form const& other)
{
    auto tmp = *this - other;
    std::swap(tmp, *this);
    return *this;
}

affine_form affine_form::operator-() const
{
    affine_form tmp(*this);
    tmp *= -1.0;
    return tmp;
}

affine_form affine_form::operator*(affine_form const& other) const
{
    if (length_ == 0 && other.length_ == 0) {
        return affine_form(context_, center_ * other.center_);
    }

    if (length_ == 0) {
        affine_form f(other);
        f *= center_;
        return f;
    }

    if (other.length_ == 0) {
        affine_form f(*this);
        f *= other.center_;
        return f;
    }

    std::vector<size_t> idx;
    idx.reserve(length_ + other.length_);

    std::set_union(indices_.begin(), indices_.end(), 
                   other.indices_.begin(), other.indices_.end(),
                   std::back_inserter(idx));

    std::vector<double> dev(idx.size());

    double common_center = 0;    // common term center
    double common_deviation = 0; // common term deviation

    auto c1 = center_;
    auto c2 = other.center_;

    auto r1 = radius_;
    auto r2 = other.radius_;

    size_t i = 0, j = 0;
    for (size_t k = 0; k < dev.size(); ++k) {
        auto x1 = indices_[i];
        auto x2 = other.indices_[j];

        auto d1 = deviations_[i];
        auto d2 = other.deviations_[j];

        if (k == length_ || indices_[i] != idx[k]) {
            dev[k] = c1 * d2;
            ++j;
            continue;
        }
        if (k == other.length_ || other.indices_[j] != idx[k]) {
            dev[k] = c2 * d1;
            ++i;
            continue;
        }
        dev[k] = c1 * d2 + c2 * d1;
        common_center += d1 * d2;
        common_deviation += std::fabs(d1 * d2);

        ++i;
        ++j;
    }

    common_center /= 2;
    common_deviation /= 2;

    auto delta = r1 * r2;

    // increment global index
    idx.push_back(context().increment_last());

    // TODO: research/implement SECANT approximation type
    dev.push_back(delta - common_deviation);

    return affine_form(context(), c1 * c2 + common_center, dev, idx, context().last_index()); 
}
} // namespace

