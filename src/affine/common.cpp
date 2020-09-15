#include "affine.hpp"

namespace pappus {
// common operations:
// - equality comparison
// - assignment
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
        return c > limits::eps;
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
    if (this != &other) {
        affine_form tmp(other.context(), other.center(), other.deviations_, other.indices_, other.length_);
    }
    return *this;
}

affine_form& affine_form::operator=(affine_form other)
{
    swap(other);
    return *this;
}
} // namespace
