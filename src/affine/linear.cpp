#include "affine.hpp"

namespace pappus {
// for an affine form, the following linear operations are defined:
// - addition, subtraction between two affine forms
// - addition, subtraction between an affine form and a scalar
// - multiplication between an affine form and a scalar

// affine-affine
affine_form affine_form::operator+(affine_form const& other) const
{
    using op = binary_op<opcode::add>;
    if (auto res = affine_form::handle_special_cases<op>(*this, other); res.has_value()) {
        return res.value();
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
    for (size_t k = 0; k < dev.size(); ++k) {
        if (k == length_ || indices_[i] != idx[k]) {
            dev[k] = other.deviations_[j++];
            continue;
        }
        if (k == other.length_ || other.indices_[j] != idx[k]) {
            dev[k] = deviations_[i++];
            continue;
        }
        dev[k] = deviations_[i++] + other.deviations_[j++];
    }

    return affine_form(context_, center_ + other.center_, dev, idx, idx.size());
}

// not the most efficient approach but reduces code duplication
// and the probability for mistakes and bugs
affine_form& affine_form::operator+=(affine_form const& other)
{
    auto tmp = *this + other;
    swap(tmp);
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
    swap(tmp);
    return *this;
}

// affine-scalar
affine_form affine_form::operator+(double v) const
{
    affine_form f(*this);
    f += v;
    return f;
}

affine_form affine_form::operator-(double v) const
{
    affine_form f(*this);
    f -= v;
    return f;
}

affine_form affine_form::operator*(double v) const
{
    affine_form f(*this);
    f *= v;
    return f;
}

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
    view::as_array(deviations_) *= v;
    center_ *= v;
    radius_ *= std::fabs(v);
    return *this;
}

affine_form& affine_form::operator/=(double v)
{
    return operator*=(1.0 / v);
}

affine_form affine_form::operator-() const
{
    affine_form tmp(*this);
    tmp *= -1.0;
    return tmp;
}
} // namespace
