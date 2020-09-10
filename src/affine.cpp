#include "affine.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <optional>

namespace pappus {

constexpr double EPS = std::numeric_limits<double>::epsilon();
constexpr double AAF_MINRAD = 1e-10;

namespace detail {
    enum opcode { add,
        multiply };

    template <opcode x = opcode::add>
    struct binary_op {
        template <typename T, typename U>
        auto operator()(T t, U u) { return t + u; }
    };

    template <>
    struct binary_op<opcode::multiply> {
        template <typename T, typename U>
        auto operator()(T t, U u) { return t * u; }
    };

    template <typename OP>
    std::optional<affine_form> handle_special_cases(affine_form const& lhs, affine_form const& rhs)
    {
        if (lhs.length() == 0 && rhs.length() == 0) {
            return affine_form(lhs.context(), OP()(lhs.center(), rhs.center()));
        }

        if (lhs.length() == 0) {
            return OP()(rhs, lhs.center());
        }

        if (rhs.length() == 0) {
            return OP()(lhs, rhs.center());
        }
        return std::nullopt;
    }
}

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
        return c > EPS;
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

// affine_form arithmetic
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
    eigen_array(deviations_) *= v;
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

affine_form affine_form::operator+(affine_form const& other) const
{
    using op = detail::binary_op<detail::opcode::add>;
    if (auto res = detail::handle_special_cases<op>(*this, other); res.has_value()) {
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

affine_form affine_form::operator-(affine_form const& other) const
{
    auto tmp(other);
    tmp *= -1.0;
    return *this + tmp;
}

affine_form affine_form::operator*(affine_form const& other) const
{
    using op = detail::binary_op<detail::opcode::multiply>;
    if (auto res = detail::handle_special_cases<op>(*this, other); res.has_value()) {
        return res.value();
    }

    std::vector<size_t> idx;
    idx.reserve(length_ + other.length_);

    std::set_union(indices_.begin(), indices_.end(),
        other.indices_.begin(), other.indices_.end(),
        std::back_inserter(idx));

    std::vector<double> dev(idx.size());

    double common_center = 0; // common term center
    double common_deviation = 0; // common term deviation

    auto c1 = center_;
    auto c2 = other.center_;

    auto r1 = radius_;
    auto r2 = other.radius_;

    size_t i = 0, j = 0;
    for (size_t k = 0; k < dev.size(); ++k) {
        if (k == length_ || indices_[i] != idx[k]) {
            dev[k] = c1 * other.deviations_[j];
            ++j;
            continue;
        }
        if (k == other.length_ || other.indices_[j] != idx[k]) {
            dev[k] = c2 * deviations_[i];
            ++i;
            continue;
        }

        auto x1 = indices_[i];
        auto x2 = other.indices_[j];

        auto d1 = deviations_[i];
        auto d2 = other.deviations_[j];

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

    return affine_form(context(), c1 * c2 + common_center, dev, idx, idx.size());
}

affine_form affine_form::operator/(affine_form const& other) const
{
    if (this == &other) {
        return affine_form(context(), 1.0);
    }
    return *this * other.inv();
}

// not the most efficient approach but reduces code duplication
// and the probability for mistakes and bugs
affine_form& affine_form::operator+=(affine_form const& other)
{
    auto tmp = *this + other;
    swap(tmp);
    return *this;
}

affine_form& affine_form::operator-=(affine_form const& other)
{
    auto tmp = *this - other;
    swap(tmp);
    return *this;
}

affine_form& affine_form::operator*=(affine_form const& other)
{
    auto tmp = *this * other;
    swap(tmp);
    return *this;
}

// other operations
affine_form affine_form::inv() const
{
    if (length() == 0) {
        return affine_form(context(), 1.0 / center_);
    }

    auto c = center();
    auto r = radius();
    auto a = c - r;
    auto b = c + r;

    // protect against division by zero
    assert(a > EPS && b > EPS);

    auto fa = 1 / a;
    auto fb = 1 / b;

    double alpha, delta, dzeta;

    switch (context().approximation_mode()) {
    case approximation_mode::CHEBYSHEV: {
        auto u = std::sqrt(a * b);
        alpha = -fa * fb;

        if (a > 0) {
            delta = +0.5 * (fa + fb - 2.0 / u);
            dzeta = fa + fb - delta;
        } else {
            delta = -0.5 * (fa + fb + 2.0 / u);
            dzeta = fa + fb + delta;
        }
        break;
    }
    case approximation_mode::MINRANGE: {
        double y_a, y_b;
        // Derivative of 1/x is -1/x*x
        if (a > 0.0) {
            alpha = -fb / b;
            // y_a = fa - alpha*a;
            // y_b = fb - alpha*b = 2.0*fb;
            y_a = fa - alpha * a;
            y_b = 2.0 * fb;
        } else {
            alpha = -fa / a;
            // y_a = fa - alpha*a = 2.0*fa;
            // y_b = fb - alpha*b;
            y_a = 2.0 * fa;
            y_b = fb - alpha * b;
        }

        delta = 0.5 * (y_a - y_b);
        dzeta = 0.5 * (y_a + y_b);
        break;
    }
    case approximation_mode::SECANT: {
        if (r > AAF_MINRAD) {
            alpha = (fb - fa) / (b - a);
        } else {
            alpha = -fa * fb;
        }
        dzeta = fa - alpha * a;
        delta = 0;
        break;
    }
    };

    std::vector<size_t> idx(length_ + 1);
    std::vector<double> dev(length_ + 1);

    // zi = alpha * xi
    eigen_array(idx).segment(0, length_) = eigen_array(indices_).segment(0, length_);
    eigen_array(dev).segment(0, length_) = eigen_array(deviations_).segment(0, length_) * alpha;

    // compute the error in a new deviation symbol zk = delta
    idx[length_] = context().increment_last();
    dev[length_] = delta;

    return affine_form(context(), alpha * c + dzeta, std::move(dev), std::move(idx), idx.size());
}

} // namespace
