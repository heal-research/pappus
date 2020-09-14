#include "affine.hpp"
#include "context.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <optional>

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

affine_form affine_form::operator-(affine_form const& other) const
{
    auto tmp(other);
    tmp *= -1.0;
    return *this + tmp;
}

affine_form affine_form::operator*(affine_form const& other) const
{
    using op = binary_op<opcode::multiply>;
    if (auto res = handle_special_cases<op>(*this, other); res.has_value()) {
        return res.value();
    }

    std::vector<size_t> idx;
    idx.reserve(length_ + other.length_);

    std::set_union(indices_.begin(), indices_.end(),
        other.indices_.begin(), other.indices_.end(),
        std::back_inserter(idx));

    std::vector<double> dev(idx.size());

    double common_term_center = 0; // common term center
    double common_term_deviation = 0; // common term deviation

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

        auto d1 = deviations_[i];
        auto d2 = other.deviations_[j];

        dev[k] = c1 * d2 + c2 * d1;
        common_term_center += d1 * d2;
        common_term_deviation += std::fabs(d1 * d2);

        ++i;
        ++j;
    }

    common_term_center /= 2;
    common_term_deviation /= 2;

    auto delta = r1 * r2;

    // increment global index
    idx.push_back(context().increment_last());

    if (context().approximation_mode() == approximation_mode::SECANT) {
        delta -= common_term_deviation;
        dev.push_back(0.0);
        auto r = view::as_array(dev).abs().sum();
        double fac = r < limits::eps ? 1.0 : 1.0 + delta / r;
        view::as_array(dev) *= fac;
    } else {
        dev.push_back(delta - common_term_deviation);
    }
    return affine_form(context(), c1 * c2 + common_term_center, dev, idx, idx.size());
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

    auto fa = 1 / a; 
    auto fb = 1 / b;

    double alpha = 0;
    double delta = 0;
    double dzeta = 0;

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
        double ya, yb;
        // Derivative of 1/x is -1/x*x
        if (a > 0.0) {
            alpha = -fb / b;
            // ya = fa - alpha*a;
            // yb = fb - alpha*b = 2.0*fb;
            ya = fa - alpha * a;
            yb = 2.0 * fb;
        } else {
            alpha = -fa / a;
            // ya = fa - alpha*a = 2.0*fa;
            // yb = fb - alpha*b;
            ya = 2.0 * fa;
            yb = fb - alpha * b;
        }

        delta = 0.5 * (ya - yb);
        dzeta = 0.5 * (ya + yb);
        break;
    }
    case approximation_mode::SECANT: {
        alpha = r > limits::minrad ? (fb - fa ) / (b - a) 
                                   : -fa * fb;
        dzeta = fa - alpha * a;
        delta = 0;
        break;
    }
    };

    std::vector<size_t> idx(length_ + 1);
    std::vector<double> dev(length_ + 1);

    // zi = alpha * xi
    view::as_array(idx).segment(0, length_) = view::as_array(indices_).segment(0, length_);
    view::as_array(dev).segment(0, length_) = view::as_array(deviations_).segment(0, length_) * alpha;

    // compute the error in a new deviation symbol zk = delta
    idx[length_] = context().increment_last();
    dev[length_] = delta;

    return affine_form(context(), alpha * c + dzeta, std::move(dev), std::move(idx), idx.size());
}

affine_form affine_form::operator^(int exponent) const
{
    if (length_ == 0)
        return affine_form(context(), std::pow(center(), exponent));

    if (exponent == 0)
        return affine_form(context(), 1.0);

    if (exponent == 1)
        return *this;

    if (exponent == -1)
        return this->inv();

    auto c = center();
    auto r = radius();

    auto a = c - r;
    auto b = c + r;

    auto fa = std::pow(a, exponent);
    auto fb = std::pow(b, exponent);

    double alpha = 0.0;

    double delta;
    double dzeta;

    switch (context().approximation_mode()) {
    case approximation_mode::CHEBYSHEV: {
        alpha = r > limits::minrad ? (fb - fa) / (b - a)
                                   : exponent * fa / a;

        // we have two points having the slope alpha
        auto e = std::fabs(alpha / exponent);
        auto x1 = -std::pow(e, 1.0 / (exponent - 1.0));
        auto x2 = -x1;

        auto y1 = x1 < a ? fa - alpha * a
                         : std::pow(x1, exponent) - alpha * x1;

        auto y2 = x2 > b ? fb - alpha * b
                         : std::pow(x2, exponent) - alpha * x2;

        delta = (y1 - y2) * 0.5;
        dzeta = (y1 + y2) * 0.5;
        break;
    }
    case approximation_mode::MINRANGE: {
        // special case: 0.0 in [a,b] : alpha = f'(0.0)
        // exp > 0 and exp even and [a,b] > 0 : MINRANGE: alpha = f'(a), y0_b > y0_a : e1 = 1 e2 = 1 e3 = 1
        // exp > 0 and exp odd  and [a,b] > 0 : MINRANGE: alpha = f'(a), y0_b > y0_a : e1 = 1 e2 = 0 e3 = 1
        // exp > 0 and exp even and [a,b] < 0 : MINRANGE: alpha = f'(b), y0_a > y0_b : e1 = 1 e2 = 1 e3 = 0
        // exp > 0 and exp odd  and [a,b] < 0 : MINRANGE: alpha = f'(b), y0_b > y0_a : e1 = 1 e2 = 0 e3 = 0
        // exp < 0 and exp even and [a,b] > 0 : MINRANGE: alpha = f'(b), y0_a > y0_b : e1 = 0 e2 = 1 e3 = 1
        // exp < 0 and exp odd  and [a,b] > 0 : MINRANGE: alpha = f'(b), y0_a > y0_b : e1 = 0 e2 = 0 e3 = 1
        // exp < 0 and exp even and [a,b] < 0 : MINRANGE: alpha = f'(a), y0_b > y0_a : e1 = 0 e2 = 1 e3 = 0
        // exp < 0 and exp odd  and [a,b] < 0 : MINRANGE: alpha = f'(a), y0_a > y0_b : e1 = 0 e2 = 0 e3 = 0
        if (a * b < 0) {
            alpha = 0.0;
            if (exponent % 2 == 0) {
                delta = 0.5 * std::max(fa, fb);
                dzeta = delta;
            } else {
                delta = 0.5 * (fb - fa);
                dzeta = 0.5 * (fb + fa);
            }
        } else {
            alpha = std::signbit(a) == std::signbit(exponent) ? exponent * fa / a
                                                              : exponent * fb / b;

            auto ya = fa - alpha * a;
            auto yb = fb - alpha * b;

            delta = 0.5 * std::fabs(ya - yb);
            dzeta = 0.5 * (fa + fb);
        }
        break;
    }
    case approximation_mode::SECANT: {
        alpha = r > limits::minrad ? (fb - fa) / (b - a)
                                   : exponent * fa / a;

        delta = 0.0;
        dzeta = fa - alpha * a;
        break;
    }
    } // switch

    std::vector<size_t> idx(length_ + 1);
    std::vector<double> dev(length_ + 1);

    view::as_array(idx).segment(0, length_) = view::as_array(indices_).segment(0, length_);
    view::as_array(dev).segment(0, length_) = view::as_array(deviations_).segment(0, length_) * alpha;

    idx[length_] = context().increment_last();
    dev[length_] = delta;

    return affine_form(context(), alpha * c + dzeta, std::move(dev), std::move(idx), idx.size());
}

} // namespace
