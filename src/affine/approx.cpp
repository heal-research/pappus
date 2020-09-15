#include "affine.hpp"

namespace pappus {
// non-affine operations approximate the original function
// and introduce an additional error term (epsilon)
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
        auto f1 = i < indices_.size() && idx[k] == indices_[i];
        auto f2 = j < other.indices_.size() && idx[k] == other.indices_[j];

        if (f1 && f2) {
            auto d1 = deviations_[i];
            auto d2 = other.deviations_[j];

            dev[k] = c1 * d2 + c2 * d1;
            common_term_center += d1 * d2;
            common_term_deviation += std::fabs(d1 * d2);

            ++i;
            ++j;
        } else if (f1) {
            dev[k] = c2 * deviations_[i];
            ++i;
        } else if (f2) {
            dev[k] = c1 * other.deviations_[j];
            ++j;
        }
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

affine_form& affine_form::operator*=(affine_form const& other)
{
    auto tmp = *this * other;
    swap(tmp);
    return *this;
}

affine_form affine_form::operator/(affine_form const& other) const
{
    if (this == &other) {
        return affine_form(context(), 1.0);
    }
    return *this * other.inv();
}

affine_form& affine_form::operator/=(affine_form const& other)
{
    auto tmp = *this / other;
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
        alpha = r > limits::minrad ? (fb - fa) / (b - a)
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

affine_form affine_form::operator^(double exponent) const
{
    if (exponent == 1.0)
        return *this;

    if (exponent == 0.0)
        return affine_form(context(), 1.0);

    auto alpha = 0.0;
    auto beeta = 0.0;
    auto gamma = 0.0;

    auto fMin = std::pow(min(), exponent);
    auto fMax = std::pow(max(), exponent);

    bool exponent_in_01 = affine_interval(0, 1).contains(exponent);

    if (context().approximation_mode() == approximation_mode::CHEBYSHEV) {
        beeta = (fMax - fMin) / (2 * radius());
        auto x2 = std::pow(beeta / exponent, 1 / (exponent - 1));
        alpha = 0.5 * (-beeta * (min() + x2) + fMin + std::pow(x2, exponent));
        gamma = 0.5 * (beeta * (x2 - min()) + fMin - std::pow(x2, exponent));
        if (exponent_in_01)
            gamma = -gamma;
    } else if (context().approximation_mode() == approximation_mode::MINRANGE) {
        beeta = exponent * std::pow(exponent_in_01 ? max() : min(), exponent - 1);
        alpha = 0.5 * (-beeta * 2 * center() + fMin + fMax);
        gamma = 0.5 * (-beeta * 2 * radius() - fMin + fMax);
    }

    auto idx = indices_;
    auto dev = deviations_;
    view::as_array(dev) *= beeta;

    if (gamma != 0) {
        idx.push_back(context().increment_last());
        dev.push_back(gamma);
    }

    return affine_form(context(), beeta * center() + alpha, std::move(dev), std::move(idx), idx.size());
}

affine_form affine_form::operator^(affine_form const& other) const
{
    // precondition
    if (min() < 0) {
        throw std::invalid_argument("affine_form::operator^: exponentiation of negative numbers is only possible for integer exponents.");
    }

    if (length() == 0 && other.length() == 0)
        return affine_form(context(), std::pow(center(), other.center()));

    if (length() == 0)
        return center() ^ other;

    if (other.length() == 0)
        return *this ^ other.center();

    // evaluate the edge of the polygon defined by the base and the exponent
    std::vector<size_t> idx;
    idx.reserve(length() + other.length());

    // merge indices and construct deviations vector
    std::set_union(indices_.begin(), indices_.end(),
        other.indices_.begin(), other.indices_.end(),
        std::back_inserter(idx));

    std::vector<double> dev_b(idx.size()); // base deviations
    std::vector<double> dev_e(idx.size()); // exponent deviations
    std::vector<double> eps(idx.size());

    size_t i = 0, j = 0;
    for (size_t k = 0; k < idx.size(); ++k) {
        double v = 0;
        if (i < length() && indices_[i] == idx[k]) {
            v = dev_b[k] = deviations_[i++];
        }
        if (j < other.length() && other.indices_[j] == idx[k]) {
            v = dev_e[k] = other.deviations_[j++];
        }
        eps[k] = (0 < v) - (v < 0); // this computes the signum, which I think was the intention
    }

    // Find minimum and maximum distance between taylor series and exponentiation
    auto x1 = (view::as_array(eps) * view::as_array(dev_b)).sum() + center();
    auto y1 = (view::as_array(eps) * view::as_array(dev_e)).sum() + other.center();

    // two-dimensional Taylor series at the center value
    auto fc = std::pow(center(), other.center());
    auto fx = other.center() * std::pow(center(), other.center() - 1);
    auto fy = fc * std::log(center());

    size_t last_eps = 0;

    auto dmin = std::numeric_limits<double>::max();
    auto dmax = std::numeric_limits<double>::min();

    // why does the outer loop iterate up to 2 * idx.size() ?
    for (i = 0; i < 2 * idx.size(); ++i) {
        // find the outmost segment
        auto phi0 = std::atan2(other.center() - y1, center() - x1);
        auto phi_max = 0;

        for (j = 0; j < idx.size(); ++j) {
            auto v1 = eps[j] * dev_b[j];
            auto v2 = eps[j] * dev_e[j];

            auto phi = std::atan2(-2 * v2, -2 * v1) - phi0;

            if (phi_max < phi && affine_interval(0, M_PI).contains(phi)) {
                phi_max = phi;
                last_eps = j;
            }
        }

        auto v1 = eps[last_eps] * dev_b[last_eps];
        auto v2 = eps[last_eps] * dev_e[last_eps];

        auto x2 = x1 - v1, x3 = x2 - v1;
        auto y2 = y1 - v2, y3 = y2 - v2;

        auto d1 = fc + fx * (x1 - center()) + fy * (y1 - other.center()) - std::pow(x1, y1);
        auto d2 = fc + fx * (x2 - center()) + fy * (y2 - other.center()) - std::pow(x2, y2);
        auto d3 = fc + fx * (x3 - center()) + fy * (y3 - other.center()) - std::pow(x3, y3);

        auto a = std::min(d1, d3);
        auto b = std::max(d1, d3);

        dmin = std::min(dmin, a);
        dmax = std::max(dmax, b);

        if (!affine_interval(a, b).contains(d2)) {
            if (dev_b[last_eps] == 0) {
                auto x1log = std::log(x1);
                if (fy / x1log > 0) {
                    auto dyc = std::log(fy / x1log) / x1log;
                    d2 = fc + fx * (x1 - center()) + fy * (dyc - other.center()) - std::pow(x1, dyc);
                }
            } else if (dev_e[last_eps] == 0) {
                if (fx / y1 > 0) {
                    auto dyc = std::pow(fx / y1, 1 / (y1 - 1));
                    d2 = fc + fx * (dyc - center()) + fy * (y1 - other.center()) - std::pow(dyc, y1);
                }
            } else {
                // information may have been lost if d2 < dmin or d2 > dmax
            }

            if (d3 < d2) {
                dmax = std::max(dmax, d2);
            }
            if (d1 > d2) {
                dmin = std::min(dmin, d2);
            }
        }

        x1 = x3;
        y1 = y3;

        eps[last_eps] = -eps[last_eps];
    }

    auto alpha = (dmax + dmin) / 2;
    auto gamma = (dmax - dmin) / 2;

    std::cout << "dmin: " << dmin << ", dmax: " << dmax
              << ", alpha: " << alpha << ", gamma: " << gamma << "\n";

    std::vector<double> dev(idx.size());
    view::as_array(dev) = fx * view::as_array(dev_b) + fy * view::as_array(dev_e);

    idx.push_back(context().increment_last());
    dev.push_back(gamma);

    return affine_form(context(), std::pow(center(), other.center()) + alpha, dev, idx, idx.size());
}
} // namespace
