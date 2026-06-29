#ifndef PAPPUS_AF_HPP
#define PAPPUS_AF_HPP

#include <algorithm>
#include <cmath>
#include <concepts>
#include <functional>
#include <numeric>
#include <optional>
#include <stdexcept>

#include <gch/small_vector.hpp>

#include "context.hpp"
#include "interval/interval.hpp"

#ifndef PROTECTED_DIVISION
#define PROTECTED_DIVISION 1
#endif

namespace pappus {

template<std::floating_point T = double>
struct limits {
    static constexpr T eps    = std::numeric_limits<T>::epsilon();
    static constexpr T minrad = T(1e-10);
};

template<std::floating_point T = double>
class affine_form {
public:
    struct term {
        size_t index;
        T      value;
    };
    using terms_t = gch::small_vector<term, 16>;

    explicit affine_form(affine_context const& ctx, interval<T> const& iv)
        : context_(ctx)
        , center_(iv.mid())
        , radius_(0)
    {
        if (iv.radius() != T(0))
            terms_.push_back({context().increment_last(), iv.radius()});
        update_radius();
    }

    explicit affine_form(affine_context const& ctx, T v)
        : context_(ctx)
        , center_(v)
        , radius_(0)
    {
    }

    affine_form(const affine_form& other) = default;
    affine_form(affine_form&& other)      = default;

    void swap(affine_form& other)
    {
        std::swap(center_, other.center_);
        std::swap(radius_, other.radius_);
        terms_.swap(other.terms_);
    }

    T operator[](std::size_t i) const
    {
        return i == 0 ? center_ : terms_[i - 1].value;
    }

    affine_context& context() { return context_; }
    affine_context const& context() const { return context_; }

    T center() const { return center_; }
    T radius() const { return radius_; }
    T max()    const { return center() + radius(); }
    T min()    const { return center() - radius(); }

    T abs_max() const
    {
        return std::fmax(std::fabs(min()), std::fabs(max()));
    }

    T abs_min() const
    {
        auto lo = min(), hi = max();
        if (std::signbit(lo) == std::signbit(hi))
            return std::fmin(std::fabs(lo), std::fabs(hi));
        return T(0);
    }

    size_t size()       const { return terms_.size(); }
    size_t length()     const { return terms_.size(); }
    size_t last_index() const
    {
        if (terms_.empty())
            throw std::logic_error("affine_form::last_index: constant form has no noise terms");
        return terms_.back().index;
    }

    bool operator<(const affine_form& other)  const { return max() < other.max(); }
    bool operator<=(const affine_form& other) const { return max() <= other.max(); }
    bool operator>(const affine_form& other)  const { return !(*this <= other); }
    bool operator>=(const affine_form& other) const { return !(*this < other); }

    bool operator==(const affine_form& other) const
    {
        if (terms_.size() != other.terms_.size()) return false;

        auto not_equal = [](T x, T y) -> bool {
            auto a = std::fabs(x), b = std::fabs(y), c = std::fabs(x - y);
            if (!(a < T(1) && b < T(1))) c /= (a + b);
            return c > limits<T>::eps;
        };

        if (not_equal(center_, other.center_)) return false;
        for (size_t i = 0; i < terms_.size(); ++i) {
            if (terms_[i].index != other.terms_[i].index)           return false;
            if (not_equal(terms_[i].value, other.terms_[i].value))  return false;
        }
        return true;
    }

    interval<T> to_interval() const
    {
        return interval<T>(min(), max());
    }

    affine_form& operator=(T v)
    {
        center_ = v;
        radius_ = 0;
        terms_.clear();
        return *this;
    }

    affine_form& operator=(affine_form other)
    {
        swap(other);
        return *this;
    }

    // ---------------------------------------------------------------------------
    // linear operations (affine-scalar)
    // ---------------------------------------------------------------------------

    affine_form operator+(T v) const
    {
        affine_form f(*this);
        f += v;
        return f;
    }

    affine_form operator-(T v) const
    {
        affine_form f(*this);
        f -= v;
        return f;
    }

    affine_form operator*(T v) const
    {
        affine_form f(*this);
        f *= v;
        return f;
    }

    affine_form operator/(T v) const
    {
        affine_form f(*this);
        f /= v;
        return f;
    }

    affine_form& operator+=(T v)
    {
        center_ += v;
        return *this;
    }

    affine_form& operator-=(T v)
    {
        center_ -= v;
        return *this;
    }

    affine_form& operator*=(T v)
    {
        for (auto& t : terms_) t.value *= v;
        center_ *= v;
        radius_ *= std::fabs(v);
        return *this;
    }

    affine_form& operator/=(T v)
    {
        return operator*=(T(1) / v);
    }

    affine_form operator-() const
    {
        affine_form tmp(*this);
        tmp *= T(-1);
        return tmp;
    }

    // ---------------------------------------------------------------------------
    // linear operations (affine-affine)
    // ---------------------------------------------------------------------------

    affine_form operator+(affine_form const& other) const
    {
        ensure_same_context(other);
        auto terms = merge_terms(terms_, other.terms_,
            [](term const& t){ return t.value; },
            [](term const& t){ return t.value; },
            [](term const& l, term const& r){ return l.value + r.value; });
        return affine_form(context(), center_ + other.center_, std::move(terms));
    }

    affine_form& operator+=(affine_form const& other)
    {
        auto tmp = *this + other;
        swap(tmp);
        return *this;
    }

    affine_form operator-(affine_form const& other) const
    {
        ensure_same_context(other);
        auto terms = merge_terms(terms_, other.terms_,
            [](term const& t){ return t.value; },
            [](term const& t){ return -t.value; },
            [](term const& l, term const& r){ return l.value - r.value; });
        return affine_form(context(), center_ - other.center_, std::move(terms));
    }

    affine_form& operator-=(affine_form const& other)
    {
        auto tmp = *this - other;
        swap(tmp);
        return *this;
    }

    affine_form operator*(affine_form const& other) const
    {
        ensure_same_context(other);
        if (terms_.empty() && other.terms_.empty())
            return affine_form(context(), center_ * other.center_);
        if (terms_.empty()) return other * center_;
        if (other.terms_.empty()) return *this * other.center_;

        auto c1 = center_, c2 = other.center_;
        T common_center = T(0), common_deviation = T(0);

        auto terms = merge_terms(terms_, other.terms_,
            [c2](term const& t){ return c2 * t.value; },
            [c1](term const& t){ return c1 * t.value; },
            [&](term const& l, term const& r) {
                common_center    += l.value * r.value;
                common_deviation += std::fabs(l.value * r.value);
                return c1 * r.value + c2 * l.value;
            });

        common_center    /= T(2);
        common_deviation /= T(2);
        auto delta = radius_ * other.radius_;

        if (context().approximation_mode() == approximation_mode::SECANT) {
            delta -= common_deviation;
            auto r = std::transform_reduce(terms.begin(), terms.end(), T(0),
                std::plus<T>{}, [](term const& t){ return std::fabs(t.value); });
            T fac = r < limits<T>::eps ? T(1) : T(1) + delta / r;
            for (auto& t : terms) t.value *= fac;
            terms.push_back({context().increment_last(), T(0)});
        } else {
            terms.push_back({context().increment_last(), delta - common_deviation});
        }
        return affine_form(context(), c1 * c2 + common_center, std::move(terms));
    }

    affine_form& operator*=(affine_form const& other)
    {
        auto tmp = *this * other;
        swap(tmp);
        return *this;
    }

    affine_form operator/(affine_form const& other) const
    {
        if (this == &other) return affine_form(context(), T(1));
        return *this * other.inv();
    }

    affine_form& operator/=(affine_form const& other)
    {
        auto tmp = *this / other;
        swap(tmp);
        return *this;
    }

    // ---------------------------------------------------------------------------
    // non-linear operations
    // ---------------------------------------------------------------------------

    affine_form inv() const
    {
        // Note: unlike interval::inv() which returns [-inf, +inf] for an interval
        // strictly containing zero (and a half-infinite interval for one that
        // touches zero), affine_form::inv() throws whenever [min, max] contains
        // zero. Affine forms cannot represent unbounded values, and any
        // non-trivial interval that contains zero yields an unbounded inverse.
        if (terms_.empty()) {
            if (center_ == T(0))
                throw std::invalid_argument("affine_form::inv: zero is not invertible");
            return affine_form(context(), T(1) / center_);
        }

        auto c = center(), r = radius();
        auto a = c - r, b = c + r;
        if (a <= T(0) && T(0) <= b)
            throw std::invalid_argument("affine_form::inv: interval containing zero is not invertible");
        auto fa = T(1) / a, fb = T(1) / b;

        T alpha = 0, delta = 0, dzeta = 0;

        switch (context().approximation_mode()) {
        case approximation_mode::CHEBYSHEV: {
            auto u = std::sqrt(a * b);
            alpha = -fa * fb;
            if (a > T(0)) {
                delta = +T(0.5) * (fa + fb - T(2) / u);
                dzeta =  fa + fb - delta;
            } else {
                delta = -T(0.5) * (fa + fb + T(2) / u);
                dzeta =  fa + fb + delta;
            }
            break;
        }
        case approximation_mode::MINRANGE: {
            T ya, yb;
            if (a > T(0)) {
                alpha = -fb / b;
                ya = fa - alpha * a;
                yb = T(2) * fb;
            } else {
                alpha = -fa / a;
                ya = T(2) * fa;
                yb = fb - alpha * b;
            }
            delta = T(0.5) * (ya - yb);
            dzeta = T(0.5) * (ya + yb);
            break;
        }
        case approximation_mode::SECANT: {
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : -fa * fb;
            dzeta = fa - alpha * a;
            delta = T(0);
            break;
        }
        }

        auto result_terms = terms_;
        for (auto& t : result_terms) t.value *= alpha;
        result_terms.push_back({context().increment_last(), delta});
        return affine_form(context(), alpha * c + dzeta, std::move(result_terms));
    }

    affine_form pow(int exponent) const
    {
        if (terms_.empty())
            return affine_form(context(), std::pow(center(), exponent));
        if (exponent == 0) return affine_form(context(), T(1));
        if (exponent == 1) return *this;
        if (exponent == -1) return this->inv();

        auto c = center(), r = radius();
        auto a = c - r, b = c + r;
        auto fa = std::pow(a, exponent), fb = std::pow(b, exponent);

        T alpha = T(0), delta, dzeta;

        switch (context().approximation_mode()) {
        case approximation_mode::CHEBYSHEV: {
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(exponent) * fa / a;
            auto e  = std::fabs(alpha / exponent);
            auto x1 = -std::pow(e, T(1) / T(exponent - 1));
            auto x2 = -x1;
            auto y1 = x1 < a ? fa - alpha * a : std::pow(x1, exponent) - alpha * x1;
            auto y2 = x2 > b ? fb - alpha * b : std::pow(x2, exponent) - alpha * x2;
            delta = (y1 - y2) * T(0.5);
            dzeta = (y1 + y2) * T(0.5);
            break;
        }
        case approximation_mode::MINRANGE: {
            if (a * b < T(0)) {
                alpha = T(0);
                if (exponent % 2 == 0) {
                    delta = T(0.5) * std::max(fa, fb);
                    dzeta = delta;
                } else {
                    delta = T(0.5) * (fb - fa);
                    dzeta = T(0.5) * (fb + fa);
                }
            } else {
                alpha = std::signbit(a) == std::signbit(exponent)
                            ? T(exponent) * fa / a
                            : T(exponent) * fb / b;
                auto ya = fa - alpha * a, yb = fb - alpha * b;
                delta = T(0.5) * std::fabs(ya - yb);
                dzeta = T(0.5) * (fa + fb);
            }
            break;
        }
        case approximation_mode::SECANT: {
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(exponent) * fa / a;
            delta = T(0);
            dzeta = fa - alpha * a;
            break;
        }
        }

        auto result_terms = terms_;
        for (auto& t : result_terms) t.value *= alpha;
        result_terms.push_back({context().increment_last(), delta});
        return affine_form(context(), alpha * c + dzeta, std::move(result_terms));
    }

    affine_form pow(T exponent) const
    {
        if (terms_.empty()) {
            if (center_ < T(0) && !(std::isfinite(exponent) && std::trunc(exponent) == exponent))
                throw std::invalid_argument("affine_form::pow: fractional exponent requires nonnegative base");
            if (center_ == T(0) && exponent < T(0))
                throw std::invalid_argument("affine_form::pow: negative exponent requires strictly positive base");
            return affine_form(context(), std::pow(center_, exponent));
        }
        if (exponent == T(1)) return *this;
        if (exponent == T(0)) return affine_form(context(), T(1));

        auto is_integer_exponent = std::isfinite(exponent) && std::trunc(exponent) == exponent;
        if (min() < T(0))
            if (!is_integer_exponent)
                throw std::invalid_argument("affine_form::pow: fractional exponent requires nonnegative base");
        if (exponent < T(0) && min() <= T(0))
            throw std::invalid_argument("affine_form::pow: negative exponent requires strictly positive base");

        T alpha = T(0), beta = T(0), gamma = T(0);
        auto fMin = std::pow(min(), exponent);
        auto fMax = std::pow(max(), exponent);
        bool exp_in_01 = interval<T>(T(0), T(1)).contains(exponent);

        switch (context().approximation_mode()) {
        case approximation_mode::CHEBYSHEV: {
            beta  = (fMax - fMin) / (T(2) * radius());
            auto x2 = std::pow(beta / exponent, T(1) / (exponent - T(1)));
            alpha = T(0.5) * (-beta * (min() + x2) + fMin + std::pow(x2, exponent));
            gamma = T(0.5) * (-beta * (min() - x2) + fMin - std::pow(x2, exponent));
            if (exp_in_01) gamma = -gamma;
            break;
        }
        case approximation_mode::MINRANGE: {
            beta  = exponent * std::pow(exp_in_01 ? max() : min(), exponent - T(1));
            alpha = T(0.5) * (-beta * T(2) * center() + fMin + fMax);
            gamma = T(0.5) * (-beta * T(2) * radius() - fMin + fMax);
            break;
        }
        case approximation_mode::SECANT: {
            // Linear interpolation through (min, fMin) and (max, fMax); no noise term.
            beta  = radius() > limits<T>::minrad ? (fMax - fMin) / (T(2) * radius()) : exponent * std::pow(min(), exponent - T(1));
            alpha = T(0.5) * (-beta * T(2) * center() + fMin + fMax);
            gamma = T(0);
            break;
        }
        }

        auto result_terms = terms_;
        for (auto& t : result_terms) t.value *= beta;
        if (gamma != T(0))
            result_terms.push_back({context().increment_last(), gamma});
        return affine_form(context(), beta * center() + alpha, std::move(result_terms));
    }

    affine_form pow(affine_form const& other) const
    {
        ensure_same_context(other);
        if (min() < T(0))
            throw std::invalid_argument("affine_form::pow: exponentiation of negative base requires integer exponent");

        if (terms_.empty() && other.terms_.empty())
            return affine_form(context(), std::pow(center(), other.center()));
        if (terms_.empty()) return pow(center(), other);
        if (other.terms_.empty()) return this->pow(other.center());

        // Build sorted union of indices
        gch::small_vector<size_t, 16> idx;
        idx.reserve(terms_.size() + other.terms_.size());
        {
            auto li = terms_.begin(), ri = other.terms_.begin();
            while (li != terms_.end() && ri != other.terms_.end()) {
                if (li->index < ri->index)      { idx.push_back(li->index); ++li; }
                else if (ri->index < li->index) { idx.push_back(ri->index); ++ri; }
                else                            { idx.push_back(li->index); ++li; ++ri; }
            }
            for (; li != terms_.end();       ++li) idx.push_back(li->index);
            for (; ri != other.terms_.end(); ++ri) idx.push_back(ri->index);
        }

        gch::small_vector<T, 16> dev_b(idx.size(), T(0)), dev_e(idx.size(), T(0)), eps_v(idx.size(), T(0));
        {
            auto li = terms_.begin(), ri = other.terms_.begin();
            for (size_t k = 0; k < idx.size(); ++k) {
                T v = T(0);
                if (li != terms_.end()       && li->index == idx[k]) { v = dev_b[k] = (li++)->value; }
                if (ri != other.terms_.end() && ri->index == idx[k]) { v = dev_e[k] = (ri++)->value; }
                eps_v[k] = T((v < T(0)) - (T(0) < v));
            }
        }

        auto x1 = std::inner_product(eps_v.begin(), eps_v.end(), dev_b.begin(), T(0)) + center();
        auto y1 = std::inner_product(eps_v.begin(), eps_v.end(), dev_e.begin(), T(0)) + other.center();

        auto fc = std::pow(center(), other.center());
        auto fx = other.center() * std::pow(center(), other.center() - T(1));
        auto fy = fc * std::log(center());

        size_t last_eps = 0;
        auto dmin = std::numeric_limits<T>::max();
        auto dmax = std::numeric_limits<T>::lowest();

        for (size_t i = 0; i < 2 * idx.size(); ++i) {
            auto phi0   = std::atan2(other.center() - y1, center() - x1);
            auto phi_max = T(0);

            for (size_t j = 0; j < idx.size(); ++j) {
                auto phi = std::atan2(dev_e[j] == T(0) ? T(0) : -T(2) * eps_v[j] * dev_e[j],
                                      dev_b[j] == T(0) ? T(0) : -T(2) * eps_v[j] * dev_b[j]);
                phi -= phi0;
                if (phi_max < phi && interval<T>(T(0), fp::pi_v<T>).contains(phi)) {
                    phi_max  = phi;
                    last_eps = j;
                }
            }

            auto v1 = eps_v[last_eps] * dev_b[last_eps];
            auto v2 = eps_v[last_eps] * dev_e[last_eps];
            auto x2 = x1 - v1, x3 = x2 - v1;
            auto y2 = y1 - v2, y3 = y2 - v2;

            auto d1 = fc + fx * (x1 - center()) + fy * (y1 - other.center()) - std::pow(x1, y1);
            auto d2 = fc + fx * (x2 - center()) + fy * (y2 - other.center()) - std::pow(x2, y2);
            auto d3 = fc + fx * (x3 - center()) + fy * (y3 - other.center()) - std::pow(x3, y3);

            auto aa = std::min(d1, d3), bb = std::max(d1, d3);
            dmin = std::min(dmin, aa);
            dmax = std::max(dmax, bb);

            if (!interval<T>(aa, bb).contains(d2)) {
                if (dev_b[last_eps] == T(0)) {
                    auto x1log = std::log(x1);
                    if (fy / x1log > T(0)) {
                        auto dyc = std::log(fy / x1log) / x1log;
                        d2 = fc + fx * (x1 - center()) + fy * (dyc - other.center()) - std::pow(x1, dyc);
                    }
                } else if (dev_e[last_eps] == T(0)) {
                    if (fx / y1 > T(0)) {
                        auto dyc = std::pow(fx / y1, T(1) / (y1 - T(1)));
                        d2 = fc + fx * (dyc - center()) + fy * (y1 - other.center()) - std::pow(dyc, y1);
                    }
                }
                dmin = std::min(dmin, std::min(aa, d2));
                dmax = std::max(dmax, std::max(bb, d2));
            }

            x1 = x3; y1 = y3;
            eps_v[last_eps] = -eps_v[last_eps];
        }

        auto alpha = (dmax + dmin) / T(2);
        auto gamma = (dmax - dmin) / T(2);

        terms_t result_terms;
        result_terms.reserve(idx.size() + 1);
        for (size_t k = 0; k < idx.size(); ++k)
            result_terms.push_back({idx[k], fx * dev_b[k] + fy * dev_e[k]});
        result_terms.push_back({context().increment_last(), gamma});

        return affine_form(context(), std::pow(center(), other.center()) + alpha,
                           std::move(result_terms));
    }

    static affine_form pow(T base, affine_form const& exponent)
    {
        if (exponent.terms_.empty())
            return affine_form(exponent.context(), std::pow(base, exponent.center()));
        if (base == T(1)) {
            auto result_terms = exponent.terms_;
            result_terms.push_back({exponent.context().increment_last(), T(0)});
            return affine_form(exponent.context(), T(1), std::move(result_terms));
        }
        if (base == T(0))
            throw std::invalid_argument("affine_form::pow: base cannot be zero");

        T alpha = T(0), beta = T(0), gamma = T(0);
        auto fMin = std::pow(base, exponent.min());
        auto fMax = std::pow(base, exponent.max());

        if (exponent.context().approximation_mode() == approximation_mode::MINRANGE) {
            beta  = fMin * std::log(base);
            alpha = -beta * T(2) * exponent.center() + fMin + fMax;
            gamma = -beta * T(2) * exponent.radius() - fMin + fMax;
        } else { // CHEBYSHEV
            auto b = std::log(base);
            beta  = (fMax - fMin) / (exponent.max() - exponent.min());
            auto x2 = std::log(beta / b) / b;
            alpha = T(0.5) * (-beta * (exponent.min() + x2) + fMin + std::pow(base, x2));
            gamma = T(0.5) * (-beta * (exponent.min() - x2) + fMin - std::pow(base, x2));
        }

        auto result_terms = exponent.terms_;
        for (auto& t : result_terms) t.value *= beta;
        result_terms.push_back({exponent.context().increment_last(), gamma});
        return affine_form(exponent.context(), beta * exponent.center() + alpha,
                           std::move(result_terms));
    }

    static affine_form pow(affine_form const& af, T v) { return af.pow(v); }

    // abs(x): V-shaped, non-differentiable at 0.
    // - domain entirely non-negative: abs(x) = x (identity)
    // - domain entirely non-positive: abs(x) = -x (negation)
    // - domain crosses 0: secant from (a, |a|) to (b, |b|) with alpha = (|b|-|a|)/(b-a),
    //   delta covers the V's deviation from the secant. For a symmetric domain [-r, r],
    //   this degenerates to alpha=0, dzeta=r, delta=0 (a constant r), which is the
    //   tightest possible affine enclosure.
    affine_form abs() const
    {
        if (terms_.empty())
            return affine_form(context(), std::fabs(center()));

        auto a = min(), b = max();
        if (a >= T(0)) return *this;             // entirely non-negative
        if (b <= T(0)) return -*this;            // entirely non-positive

        // domain crosses zero
        auto fa = std::fabs(a), fb = std::fabs(b);
        T alpha, dzeta, delta;
        switch (context().approximation_mode()) {
        case approximation_mode::SECANT:
            alpha = (fb - fa) / (b - a);
            dzeta = fa - alpha * a;
            delta = T(0);
            break;
        default: {
            // Chebyshev / MINRANGE: the V-shape f(x) = |x| on [a, 0] (decreasing to 0)
            // and [0, b] (increasing from 0). The secant connects (a, fa) and (b, fb).
            // The maximum deviation of |x| from the secant occurs at x=0 where
            // f(0)=0 but secant(0) = dzeta. So delta = |dzeta| and we split the
            // difference to get a symmetric enclosure.
            alpha = (fb - fa) / (b - a);
            auto secant_at_zero = fa - alpha * a; // = fb - alpha * b
            delta = std::fabs(secant_at_zero);
            dzeta = secant_at_zero - delta;       // shift so the enclosure is centered
            break;
        }
        }
        return apply_unary(center(), alpha, dzeta, delta);
    }

    affine_form cbrt() const
    {
        if (terms_.empty())
            return affine_form(context(), std::cbrt(center()));

        auto a = min(), b = max();
        auto c = center(), r = radius();
        auto fa = std::cbrt(a), fb = std::cbrt(b);
        T alpha, dzeta, delta;
        switch (context().approximation_mode()) {
        case approximation_mode::SECANT:
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) / (T(3) * std::cbrt(c * c));
            dzeta = fa - alpha * a;
            delta = T(0);
            break;
        default: {
            // cbrt is concave on (0, inf) and convex on (-inf, 0).
            // For domains not crossing 0, use the Chebyshev approach:
            // tangent point where f'(x*) = alpha, i.e. x* = (1/(3*alpha))^(3/2)
            // For simplicity and soundness, use the secant + max deviation.
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) / (T(3) * std::cbrt(c * c));
            auto g = [&](T x) { return std::cbrt(x) - alpha * x; };
            auto ga = g(a), gb = g(b);
            auto g_max = std::max(ga, gb), g_min = std::min(ga, gb);
            // Check the critical point where f'(x) = alpha: 1/(3*cbrt(x^2)) = alpha
            // => cbrt(x^2) = 1/(3*alpha) => x^2 = 1/(3*alpha)^3 => x = ±1/(3*alpha)^(3/2)
            if (alpha != T(0)) {
                auto x_crit = T(1) / std::pow(T(3) * alpha, T(1.5));
                if (x_crit > a && x_crit < b) {
                    auto gx = g(x_crit);
                    g_max = std::max(g_max, gx);
                    g_min = std::min(g_min, gx);
                }
                x_crit = -x_crit;
                if (x_crit > a && x_crit < b) {
                    auto gx = g(x_crit);
                    g_max = std::max(g_max, gx);
                    g_min = std::min(g_min, gx);
                }
            }
            delta = T(0.5) * (g_max - g_min);
            dzeta = T(0.5) * (g_max + g_min);
            break;
        }
        }
        return apply_unary(c, alpha, dzeta, delta);
    }

    affine_form log1p() const
    {
        if (terms_.empty()) {
            if (center_ <= T(-1))
                throw std::invalid_argument("affine_form::log1p: argument <= -1");
            return affine_form(context(), std::log1p(center_));
        }
        auto a = min(), b = max();
        if (a <= T(-1))
            throw std::invalid_argument("affine_form::log1p: interval contains values <= -1");
        auto c = center(), r = radius();
        auto fa = std::log1p(a), fb = std::log1p(b);
        T alpha, dzeta, delta;
        switch (context().approximation_mode()) {
        case approximation_mode::CHEBYSHEV: {
            // log1p is concave, same shape as log shifted by 1.
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) / (T(1) + c);
            auto ga = fa - alpha * a;
            // tangent point: f'(x*) = alpha => 1/(1+x*) = alpha => x* = 1/alpha - 1
            auto gx = -std::log(alpha) - T(1);
            delta = T(0.5) * (gx - ga);    // concave: f above secant, gx > ga
            dzeta = T(0.5) * (gx + ga);
            break;
        }
        case approximation_mode::MINRANGE: {
            alpha = T(1) / (T(1) + b);
            auto ya = fa - alpha * a, yb = fb - alpha * b;
            delta = T(0.5) * std::fabs(yb - ya);
            dzeta = T(0.5) * (ya + yb);
            break;
        }
        default:
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) / (T(1) + c);
            dzeta = fa - alpha * a;
            delta = T(0);
        }
        return apply_unary(c, alpha, dzeta, delta);
    }

    // floor/ceil are step functions — not differentiable. The best affine
    // enclosure uses the secant (linear interpolation between endpoints)
    // with a delta covering the maximum step deviation.
    affine_form floor() const
    {
        if (terms_.empty())
            return affine_form(context(), std::floor(center()));

        auto a = min(), b = max();
        auto fa = std::floor(a), fb = std::floor(b);
        // Secant slope; for a domain within a single floor interval, slope = 0.
        T alpha = (b > a) ? (fb - fa) / (b - a) : T(0);
        // Maximum deviation: the floor function differs from the secant by at
        // most 1 unit. Use the secant with delta = max deviation.
        auto g = [&](T x) { return std::floor(x) - alpha * x; };
        // The extrema of g occur at the endpoints and at integer boundaries.
        // For simplicity, use endpoint deviation + 1 (worst-case step).
        auto ga = g(a), gb = g(b);
        auto g_max = std::max(ga, gb), g_min = std::min(ga, gb);
        // floor can deviate by up to 1 from any linear approximation within [a,b]
        g_max += T(1);
        g_min -= T(0);
        T delta = T(0.5) * (g_max - g_min);
        T dzeta = T(0.5) * (g_max + g_min);
        return apply_unary(center(), alpha, dzeta, delta);
    }

    affine_form ceil() const
    {
        if (terms_.empty())
            return affine_form(context(), std::ceil(center()));

        auto a = min(), b = max();
        auto fa = std::ceil(a), fb = std::ceil(b);
        T alpha = (b > a) ? (fb - fa) / (b - a) : T(0);
        auto g = [&](T x) { return std::ceil(x) - alpha * x; };
        auto ga = g(a), gb = g(b);
        auto g_max = std::max(ga, gb), g_min = std::min(ga, gb);
        g_max += T(0);
        g_min -= T(1);
        T delta = T(0.5) * (g_max - g_min);
        T dzeta = T(0.5) * (g_max + g_min);
        return apply_unary(center(), alpha, dzeta, delta);
    }

    affine_form sqrt() const
    {
        if (terms_.empty())
            return affine_form(context(), std::sqrt(center()));

        auto a = min(), b = max();
        if (a < T(0))
            throw std::runtime_error("affine_form::sqrt: negative argument");

        auto fa = std::sqrt(a), fb = std::sqrt(b);
        T alpha, dzeta, delta;

        switch (context().approximation_mode()) {
        case approximation_mode::CHEBYSHEV: {
            alpha = T(1) / (fa + fb);
            auto t = T(0.25) * fa * fb * alpha;
            auto u = T(0.125) * (a + b) * alpha;
            dzeta = u + T(3) * t;
            delta = u - t;
            break;
        }
        case approximation_mode::MINRANGE: {
            alpha = T(1) / (T(2) * fb);
            delta = T(0.5) * alpha * (a - T(2) * fa * fb + b);
            dzeta = T(0.5) * fb - delta;
            break;
        }
        case approximation_mode::SECANT: {
            alpha = radius() > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) / (T(2) * fb);
            dzeta = fa - alpha * a;
            delta = T(0);
            break;
        }
        }

        auto result_terms = terms_;
        for (auto& t : result_terms) t.value *= alpha;
        result_terms.push_back({context().increment_last(), delta});
        return affine_form(context(), center() * alpha + dzeta, std::move(result_terms));
    }

    affine_form isqrt() const
    {
        if (terms_.empty())
            return affine_form(context(), T(1) / std::sqrt(center()));

        auto a = min(), b = max();
        if (a < T(0) || b < T(0))
            throw std::runtime_error("affine_form::isqrt: negative argument");

        auto fa = T(1) / std::sqrt(a), fb = T(1) / std::sqrt(b);
        T alpha, dzeta, delta;

        switch (context().approximation_mode()) {
        case approximation_mode::CHEBYSHEV: {
            alpha = radius() > limits<T>::minrad
                        ? (fb - fa) / (b - a)
                        : T(-0.5) * fb * fb * fb;
            auto x = std::pow(T(0.5) * (a / fb + b / fa), T(2) / T(3));
            delta = T(0.5) * (fa + alpha * (x - a) - T(1) / std::sqrt(x));
            dzeta = fa - alpha * a - delta;
            break;
        }
        case approximation_mode::MINRANGE: {
            alpha = T(-0.5) * fb * fb * fb;
            delta = T(0.5) * (fa - fb + alpha * (b - a));
            dzeta = fa - alpha * a - delta;
            break;
        }
        case approximation_mode::SECANT: {
            alpha = radius() > limits<T>::minrad
                        ? (fb - fa) / (b - a)
                        : T(-0.5) * fb * fb * fb;
            dzeta = fa - alpha * a;
            delta = T(0);
            break;
        }
        }

        auto result_terms = terms_;
        for (auto& t : result_terms) t.value *= alpha;
        result_terms.push_back({context().increment_last(), delta});
        return affine_form(context(), center() * alpha + dzeta, std::move(result_terms));
    }

    // ---------------------------------------------------------------------------
    // transcendental functions
    // ---------------------------------------------------------------------------

    affine_form exp() const
    {
        if (terms_.empty())
            return affine_form(context(), std::exp(center_));
        auto c = center(), r = radius(), a = c - r, b = c + r;
        auto fa = std::exp(a), fb = std::exp(b);
        T alpha, dzeta, delta;
        switch (context().approximation_mode()) {
        case approximation_mode::CHEBYSHEV: {
            // convex: tangent point x* = log(alpha), g(x*) = alpha*(1-log(alpha))
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : fa;
            auto ga = fa - alpha * a;
            auto gx = alpha * (T(1) - std::log(alpha));
            delta = T(0.5) * (ga - gx);    // convex: secant above f, ga > gx
            dzeta = T(0.5) * (ga + gx);
            break;
        }
        case approximation_mode::MINRANGE: {
            alpha = fa;                     // f'(a) = exp(a), min slope (convex)
            auto ya = fa - alpha * a, yb = fb - alpha * b;
            delta = T(0.5) * (yb - ya);
            dzeta = T(0.5) * (ya + yb);
            break;
        }
        default:
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : fa;
            dzeta = fa - alpha * a; delta = T(0);
        }
        return apply_unary(c, alpha, dzeta, delta);
    }

    affine_form log() const
    {
        if (terms_.empty()) {
            if (center_ <= T(0))
                throw std::invalid_argument("affine_form::log: non-positive argument");
            return affine_form(context(), std::log(center_));
        }
        auto a = min(), b = max();
        if (a <= T(0))
            throw std::invalid_argument("affine_form::log: interval contains non-positive values");
        auto c = center(), r = radius();
        auto fa = std::log(a), fb = std::log(b);
        T alpha, dzeta, delta;
        switch (context().approximation_mode()) {
        case approximation_mode::CHEBYSHEV: {
            // concave: tangent point x* = 1/alpha, g(x*) = -log(alpha) - 1
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) / c;
            auto ga = fa - alpha * a;
            auto gx = -std::log(alpha) - T(1);
            delta = T(0.5) * (gx - ga);    // concave: f above secant, gx > ga
            dzeta = T(0.5) * (gx + ga);
            break;
        }
        case approximation_mode::MINRANGE: {
            alpha = T(1) / b;              // f'(b) = 1/b, min slope (concave, f' decreasing)
            auto ya = fa - alpha * a, yb = fb - alpha * b;
            delta = T(0.5) * std::fabs(yb - ya);
            dzeta = T(0.5) * (ya + yb);
            break;
        }
        default:
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) / c;
            dzeta = fa - alpha * a; delta = T(0);
        }
        return apply_unary(c, alpha, dzeta, delta);
    }

    affine_form sin() const
    {
        if (terms_.empty())
            return affine_form(context(), std::sin(center_));
        auto c = center(), r = radius(), a = c - r, b = c + r;
        if (b - a >= fp::two_pi_v<T>)
            return affine_form(context(), interval<T>(T(-1), T(1)));
        auto fa = std::sin(a), fb = std::sin(b);
        T alpha, dzeta, delta;
        if (context().approximation_mode() == approximation_mode::SECANT) {
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : std::cos(c);
            dzeta = fa - alpha * a; delta = T(0);
        } else {
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : std::cos(c);
            auto g   = [&](T x) { return std::sin(x) - alpha * x; };
            auto ga = g(a), gb = g(b);
            auto g_max = std::max(ga, gb), g_min = std::min(ga, gb);
            if (std::fabs(alpha) <= T(1)) {
                auto phi = std::acos(alpha);           // ∈ [0, π]
                auto two_pi = fp::two_pi_v<T>;
                auto check = [&](T x) {
                    if (x > a && x < b) { auto gx = g(x); g_max = std::max(g_max, gx); g_min = std::min(g_min, gx); }
                };
                // cos(x) = alpha at x = phi + 2kπ and x = -phi + 2kπ
                auto k0 = static_cast<long>(std::ceil((a - phi) / two_pi));
                auto k1 = static_cast<long>(std::floor((b - phi) / two_pi));
                for (auto k = k0; k <= k1; ++k) check(phi + T(k) * two_pi);
                auto k2 = static_cast<long>(std::ceil((a + phi) / two_pi));
                auto k3 = static_cast<long>(std::floor((b + phi) / two_pi));
                for (auto k = k2; k <= k3; ++k) check(-phi + T(k) * two_pi);
            }
            delta = T(0.5) * (g_max - g_min);
            dzeta = T(0.5) * (g_max + g_min);
        }
        return apply_unary(c, alpha, dzeta, delta);
    }

    affine_form cos() const
    {
        if (terms_.empty())
            return affine_form(context(), std::cos(center_));
        auto c = center(), r = radius(), a = c - r, b = c + r;
        if (b - a >= fp::two_pi_v<T>)
            return affine_form(context(), interval<T>(T(-1), T(1)));
        auto fa = std::cos(a), fb = std::cos(b);
        T alpha, dzeta, delta;
        if (context().approximation_mode() == approximation_mode::SECANT) {
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : -std::sin(c);
            dzeta = fa - alpha * a; delta = T(0);
        } else {
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : -std::sin(c);
            auto g   = [&](T x) { return std::cos(x) - alpha * x; };
            auto ga = g(a), gb = g(b);
            auto g_max = std::max(ga, gb), g_min = std::min(ga, gb);
            if (std::fabs(alpha) <= T(1)) {
                // -sin(x) = alpha → sin(x) = -alpha → x = arcsin(-alpha) + 2kπ or (π-arcsin(-alpha)) + 2kπ
                auto phi = std::asin(-alpha);          // ∈ [-π/2, π/2]
                auto pi_mphi = fp::pi_v<T> - phi;
                auto two_pi  = fp::two_pi_v<T>;
                auto check = [&](T x) {
                    if (x > a && x < b) { auto gx = g(x); g_max = std::max(g_max, gx); g_min = std::min(g_min, gx); }
                };
                auto k0 = static_cast<long>(std::ceil((a - phi) / two_pi));
                auto k1 = static_cast<long>(std::floor((b - phi) / two_pi));
                for (auto k = k0; k <= k1; ++k) check(phi + T(k) * two_pi);
                auto k2 = static_cast<long>(std::ceil((a - pi_mphi) / two_pi));
                auto k3 = static_cast<long>(std::floor((b - pi_mphi) / two_pi));
                for (auto k = k2; k <= k3; ++k) check(pi_mphi + T(k) * two_pi);
            }
            delta = T(0.5) * (g_max - g_min);
            dzeta = T(0.5) * (g_max + g_min);
        }
        return apply_unary(c, alpha, dzeta, delta);
    }

    affine_form tan() const
    {
        if (terms_.empty())
            return affine_form(context(), std::tan(center_));
        auto a = min(), b = max();
        // Reject if [a,b] contains an asymptote (k + 1/2)π
        {
            auto half_pi = fp::half_pi_v<T>;
            auto pi      = fp::pi_v<T>;
            if (static_cast<long>(std::floor((a + half_pi) / pi)) !=
                static_cast<long>(std::floor((b + half_pi) / pi)))
                throw std::domain_error("affine_form::tan: interval crosses an asymptote");
        }
        auto c = center(), r = radius();
        auto fa = std::tan(a), fb = std::tan(b);
        T alpha, dzeta, delta;
        if (context().approximation_mode() == approximation_mode::SECANT) {
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) + fa * fa;
            dzeta = fa - alpha * a; delta = T(0);
        } else {
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) + fa * fa;
            auto g   = [&](T x) { auto t = std::tan(x); return t - alpha * x; };
            auto ga = g(a), gb = g(b);
            auto g_max = std::max(ga, gb), g_min = std::min(ga, gb);
            // sec²(x*) = alpha → cos²(x*) = 1/alpha (valid when alpha ≥ 1)
            if (alpha >= T(1)) {
                auto phi = std::acos(T(1) / std::sqrt(alpha));  // ∈ [0, π/2]
                auto two_pi = fp::two_pi_v<T>;
                auto check = [&](T x) {
                    if (x > a && x < b) { auto gx = g(x); g_max = std::max(g_max, gx); g_min = std::min(g_min, gx); }
                };
                for (T off : {phi, -phi}) {
                    auto k0 = static_cast<long>(std::ceil((a - off) / two_pi));
                    auto k1 = static_cast<long>(std::floor((b - off) / two_pi));
                    for (auto k = k0; k <= k1; ++k) check(off + T(k) * two_pi);
                }
            }
            delta = T(0.5) * (g_max - g_min);
            dzeta = T(0.5) * (g_max + g_min);
        }
        return apply_unary(c, alpha, dzeta, delta);
    }

    affine_form sinh() const
    {
        if (terms_.empty())
            return affine_form(context(), std::sinh(center_));
        auto c = center(), r = radius(), a = c - r, b = c + r;
        auto fa = std::sinh(a), fb = std::sinh(b);
        T alpha, dzeta, delta;
        switch (context().approximation_mode()) {
        case approximation_mode::CHEBYSHEV: {
            // f'(x) = cosh(x) ≥ 1; secant slope ≥ 1 by MVT
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : std::cosh(c);
            auto g   = [&](T x) { return std::sinh(x) - alpha * x; };
            auto g_max = std::max(g(a), g(b)), g_min = std::min(g(a), g(b));
            // cosh(x*) = alpha → x* = ±acosh(alpha)
            if (alpha >= T(1)) {
                auto xs = std::acosh(alpha);
                for (T x : {xs, -xs}) {
                    if (x > a && x < b) { auto gx = g(x); g_max = std::max(g_max, gx); g_min = std::min(g_min, gx); }
                }
            }
            delta = T(0.5) * (g_max - g_min);
            dzeta = T(0.5) * (g_max + g_min);
            break;
        }
        case approximation_mode::MINRANGE: {
            // cosh(x) ≥ 1 with minimum at x=0; use slope at point closest to 0
            T x_min = (a >= T(0)) ? a : (b <= T(0)) ? b : T(0);
            alpha = std::cosh(x_min);
            auto ya = fa - alpha * a, yb = fb - alpha * b;
            delta = T(0.5) * std::fabs(yb - ya);
            dzeta = T(0.5) * (ya + yb);
            break;
        }
        default:
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : std::cosh(c);
            dzeta = fa - alpha * a; delta = T(0);
        }
        return apply_unary(c, alpha, dzeta, delta);
    }

    affine_form cosh() const
    {
        if (terms_.empty())
            return affine_form(context(), std::cosh(center_));
        auto c = center(), r = radius(), a = c - r, b = c + r;
        auto fa = std::cosh(a), fb = std::cosh(b);
        T alpha, dzeta, delta;
        switch (context().approximation_mode()) {
        // cosh is convex everywhere (f'' = cosh > 0); tangent: sinh(x*) = alpha → x* = asinh(alpha)
        case approximation_mode::CHEBYSHEV:
        case approximation_mode::MINRANGE: {
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : std::sinh(c);
            auto ga = fa - alpha * a, gb = fb - alpha * b;
            auto xs = std::asinh(alpha);                // unique tangent point
            auto gx = std::sqrt(T(1) + alpha * alpha) - alpha * xs;  // cosh(asinh(alpha)) = sqrt(1+alpha²)
            auto g_max = std::max({ga, gb, gx}), g_min = std::min({ga, gb, gx});
            delta = T(0.5) * (g_max - g_min);
            dzeta = T(0.5) * (g_max + g_min);
            break;
        }
        default:
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : std::sinh(c);
            dzeta = fa - alpha * a; delta = T(0);
        }
        return apply_unary(c, alpha, dzeta, delta);
    }

    affine_form tanh() const
    {
        if (terms_.empty())
            return affine_form(context(), std::tanh(center_));
        auto c = center(), r = radius(), a = c - r, b = c + r;
        auto fa = std::tanh(a), fb = std::tanh(b);
        T alpha, dzeta, delta;
        switch (context().approximation_mode()) {
        case approximation_mode::CHEBYSHEV: {
            // f'(x) = sech²(x) = 1 - tanh²(x); tangent: tanh(x*) = ±sqrt(1-alpha)
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) - fa * fa;
            auto g   = [&](T x) { return std::tanh(x) - alpha * x; };
            auto g_max = std::max(g(a), g(b)), g_min = std::min(g(a), g(b));
            if (alpha >= T(0) && alpha <= T(1)) {
                auto t = std::sqrt(T(1) - alpha);
                for (T tv : {t, -t}) {
                    auto xs = std::atanh(tv);
                    if (xs > a && xs < b) { auto gx = tv - alpha * xs; g_max = std::max(g_max, gx); g_min = std::min(g_min, gx); }
                }
            }
            delta = T(0.5) * (g_max - g_min);
            dzeta = T(0.5) * (g_max + g_min);
            break;
        }
        case approximation_mode::MINRANGE: {
            // sech² = 1-tanh² is maximized at 0, decreasing in |x|; min slope at largest |endpoint|
            auto da = T(1) - fa * fa, db = T(1) - fb * fb;
            alpha = std::min(da, db);
            auto ya = fa - alpha * a, yb = fb - alpha * b;
            delta = T(0.5) * std::fabs(yb - ya);
            dzeta = T(0.5) * (ya + yb);
            break;
        }
        default:
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) - fa * fa;
            dzeta = fa - alpha * a; delta = T(0);
        }
        return apply_unary(c, alpha, dzeta, delta);
    }

    affine_form asin() const
    {
        if (terms_.empty()) {
            if (center_ < T(-1) || center_ > T(1))
                throw std::domain_error("affine_form::asin: argument out of [-1, 1]");
            return affine_form(context(), std::asin(center_));
        }
        auto a = min(), b = max();
        if (a < T(-1) || b > T(1))
            throw std::domain_error("affine_form::asin: interval not contained in [-1, 1]");
        auto c = center(), r = radius();
        auto fa = std::asin(a), fb = std::asin(b);
        T alpha, dzeta, delta;
        switch (context().approximation_mode()) {
        case approximation_mode::CHEBYSHEV: {
            // f'(x) = 1/sqrt(1-x²); tangent: x* = ±sqrt(1 - 1/alpha²) (alpha ≥ 1)
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) / std::sqrt(T(1) - c * c);
            auto g   = [&](T x) { return std::asin(x) - alpha * x; };
            auto g_max = std::max(g(a), g(b)), g_min = std::min(g(a), g(b));
            if (alpha > T(1)) {
                auto xs_pos = std::sqrt(T(1) - T(1) / (alpha * alpha));
                for (T xs : {xs_pos, -xs_pos}) {
                    if (xs > a && xs < b) { auto gx = g(xs); g_max = std::max(g_max, gx); g_min = std::min(g_min, gx); }
                }
            }
            delta = T(0.5) * (g_max - g_min);
            dzeta = T(0.5) * (g_max + g_min);
            break;
        }
        case approximation_mode::MINRANGE: {
            // f'(x) = 1/sqrt(1-x²); min slope at endpoint with larger |x|
            auto da = T(1) / std::sqrt(T(1) - a * a), db = T(1) / std::sqrt(T(1) - b * b);
            alpha = std::min(da, db);
            auto ya = fa - alpha * a, yb = fb - alpha * b;
            delta = T(0.5) * std::fabs(yb - ya);
            dzeta = T(0.5) * (ya + yb);
            break;
        }
        default:
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) / std::sqrt(T(1) - c * c);
            dzeta = fa - alpha * a; delta = T(0);
        }
        return apply_unary(c, alpha, dzeta, delta);
    }

    affine_form acos() const
    {
        if (terms_.empty()) {
            if (center_ < T(-1) || center_ > T(1))
                throw std::domain_error("affine_form::acos: argument out of [-1, 1]");
            return affine_form(context(), std::acos(center_));
        }
        auto a = min(), b = max();
        if (a < T(-1) || b > T(1))
            throw std::domain_error("affine_form::acos: interval not contained in [-1, 1]");
        auto c = center(), r = radius();
        auto fa = std::acos(a), fb = std::acos(b);
        T alpha, dzeta, delta;
        switch (context().approximation_mode()) {
        case approximation_mode::CHEBYSHEV: {
            // f'(x) = -1/sqrt(1-x²); tangent: x* = ±sqrt(1 - 1/alpha²) (alpha ≤ -1)
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(-1) / std::sqrt(T(1) - c * c);
            auto g   = [&](T x) { return std::acos(x) - alpha * x; };
            auto g_max = std::max(g(a), g(b)), g_min = std::min(g(a), g(b));
            if (alpha <= T(-1)) {
                auto xs_pos = std::sqrt(T(1) - T(1) / (alpha * alpha));
                for (T xs : {xs_pos, -xs_pos}) {
                    if (xs > a && xs < b) { auto gx = g(xs); g_max = std::max(g_max, gx); g_min = std::min(g_min, gx); }
                }
            }
            delta = T(0.5) * (g_max - g_min);
            dzeta = T(0.5) * (g_max + g_min);
            break;
        }
        case approximation_mode::MINRANGE: {
            // f'(x) = -1/sqrt(1-x²) ≤ 0; min slope (most negative) at endpoint with larger |x|
            auto da = T(-1) / std::sqrt(T(1) - a * a), db = T(-1) / std::sqrt(T(1) - b * b);
            alpha = std::min(da, db);
            auto ya = fa - alpha * a, yb = fb - alpha * b;
            delta = T(0.5) * std::fabs(yb - ya);
            dzeta = T(0.5) * (ya + yb);
            break;
        }
        default:
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(-1) / std::sqrt(T(1) - c * c);
            dzeta = fa - alpha * a; delta = T(0);
        }
        return apply_unary(c, alpha, dzeta, delta);
    }

    affine_form atan() const
    {
        if (terms_.empty())
            return affine_form(context(), std::atan(center_));
        auto c = center(), r = radius(), a = c - r, b = c + r;
        auto fa = std::atan(a), fb = std::atan(b);
        T alpha, dzeta, delta;
        switch (context().approximation_mode()) {
        case approximation_mode::CHEBYSHEV: {
            // f'(x) = 1/(1+x²); tangent: x* = ±sqrt(1/alpha - 1) (0 < alpha ≤ 1)
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) / (T(1) + c * c);
            auto g   = [&](T x) { return std::atan(x) - alpha * x; };
            auto g_max = std::max(g(a), g(b)), g_min = std::min(g(a), g(b));
            if (alpha > T(0) && alpha <= T(1)) {
                auto xs_pos = std::sqrt(T(1) / alpha - T(1));
                for (T xs : {xs_pos, -xs_pos}) {
                    if (xs > a && xs < b) { auto gx = g(xs); g_max = std::max(g_max, gx); g_min = std::min(g_min, gx); }
                }
            }
            delta = T(0.5) * (g_max - g_min);
            dzeta = T(0.5) * (g_max + g_min);
            break;
        }
        case approximation_mode::MINRANGE: {
            // f'(x) = 1/(1+x²) maximized at 0, decreasing in |x|; min slope at largest |endpoint|
            auto da = T(1) / (T(1) + a * a), db = T(1) / (T(1) + b * b);
            alpha = std::min(da, db);
            auto ya = fa - alpha * a, yb = fb - alpha * b;
            delta = T(0.5) * std::fabs(yb - ya);
            dzeta = T(0.5) * (ya + yb);
            break;
        }
        default:
            alpha = r > limits<T>::minrad ? (fb - fa) / (b - a) : T(1) / (T(1) + c * c);
            dzeta = fa - alpha * a; delta = T(0);
        }
        return apply_unary(c, alpha, dzeta, delta);
    }

    // ---------------------------------------------------------------------------
    // accuracy improvements
    // ---------------------------------------------------------------------------

    // Split into two independent sub-forms covering [min,mid] and [mid,max].
    std::pair<affine_form, affine_form> bisect() const
    {
        auto [lo, hi] = to_interval().split();
        return { affine_form(context(), lo), affine_form(context(), hi) };
    }

    // Condense to at most max_terms noise symbols by merging the smallest ones.
    // Preserves interval bounds (conservative).
    affine_form& condense(size_t max_terms)
    {
        if (terms_.size() <= max_terms) return *this;

        size_t n_merge = terms_.size() - max_terms + 1;

        std::nth_element(terms_.begin(), terms_.begin() + n_merge, terms_.end(),
            [](term const& a, term const& b){ return std::fabs(a.value) < std::fabs(b.value); });

        T error_sum = std::transform_reduce(terms_.begin(), terms_.begin() + n_merge,
            T(0), std::plus<T>{}, [](term const& t){ return std::fabs(t.value); });

        terms_.erase(terms_.begin(), terms_.begin() + n_merge);

        std::sort(terms_.begin(), terms_.end(),
            [](term const& a, term const& b){ return a.index < b.index; });

        // New error symbol has a fresh (largest) index — push_back preserves sort order.
        terms_.push_back({context().increment_last(), error_sum});

        update_radius();
        return *this;
    }

    // friends
    friend affine_form operator+(T v, affine_form const& af) { return af + v; }
    friend affine_form operator-(T v, affine_form const& af) { return -af + v; }
    friend affine_form operator*(T v, affine_form const& af) { return af * v; }
    friend affine_form operator/(T v, affine_form const& af) { return af.inv() * v; }

    friend std::ostream& operator<<(std::ostream& s, affine_form const& af)
    {
        s << "-------------------\n"
          << "center:     " << af.center_ << "\n"
          << "radius:     " << af.radius_ << "\n"
          << "terms:      [";
        for (size_t i = 0; i < af.terms_.size(); ++i) {
            if (i) s << ", ";
            s << "(" << af.terms_[i].index << ", " << af.terms_[i].value << ")";
        }
        s << "]\n"
          << "-------------------\n";
        return s;
    }

private:
    affine_context context_;
    T       center_;
    T       radius_;
    terms_t terms_;

    affine_form(affine_context const& ctx, T center, terms_t&& terms)
        : context_(ctx)
        , center_(center)
        , radius_(0)
        , terms_(std::move(terms))
    {
        update_radius();
    }

    affine_form apply_unary(T c, T alpha, T dzeta, T delta) const
    {
        auto result_terms = terms_;
        for (auto& t : result_terms) t.value *= alpha;
        result_terms.push_back({context().increment_last(), delta});
        return affine_form(context(), alpha * c + dzeta, std::move(result_terms));
    }

    void update_radius()
    {
        prune_zero_terms();
        radius_ = std::transform_reduce(terms_.begin(), terms_.end(), T(0),
            std::plus<T>{}, [](term const& t){ return std::fabs(t.value); });
    }

    void prune_zero_terms()
    {
        terms_.erase(std::remove_if(terms_.begin(), terms_.end(), [](term const& t) {
            return t.value == T(0);
        }), terms_.end());
    }

    void ensure_same_context(affine_form const& other) const
    {
        if (!context_.shares_state_with(other.context_))
            throw std::invalid_argument("affine_form: cannot combine forms from different affine_context instances");
    }

    template<typename LoFn, typename RoFn, typename BothFn>
    static terms_t merge_terms(terms_t const& lhs, terms_t const& rhs,
                                LoFn lo_fn, RoFn ro_fn, BothFn both_fn)
    {
        terms_t result;
        result.reserve(lhs.size() + rhs.size());
        auto li = lhs.begin(), ri = rhs.begin();
        while (li != lhs.end() && ri != rhs.end()) {
            if (li->index < ri->index) {
                result.push_back({li->index, lo_fn(*li)});
                ++li;
            } else if (ri->index < li->index) {
                result.push_back({ri->index, ro_fn(*ri)});
                ++ri;
            } else {
                result.push_back({li->index, both_fn(*li, *ri)});
                ++li; ++ri;
            }
        }
        for (; li != lhs.end(); ++li) result.push_back({li->index, lo_fn(*li)});
        for (; ri != rhs.end(); ++ri) result.push_back({ri->index, ro_fn(*ri)});
        return result;
    }
};

// ---------------------------------------------------------------------------
// evaluate_bisected: recursive bisection evaluator
// ---------------------------------------------------------------------------

template<typename F, std::floating_point T>
interval<T> evaluate_bisected(F&& f, affine_form<T> x, int depth = 1)
{
    if (depth <= 0)
        return f(x).to_interval();
    auto [lo, hi] = x.bisect();
    return evaluate_bisected(f, std::move(lo), depth - 1)
         | evaluate_bisected(f, std::move(hi), depth - 1);
}

} // namespace pappus

#endif
