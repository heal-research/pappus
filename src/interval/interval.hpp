#ifndef PAPPUS_INTERVAL_HPP
#define PAPPUS_INTERVAL_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <limits>
#include <ostream>
#include <utility>

#include <eve/module/core.hpp>
#include <eve/module/math.hpp>

#include "fp/util.hpp"
#include "fp/math.hpp"
#include "util/contracts.hpp"

namespace pappus {

// iterator and subdivision are scalar-only
template<std::floating_point T> class iterator;
template<std::floating_point T> class subdivision;

enum interval_class { M, Z, P, P0, P1, N, N0, N1, NONE };

template<typename T = double>
class interval {
public:
    interval()
        : interval(T(0), T(0))
    {
    }

    explicit interval(T v) : interval(v, v)
    {
    }

    explicit interval(T l, T h)
        : interval(check_bounds(l, h))
    {
    }

    explicit interval(std::string const& s) requires std::floating_point<T>
        : interval(fp::from_string<FE_DOWNWARD, T>(s), fp::from_string<FE_UPWARD, T>(s))
    {
    }

    explicit interval(std::string const& l, std::string const& u) requires std::floating_point<T>
        : interval(fp::from_string<FE_DOWNWARD, T>(l), fp::from_string<FE_UPWARD, T>(u))
    {
    }

    interval& operator=(const interval& other)
    {
        if (this != &other) {
            inf_ = other.inf_;
            sup_ = other.sup_;
        }
        return *this;
    }

    T inf() const { return inf_; }
    void set_lower(T value) { inf_ = value; }

    T sup() const { return sup_; }
    void set_upper(T value) { sup_ = value; }

    void set_bounds(T lower, T upper)
    {
        auto [lo, up] = check_bounds(lower, upper);
        inf_ = lo;
        sup_ = up;
    }

    void set_bounds(std::pair<T, T> bounds)
    {
        set_bounds(bounds.first, bounds.second);
    }

    std::pair<T, T> bounds() const
    {
        return {inf(), sup()};
    }

    T mid() const requires std::floating_point<T>
    {
        using num = std::numeric_limits<T>;
        if (is_empty())    return num::quiet_NaN();
        if (is_infinite()) return T(0);
        if (std::isinf(inf())) return num::lowest();
        if (std::isinf(sup())) return num::max();
        return fp::ropd<fp::op_mul>(inf(), T(0.5)) + fp::ropu<fp::op_mul>(sup(), T(0.5));
    }

    T radius() const requires std::floating_point<T>
    {
        if (is_empty()) return std::numeric_limits<T>::quiet_NaN();
        T m = mid();
        return std::fmax(fp::ropd<fp::op_sub>(m, inf()), fp::ropu<fp::op_sub>(sup(), m));
    }

    T diameter() const requires std::floating_point<T>
    {
        if (is_empty()) return std::numeric_limits<T>::quiet_NaN();
        return fp::ropu<fp::op_sub>(sup(), inf());
    }

    T mig() const requires std::floating_point<T>
    {
        if (is_empty())     return fp::nan_v<T>;
        if (contains(T(0))) return T(0);
        return std::min(fp::ropd<fp::op_abs>(inf()), fp::ropd<fp::op_abs>(sup()));
    }

    T mag() const requires std::floating_point<T>
    {
        if (is_empty()) return fp::nan_v<T>;
        return std::max(fp::ropd<fp::op_abs>(inf()), fp::ropd<fp::op_abs>(sup()));
    }

    bool contains(T v) const requires std::floating_point<T>
    {
        return inf() <= v && v <= sup();
    }

    bool contains_strict(T v) const requires std::floating_point<T>
    {
        return inf() < v && v < sup();
    }

    bool contains(interval const other) const requires std::floating_point<T>
    {
        return inf() <= other.inf() && other.sup() <= sup();
    }

    bool contains_strict(interval const other) const requires std::floating_point<T>
    {
        return inf() < other.inf() && other.sup() < sup();
    }

    bool is_finite() const requires std::floating_point<T>
    {
        return std::isfinite(inf()) && std::isfinite(sup());
    }

    bool is_infinite() const requires std::floating_point<T>
    {
        return std::isinf(inf()) && std::isinf(sup());
    }

    bool is_empty() const requires std::floating_point<T>
    {
        return !(inf() <= sup()); // negation handles NaNs
    }

    bool is_zero() const requires std::floating_point<T>
    {
        return inf() == T(0) && sup() == T(0);
    }

    bool is_symmetric() const requires std::floating_point<T>
    {
        return !is_empty() && (-inf() == sup());
    }

    bool is_positive() const requires std::floating_point<T>
    {
        return is_empty() || (inf() >= T(0));
    }

    bool is_strictly_positive() const requires std::floating_point<T>
    {
        return is_empty() || (inf() > T(0));
    }

    std::pair<interval, interval> split() const requires std::floating_point<T>
    {
        auto m = mid();
        auto left  = interval { inf(), m };
        auto right = interval { std::nextafter(m, sup()), sup() };
        return {left, right};
    }

    subdivision<T> split(size_t n) const requires std::floating_point<T>;

    // intersection — scalar only (SIMD branching on empty check is non-trivial)
    interval operator&(interval const other) const requires std::floating_point<T>
    {
        if (is_empty() || other.is_empty())
            return interval::empty();
        if (sup() < other.inf() || other.sup() < inf())
            return interval::empty();
        return interval(std::fmax(inf(), other.inf()), std::fmin(sup(), other.sup()));
    }

    interval& operator&=(interval const other) requires std::floating_point<T>
    {
        auto tmp = *this & other;
        std::swap(*this, tmp);
        return *this;
    }

    // hull — works for SIMD: eve::min/max ignore NaN (NaN = empty interval)
    interval operator|(interval const other) const
    {
        if constexpr (std::floating_point<T>) {
            if (is_empty())       return other;
            if (other.is_empty()) return *this;
            return interval(std::fmin(inf(), other.inf()), std::fmax(sup(), other.sup()));
        } else {
            return interval(std::pair<T,T>{
                eve::min(inf(), other.inf()),
                eve::max(sup(), other.sup())
            });
        }
    }

    interval& operator|=(interval const other)
    {
        auto tmp = *this | other;
        std::swap(*this, tmp);
        return *this;
    }

    bool operator==(interval const other) const requires std::floating_point<T>
    {
        return (is_empty() && other.is_empty()) ||
               (inf() == other.inf() && sup() == other.sup());
    }

    bool operator!=(interval const other) const requires std::floating_point<T>
    {
        return !(*this == other);
    }

    bool operator<(interval const other) const requires std::floating_point<T>
    {
        if (is_empty()) return !other.is_empty();
        return sup() < other.inf();
    }

    bool operator<=(interval const other) const requires std::floating_point<T>
    {
        if (is_empty()) return true;
        return sup() <= other.inf();
    }

    // arithmetic operators
    interval operator+() const { return *this; }

    interval operator+(T v) const
    {
        return interval(fp::ropd<fp::op_add>(inf(), v), fp::ropu<fp::op_add>(sup(), v));
    }

    interval operator+(interval const other) const
    {
        auto [a, b] = bounds();
        auto [c, d] = other.bounds();
        return interval(fp::ropd<fp::op_add>(a, c), fp::ropu<fp::op_add>(b, d));
    }

    interval operator-() const
    {
        return interval(-sup(), -inf());
    }

    interval operator-(T v) const
    {
        return interval(fp::ropd<fp::op_sub>(inf(), v), fp::ropu<fp::op_sub>(sup(), v));
    }

    interval operator-(interval const other) const
    {
        auto [a, b] = bounds();
        auto [c, d] = other.bounds();
        return interval(fp::ropd<fp::op_sub>(a, d), fp::ropu<fp::op_sub>(b, c));
    }

    interval operator*(T v) const
    {
        if constexpr (std::floating_point<T>) {
            if (is_empty() || std::isnan(v)) return interval::empty();
            if (v == T(0))                   return interval::zero();
            return v < T(0)
                ? interval(fp::ropd<fp::op_mul>(sup(), v), fp::ropu<fp::op_mul>(inf(), v))
                : interval(fp::ropd<fp::op_mul>(inf(), v), fp::ropu<fp::op_mul>(sup(), v));
        } else {
            auto neg_mask = v < T(0);
            return interval(std::pair<T,T>{
                eve::if_else(neg_mask, fp::ropd<fp::op_mul>(sup(), v), fp::ropd<fp::op_mul>(inf(), v)),
                eve::if_else(neg_mask, fp::ropu<fp::op_mul>(inf(), v), fp::ropu<fp::op_mul>(sup(), v))
            });
        }
    }

    interval operator*(interval const other) const
    {
        if constexpr (std::floating_point<T>) {
            if (is_empty() || other.is_empty())       return interval::empty();
            if (is_zero()  || other.is_zero())        return interval::zero();
            if (is_infinite() || other.is_infinite()) return interval::infinite();

            auto [a, b] = bounds();
            auto [c, d] = other.bounds();

            const auto D = [](auto x, auto y) { return fp::ropd<fp::op_mul>(x, y); };
            const auto U = [](auto x, auto y) { return fp::ropu<fp::op_mul>(x, y); };

            if (b <= T(0)) { // N
                if (d <= T(0))                     return interval(D(b, d), U(a, c));
                if (other.contains_strict(T(0)))   return interval(D(a, d), U(a, c));
                if (T(0) <= c)                     return interval(D(a, d), U(b, c));
            } else if (contains_strict(T(0))) { // M
                if (d <= T(0))                     return interval(D(b, c), U(a, c));
                if (other.contains_strict(T(0)))   return interval(std::fmin(D(a, d), D(b, c)), std::fmax(U(a, c), U(b, d)));
                if (T(0) <= c)                     return interval(D(a, d), U(b, d));
            } else if (T(0) <= a) { // P
                if (d <= T(0))                     return interval(D(b, c), U(a, d));
                if (other.contains_strict(T(0)))   return interval(D(b, c), U(b, d));
                if (T(0) <= c)                     return interval(D(a, c), U(b, d));
            }
            return interval::empty();
        } else {
            // SIMD: compute all 4 corner products, take lane-wise min/max
            auto [a, b] = bounds();
            auto [c, d] = other.bounds();
            const auto D = [](T x, T y) { return fp::ropd<fp::op_mul>(x, y); };
            const auto U = [](T x, T y) { return fp::ropu<fp::op_mul>(x, y); };
            return interval(std::pair<T,T>{
                eve::min(eve::min(D(a,c), D(a,d)), eve::min(D(b,c), D(b,d))),
                eve::max(eve::max(U(a,c), U(a,d)), eve::max(U(b,c), U(b,d)))
            });
        }
    }

    interval operator/(T v) const
    {
        if constexpr (std::floating_point<T>) {
            if (is_empty() || v == T(0) || std::isnan(v)) return interval::empty();
            return v < T(0)
                ? interval(fp::ropd<fp::op_div>(sup(), v), fp::ropu<fp::op_div>(inf(), v))
                : interval(fp::ropd<fp::op_div>(inf(), v), fp::ropu<fp::op_div>(sup(), v));
        } else {
            auto neg_mask = v < T(0);
            return interval(std::pair<T,T>{
                eve::if_else(neg_mask, fp::ropd<fp::op_div>(sup(), v), fp::ropd<fp::op_div>(inf(), v)),
                eve::if_else(neg_mask, fp::ropu<fp::op_div>(inf(), v), fp::ropu<fp::op_div>(sup(), v))
            });
        }
    }

    interval operator/(interval const other) const
    {
        if constexpr (std::floating_point<T>) {
            if (is_empty() || other.is_empty() || other.is_zero()) return interval::empty();
            if (is_zero())                                          return interval::zero();
            if (other.contains_strict(T(0)))                       return interval::infinite();

            auto [a, b] = bounds();
            auto [c, d] = other.bounds();

            const auto D = [](auto x, auto y) { return fp::ropd<fp::op_div>(x, y); };
            const auto U = [](auto x, auto y) { return fp::ropu<fp::op_div>(x, y); };

            if (T(0) <= c) { // P
                if (T(0) < a)   return interval(D(a, d), U(b, c));
                if (T(0) == a)  return interval(T(0), U(b, c));
                if (T(0) < b)   return interval(D(a, c), U(b, c));
                if (T(0) == b)  return interval(D(a, c), T(0));
                if (b < T(0))   return interval(D(a, c), U(b, d));
            } else if (d <= T(0)) { // N
                if (T(0) < a)   return interval(D(b, d), U(a, c));
                if (T(0) == a)  return interval(D(b, d), T(0));
                if (T(0) < b)   return interval(D(b, d), U(a, d));
                if (T(0) == b)  return interval(T(0), U(a, d));
                if (b < T(0))   return interval(D(b, c), U(a, d));
            }
            return interval::empty();
        } else {
            // SIMD: compute all 4 quotients; override with ±inf when divisor contains zero
            auto [a, b] = bounds();
            auto [c, d] = other.bounds();
            const auto D = [](T x, T y) { return fp::ropd<fp::op_div>(x, y); };
            const auto U = [](T x, T y) { return fp::ropu<fp::op_div>(x, y); };
            auto lo = eve::min(eve::min(D(a,c), D(a,d)), eve::min(D(b,c), D(b,d)));
            auto hi = eve::max(eve::max(U(a,c), U(a,d)), eve::max(U(b,c), U(b,d)));
            auto zero_div = (c <= T(0)) && (T(0) <= d);
            auto iv = eve::inf(eve::as<T>{});
            return interval(std::pair<T,T>{
                eve::if_else(zero_div, -iv, lo),
                eve::if_else(zero_div,  iv, hi)
            });
        }
    }

    interval inv() const
    {
        if constexpr (std::floating_point<T>) {
            if (is_empty() || is_zero())   return interval::empty();
            if (contains_strict(T(0)))     return interval::infinite();
        }
        return interval(fp::ropd<fp::op_div>(T(1), sup()), fp::ropu<fp::op_div>(T(1), inf()));
    }

    interval exp() const
    {
        if constexpr (std::floating_point<T>) {
            return interval(std::fmax(T(0), fp::ropd<fp::op_exp>(inf())), fp::ropu<fp::op_exp>(sup()));
        } else {
            return interval(eve::max(T(0), fp::ropd<fp::op_exp>(inf())), fp::ropu<fp::op_exp>(sup()));
        }
    }

    interval log() const
    {
        if constexpr (std::floating_point<T>) {
            auto x = *this & interval(T(0), fp::inf_v<T>);
            if (x.is_empty() || x.sup() <= T(0)) return interval::empty();
            return interval(fp::ropd<fp::op_log>(x.inf()), fp::ropu<fp::op_log>(x.sup()));
        } else {
            // Clamp both endpoints to [0, +inf] before taking log
            auto lo = eve::max(inf(), T(0));
            auto hi = eve::max(sup(), T(0));
            return interval(fp::ropd<fp::op_log>(lo), fp::ropu<fp::op_log>(hi));
        }
    }

    interval sin() const requires std::floating_point<T>
    {
        if (is_empty() || is_zero()) return *this;
        if (is_infinite() || diameter() > fp::two_pi_v<T>) return interval(T(-1), T(1));

        auto [a, b] = bounds();
        int x = fp::trig::get_quadrant(a);
        int y = fp::trig::get_quadrant(b);

        const auto D = [](auto v) { return fp::ropd<fp::op_sin>(v); };
        const auto U = [](auto v) { return fp::ropu<fp::op_sin>(v); };

        if (x == y) {
            if (diameter() > fp::pi_v<T>) return interval(T(-1), T(1));
            if (x == 0 || x == 3) return interval(D(a), U(b));
            if (x == 1 || x == 2) return interval(D(b), U(a));
        } else {
            if (x == 3 && y == 0) return interval(D(a), U(b));
            if (x == 1 && y == 2) return interval(D(b), U(a));
            if ((x == 0 || x == 3) && (y == 1 || y == 2)) return interval(std::fmin(D(a), D(b)), T(1));
            if ((x == 1 || x == 2) && (y == 0 || y == 3)) return interval(T(-1), std::fmax(U(a), U(b)));
        }
        return interval(T(-1), T(1));
    }

    interval cos() const requires std::floating_point<T>
    {
        if (is_empty()) return *this;
        if (is_infinite() || diameter() > fp::two_pi_v<T>) return interval(T(-1), T(1));

        auto [a, b] = bounds();
        auto x = fp::trig::get_quadrant(a);
        auto y = fp::trig::get_quadrant(b);

        const auto D = [](auto v) { return fp::ropd<fp::op_cos>(v); };
        const auto U = [](auto v) { return fp::ropu<fp::op_cos>(v); };

        if (x == y) {
            if (diameter() > fp::pi_v<T>) return interval(T(-1), T(1));
            if (x == 0 || x == 1) return interval(D(b), U(a));
            if (x == 2 || x == 3) return interval(D(a), U(b));
        } else {
            if (x == 2 && y == 3) return interval(D(a), U(b));
            if (x == 0 && y == 1) return interval(D(b), U(a));
            if ((x == 2 || x == 3) && (y == 0 || y == 1)) return interval(std::fmin(D(a), D(b)), T(1));
            if ((x == 0 || x == 1) && (y == 2 || y == 3)) return interval(T(-1), std::fmax(U(a), U(b)));
        }
        return interval(T(-1), T(1));
    }

    interval tan() const requires std::floating_point<T>
    {
        if (is_empty()) return *this;
        if (diameter() > fp::pi_v<T>) return interval::infinite();

        auto [a, b] = bounds();
        auto x = fp::trig::get_quadrant(a);
        auto y = fp::trig::get_quadrant(b);

        if (x % 2 <= y % 2 && x < y) return interval::infinite();
        return interval(fp::ropd<fp::op_tan>(a), fp::ropu<fp::op_tan>(b));
    }

    interval asin() const requires std::floating_point<T>
    {
        auto t = *this & interval(T(-1), T(1));
        if (t.is_empty()) return interval::empty();
        return interval(fp::ropd<fp::op_asin>(t.inf()), fp::ropu<fp::op_asin>(t.sup()));
    }

    interval acos() const requires std::floating_point<T>
    {
        auto t = *this & interval(T(-1), T(1));
        if (t.is_empty()) return interval::empty();
        return interval(fp::ropd<fp::op_acos>(t.sup()), fp::ropu<fp::op_acos>(t.inf()));
    }

    interval atan() const requires std::floating_point<T>
    {
        if (is_empty()) return interval::empty();
        return interval(fp::ropd<fp::op_atan>(inf()), fp::ropu<fp::op_atan>(sup()));
    }

    interval sinh() const requires std::floating_point<T>
    {
        if (is_empty()) return interval::empty();
        return interval(fp::ropd<fp::op_sinh>(inf()), fp::ropu<fp::op_sinh>(sup()));
    }

    interval cosh() const requires std::floating_point<T>
    {
        if (is_empty()) return interval::empty();
        return interval(fp::ropd<fp::op_cosh>(mig()), fp::ropu<fp::op_cosh>(mag()));
    }

    interval tanh() const requires std::floating_point<T>
    {
        if (is_empty()) return interval::empty();
        return interval(fp::ropd<fp::op_tanh>(inf()), fp::ropu<fp::op_tanh>(sup()));
    }

    interval square() const
    {
        if constexpr (std::floating_point<T>) {
            if (is_zero() || is_empty()) return *this;
            auto [a, b] = bounds();
            const auto D = [](auto x, auto y) { return fp::ropd<fp::op_mul>(x, y); };
            const auto U = [](auto x, auto y) { return fp::ropu<fp::op_mul>(x, y); };
            if (T(0) <= a) return interval(D(a, a), U(b, b));
            if (b <= T(0)) return interval(D(b, b), U(a, a));
            return interval(T(0), std::fmax(U(a, a), U(b, b)));
        } else {
            auto [a, b] = bounds();
            const auto D = [](T x, T y) { return fp::ropd<fp::op_mul>(x, y); };
            const auto U = [](T x, T y) { return fp::ropu<fp::op_mul>(x, y); };
            // Select correct bound for each lane based on sign
            auto lo = eve::if_else(a >= T(0), D(a,a),
                      eve::if_else(b <= T(0), D(b,b), T(0)));
            auto hi = eve::if_else(a >= T(0), U(b,b),
                      eve::if_else(b <= T(0), U(a,a),
                                   eve::max(U(a,a), U(b,b))));
            return interval(std::pair<T,T>{lo, hi});
        }
    }

    interval pow(int p) const requires std::floating_point<T>
    {
        if (is_empty()) return interval::empty();
        if (p == 0) return is_zero() ? interval::empty() : interval(T(1));
        if (p == 1) return *this;
        if (p < 0 && is_zero()) return interval::empty();

        auto [a, b] = bounds();

        if (p % 2 == 0) { // even power
            if (p > 0) {
                if (a >= T(0)) return interval(ipow_d(a, p), ipow_u(b, p));
                if (b <= T(0)) return interval(ipow_d(b, p), ipow_u(a, p));
                return interval(ipow_d(mig(), p), ipow_u(mag(), p));
            } else {
                if (a >= T(0)) return interval(ipow_d(b, p), ipow_u(a, p));
                if (b <= T(0)) return interval(ipow_d(a, p), ipow_u(b, p));
                return interval(ipow_d(mag(), p), ipow_u(mig(), p));
            }
        } else { // odd power
            if (is_infinite()) return interval::infinite();
            if (p > 0) {
                if (a == T(0)) return interval(T(0), ipow_u(b, p));
                if (b == T(0)) return interval(ipow_d(a, p), T(0));
                return interval(ipow_d(a, p), ipow_u(b, p));
            } else {
                if (a == T(0)) return interval(ipow_d(b, p), fp::inf_v<T>);
                if (b == T(0)) return interval(-fp::inf_v<T>, ipow_u(a, p));
                if (contains(T(0))) return interval::infinite();
                return interval(ipow_d(b, p), ipow_u(a, p));
            }
        }
    }

    interval pow(T p) const requires std::floating_point<T>
    {
        if (is_empty()) return interval::empty();
        if (std::fmod(p, T(1)) == T(0)) return this->pow(static_cast<int>(p));
        if (is_zero()) return p > T(0) ? interval::zero() : interval::empty();
        if (p == T(0.5)) return this->sqrt();
        return (p * this->log()).exp();
    }

    interval pow(interval const other) const requires std::floating_point<T>
    {
        auto x = (*this) & interval(T(0), fp::inf_v<T>);
        if (x.is_empty() || other.is_empty()) return interval::empty();
        return x.pow(other.inf()) | x.pow(other.sup());
    }

    interval sqrt() const
    {
        if constexpr (std::floating_point<T>) {
            auto x = *this & interval(T(0), fp::inf_v<T>);
            if (x.is_empty()) return interval::empty();
            return interval(fp::ropd<fp::op_sqrt>(x.inf()), fp::ropu<fp::op_sqrt>(x.sup()));
        } else {
            auto lo = eve::max(inf(), T(0));
            auto hi = eve::max(sup(), T(0));
            return interval(fp::ropd<fp::op_sqrt>(lo), fp::ropu<fp::op_sqrt>(hi));
        }
    }

    interval abs() const
    {
        if constexpr (std::floating_point<T>) {
            if (is_empty()) return interval::empty();
            auto a = inf(), b = sup();
            auto alo = std::fabs(a), ahi = std::fabs(b);
            if (a <= T(0) && b >= T(0)) return interval(T(0), std::fmax(alo, ahi));
            return interval(std::fmin(alo, ahi), std::fmax(alo, ahi));
        } else {
            auto lo = eve::abs(inf()), hi = eve::abs(sup());
            auto contains_zero = inf() <= T(0) && sup() >= T(0);
            auto out_lo = eve::if_else(contains_zero, T(0), eve::min(lo, hi));
            auto out_hi = eve::max(lo, hi);
            return interval(std::pair<T,T>{out_lo, out_hi});
        }
    }

    interval cbrt() const
    {
        if constexpr (std::floating_point<T>) {
            if (is_empty()) return interval::empty();
            // cbrt is monotonically increasing on all of R.
            return interval(fp::ropd<fp::op_cbrt>(inf()), fp::ropu<fp::op_cbrt>(sup()));
        } else {
            return interval(fp::ropd<fp::op_cbrt>(inf()), fp::ropu<fp::op_cbrt>(sup()));
        }
    }

    interval log1p() const
    {
        if constexpr (std::floating_point<T>) {
            auto x = *this & interval(T(-1), fp::inf_v<T>);
            if (x.is_empty() || x.sup() <= T(-1)) return interval::empty();
            auto lo = eve::max(x.inf(), T(-1));
            return interval(fp::ropd<fp::op_log1p>(lo), fp::ropu<fp::op_log1p>(x.sup()));
        } else {
            auto lo = eve::max(inf(), T(-1));
            auto hi = eve::max(sup(), T(-1));
            return interval(fp::ropd<fp::op_log1p>(lo), fp::ropu<fp::op_log1p>(hi));
        }
    }

    interval floor() const
    {
        if constexpr (std::floating_point<T>) {
            if (is_empty()) return interval::empty();
            // floor is monotonically non-decreasing: floor([a,b]) = [floor(a), floor(b)].
            // floor/ceil produce exact integers; directed rounding is a no-op.
            return interval(fp::ropd<fp::op_floor>(inf()), fp::ropu<fp::op_floor>(sup()));
        } else {
            return interval(fp::ropd<fp::op_floor>(inf()), fp::ropu<fp::op_floor>(sup()));
        }
    }

    interval ceil() const
    {
        if constexpr (std::floating_point<T>) {
            if (is_empty()) return interval::empty();
            // ceil is monotonically non-decreasing: ceil([a,b]) = [ceil(a), ceil(b)].
            return interval(fp::ropd<fp::op_ceil>(inf()), fp::ropu<fp::op_ceil>(sup()));
        } else {
            return interval(fp::ropd<fp::op_ceil>(inf()), fp::ropu<fp::op_ceil>(sup()));
        }
    }

    interval& operator+=(interval const other)
    {
        auto tmp = *this + other;
        std::swap(*this, tmp);
        return *this;
    }

    interval& operator-=(interval const other)
    {
        auto tmp = *this - other;
        std::swap(*this, tmp);
        return *this;
    }

    interval& operator*=(interval const other)
    {
        auto tmp = *this * other;
        std::swap(*this, tmp);
        return *this;
    }

    interval& operator/=(interval const other)
    {
        auto tmp = *this / other;
        std::swap(*this, tmp);
        return *this;
    }

    // friends
    friend interval operator+(T v, interval const i) { return i + v; }
    friend interval operator-(T v, interval const i) { return -i + v; }
    friend interval operator*(T v, interval const i) { return i * v; }
    friend interval operator/(T v, interval const i) { return i.inv() * v; }
    friend interval pow(T v, interval const i)       { return interval(v).pow(i); }

    template <interval_class C>
    bool is() const requires std::floating_point<T>
    {
        static_assert(C >= interval_class::M && C <= interval_class::N1, "Unknown interval class.");
        auto [a, b] = bounds();
        switch (C) {
        case M:  return a < T(0) && b > T(0);
        case Z:  return a == T(0) && b == T(0);
        case P:  return a >= T(0) && b > T(0);
        case P0: return a == T(0) && b > T(0);
        case P1: return a > T(0) && b > T(0);
        case N:  return a < T(0) && b <= T(0);
        case N0: return a < T(0) && b == T(0);
        case N1: return a < T(0) && b < T(0);
        default: return false;
        }
    }

    friend std::ostream& operator<<(std::ostream& s, interval const iv)
        requires std::floating_point<T>
    {
        s << "[" << iv.inf() << ", " << iv.sup() << "]";
        return s;
    }

    // factory constants
    static interval empty()
    {
        if constexpr (std::floating_point<T>)
            return interval(+fp::nan_v<T>, -fp::nan_v<T>);
        else {
            auto nan_val = eve::nan(eve::as<T>{});
            return interval(std::pair<T,T>{nan_val, nan_val});
        }
    }

    static interval zero() { return interval(T(+0.0), T(-0.0)); }

    static interval infinite()
    {
        if constexpr (std::floating_point<T>)
            return interval(-fp::inf_v<T>, +fp::inf_v<T>);
        else {
            auto iv = eve::inf(eve::as<T>{});
            return interval(std::pair<T,T>{-iv, +iv});
        }
    }

    static interval force_interval(T a, T b)
    {
        if (a > b) std::swap(a, b);
        return interval(a, b);
    }

    interval segment(std::size_t i, std::size_t n) const requires std::floating_point<T>
    {
        EXPECT(i < n);
        auto h = diameter() / T(n);
        auto a = fp::ropd<fp::op_add>(inf(), T(i) * h);
        auto b = fp::ropu<fp::op_add>(inf(), T(i + 1) * h);
        return interval(a, b);
    }

private:
    T inf_;
    T sup_;

    // Integer power rounded down/up via repeated squaring using op_mul.
    // Avoids std::pow(x,n) which uses exp(n*log(x)) and gives non-exact
    // results under directed rounding even for small integer arguments.
    static T ipow_d(T x, int p)
    {
        if (p < 0) return fp::ropd<fp::op_div>(T(1), ipow_u(x, -p));
        T r(1);
        for (int n = p; n > 0; n >>= 1) {
            if (n & 1) r = fp::ropd<fp::op_mul>(r, x);
            if (n > 1) x = fp::ropd<fp::op_mul>(x, x);
        }
        return r;
    }
    static T ipow_u(T x, int p)
    {
        if (p < 0) return fp::ropu<fp::op_div>(T(1), ipow_d(x, -p));
        T r(1);
        for (int n = p; n > 0; n >>= 1) {
            if (n & 1) r = fp::ropu<fp::op_mul>(r, x);
            if (n > 1) x = fp::ropu<fp::op_mul>(x, x);
        }
        return r;
    }

    static std::pair<T, T> check_bounds(T lo, T hi)
    {
        if constexpr (std::floating_point<T>) {
            // signed zero convention: +0 for left endpoints, -0 for right endpoints
            // see T. Hickey, Q. Ju, M. H. van Emden, Interval Arithmetic — Section 5.2
            EXPECT(!(lo > hi));
            if (lo == T(0)) lo = T(+0.0);
            if (hi == T(0)) hi = T(-0.0);
        } else {
            lo = eve::if_else(lo == T(0), T(+0.0), lo);
            hi = eve::if_else(hi == T(0), T(-0.0), hi);
        }
        return {lo, hi};
    }

    explicit interval(std::pair<T, T> iv)
        : inf_(iv.first)
        , sup_(iv.second)
    {
    }
};

// ---------------------------------------------------------------------------
// iterator<T> and subdivision<T> — scalar only
// ---------------------------------------------------------------------------

template<std::floating_point T>
class iterator {
public:
    using value_type        = interval<T>;
    using reference         = value_type const&;
    using pointer           = value_type const*;
    using difference_type   = std::size_t;
    using iterator_category = std::forward_iterator_tag;

    iterator(interval<T> const& iv, size_t n)
        : iv_(iv), i_(0), n_(n)
        , h_(iv_.diameter() / T(n_))
        , a_(iv_.inf())
        , b_(fp::ropu<fp::op_add>(a_, h_))
    {
        EXPECT(n_ > 0);
    }

    iterator& operator++()
    {
        ++i_;
        a_ = b_;
        b_ = (i_ == n_) ? iv_.sup() : fp::ropu<fp::op_add>(iv_.inf(), T(i_ + 1) * h_);
        return *this;
    }

    iterator operator++(int)
    {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    interval<T> operator*() const { return interval<T>(a_, b_); }

    bool operator==(iterator other) const
    {
        if (iv_ != other.iv_) return false;
        return i_ == other.i_;
    }

    bool operator!=(iterator other) const { return !(*this == other); }

    bool operator<(iterator other) const
    {
        if (iv_ != other.iv_) return false;
        return i_ < other.i_;
    }

    iterator begin() { return iterator(iv_, n_); }

    iterator end()
    {
        auto tmp = *this;
        tmp.advance_to_end();
        return tmp;
    }

private:
    void advance_to_end()
    {
        i_ = n_;
        a_ = fp::nan_v<T>;
        b_ = fp::nan_v<T>;
    }

    interval<T> const& iv_;
    std::size_t i_;
    std::size_t const n_;
    T const h_;
    T a_;
    T b_;
};

template<std::floating_point T>
class subdivision {
public:
    subdivision(interval<T> const& iv, std::size_t n) : iv_(iv), n_(n) {}

    iterator<T> begin() const { return iterator<T>(iv_, n_).begin(); }
    iterator<T> end()   const { return iterator<T>(iv_, n_).end(); }

private:
    interval<T> const& iv_;
    std::size_t const n_;
};

// split(n) definition — after subdivision<T> is complete
template<typename T>
subdivision<T> interval<T>::split(size_t n) const requires std::floating_point<T>
{
    EXPECT(n > 0);
    return subdivision<T>(*this, n);
}

// ---------------------------------------------------------------------------
// batch_evaluate_ia: SIMD interval evaluation via sub-interval packing
//
// Divides x into n_leaves sub-intervals, packs eve::cardinal_v<wide<T>> of
// them at a time into interval<wide<T>>, evaluates f on the batch, then
// reduces each SIMD result to a scalar hull contribution.
//
// f must be a generic callable: auto f(auto x) { ... } or a template.
// ---------------------------------------------------------------------------
template<typename F, std::floating_point T>
interval<T> batch_evaluate_ia(F&& f, interval<T> x, int n_leaves)
{
    using W = eve::wide<T>;
    constexpr int W_size = static_cast<int>(eve::cardinal_v<W>);

    EXPECT(n_leaves > 0);

    auto result = interval<T>::empty();

    alignas(W) std::array<T, W_size> lo_arr, hi_arr;
    int k = 0;
    for (; k + W_size <= n_leaves; k += W_size) {
        for (int i = 0; i < W_size; ++i) {
            auto seg = x.segment(static_cast<std::size_t>(k + i), static_cast<std::size_t>(n_leaves));
            lo_arr[i] = seg.inf();
            hi_arr[i] = seg.sup();
        }
        interval<W> sub(W(lo_arr.data()), W(hi_arr.data()));
        auto r = f(sub);
        result |= interval<T>(eve::minimum(r.inf()), eve::maximum(r.sup()));
    }
    // scalar tail for remaining sub-intervals
    for (; k < n_leaves; ++k)
        result |= f(x.segment(k, n_leaves));

    return result;
}

} // namespace pappus
#endif
