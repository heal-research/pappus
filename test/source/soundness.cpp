//
// Directed-rounding soundness suite (float, double, and SIMD lanes).
//
// Reference: MPFR at 256 bits. Every input endpoint is a float/double and is
// converted to MPFR exactly, so the reference value R is the true real result
// (up to 2^-256 relative, far below any 1-ULP effect). A bound is "sound" iff
// lo <= R <= hi in exact arithmetic; no tolerance of any kind is applied.
//
// For every operation and every interval shape (point, 1-ULP, 3-ULP, narrow,
// wide) the result must either be the empty interval (only when the input
// has no point in the operation's domain) or satisfy lo <= hi and enclose R
// at BOTH endpoints and at a few interior points of the input interval.
//
// Named regression tests cover the captured operon failures:
//   - cbrt of a 1-ULP-wide float interval produced lo > hi and aborted;
//   - eve::div[lower] returned bounds that were not lower bounds;
//   - transcendental bounds that were 1 ULP wide although eve's kernels err
//     by up to ~3.6 ULP.
//

#include <catch2/catch_all.hpp>

#include <mpfr.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

#include <eve/eve.hpp>
#include <eve/module/core.hpp>
#include <eve/module/math.hpp>

#include "pappus/pappus.hpp"

namespace {

using rng_t = std::mt19937_64;
constexpr double kInf = std::numeric_limits<double>::infinity();

// Intervals per (operation, shape) and per type, and the RNG seed, can be
// raised from the environment for a longer soak run:
//   PAPPUS_SOUNDNESS_PER_SHAPE=5000 PAPPUS_SOUNDNESS_SEED=7 paptest "[soundness]"
// The defaults are deterministic. Runtime is dominated by MPFR transcendental
// evaluations (~6 per interval).
std::uint64_t env_u64(char const* name, std::uint64_t dflt)
{
    char const* v = std::getenv(name);
    return v ? std::strtoull(v, nullptr, 10) : dflt;
}

int per_shape()
{
    static int const n = static_cast<int>(env_u64("PAPPUS_SOUNDNESS_PER_SHAPE", 120));
    return n;
}

rng_t make_rng(std::uint64_t tag)
{
    static std::uint64_t const seed = env_u64("PAPPUS_SOUNDNESS_SEED", 0);
    return rng_t(tag ^ (seed * 0x9E3779B97F4A7C15ULL));
}

enum Shape { POINT, ULP1, ULP3, NARROW, WIDE, N_SHAPES };
char const* const shape_name[N_SHAPES] = { "point", "1-ulp", "3-ulp", "narrow", "wide" };

// ---------------------------------------------------------------------------
// MPFR reference
// ---------------------------------------------------------------------------

struct Scratch {
    mpfr_t x, y, r, a, b;
    Scratch() { mpfr_inits2(256, x, y, r, a, b, static_cast<mpfr_ptr>(nullptr)); }
    ~Scratch() { mpfr_clears(x, y, r, a, b, static_cast<mpfr_ptr>(nullptr)); }
    Scratch(Scratch const&) = delete;
    Scratch& operator=(Scratch const&) = delete;
};

Scratch& S()
{
    static Scratch s;
    return s;
}

// lo <= r <= hi in exact arithmetic (lo/hi may be +-inf; a NaN reference
// value, i.e. an undefined result, is treated as "nothing to check").
bool encloses(double lo, double hi, mpfr_srcptr r)
{
    if (mpfr_nan_p(r)) return true;
    mpfr_set_d(S().a, lo, MPFR_RNDN);
    mpfr_set_d(S().b, hi, MPFR_RNDN);
    return mpfr_lessequal_p(S().a, r) && mpfr_lessequal_p(r, S().b);
}

struct Tally {
    std::size_t checks = 0, bad = 0;
    std::string first;
    void fail(std::string msg)
    {
        if (bad++ == 0) first = std::move(msg);
    }
};

std::string hex(double v)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, "%a", v);
    return buf;
}

std::string iv_str(double lo, double hi) { return "[" + hex(lo) + ", " + hex(hi) + "]"; }

// ---------------------------------------------------------------------------
// Random generation
// ---------------------------------------------------------------------------

double uni(rng_t& g, double a, double b) { return std::uniform_real_distribution<double>(a, b)(g); }
double sgn(rng_t& g) { return (g() & 1) ? 1.0 : -1.0; }
double mag(rng_t& g, double e0, double e1) { return std::pow(10.0, uni(g, e0, e1)); }
double smag(rng_t& g, double e0, double e1) { return sgn(g) * mag(g, e0, e1); }
// exponent bound for float / double
double E(bool f, double fl, double db) { return f ? fl : db; }

double g_real(rng_t& g, bool f)
{
    switch (g() % 5) {
    case 0: return uni(g, -10, 10);
    case 1: return smag(g, -E(f, 15, 60), E(f, 15, 60));
    case 2: return uni(g, -1, 1);
    case 3: return g() % 8 == 0 ? 0.0 : uni(g, -100, 100);
    default: return smag(g, -3, 3);
    }
}

double g_small(rng_t& g, bool f)   // cube / 4th / 5th / 7th power: no overflow
{
    switch (g() % 4) {
    case 0: return uni(g, -5, 5);
    case 1: return smag(g, -E(f, 4, 30), E(f, 4, 30));
    case 2: return uni(g, -1, 1);
    default: return g() % 10 == 0 ? 0.0 : uni(g, -2, 2);
    }
}

double g_nonzero(rng_t& g, bool f)   // 1/x, x^-2, x^-3: no overflow
{
    switch (g() % 3) {
    case 0: return smag(g, -E(f, 8, 60), E(f, 8, 60));
    case 1: return smag(g, -2, 2);
    default: return uni(g, -10, 10) + (g() & 1 ? 1e-3 : -1e-3);
    }
}

double g_pos(rng_t& g, bool f)   // strictly positive
{
    switch (g() % 4) {
    case 0: return mag(g, -E(f, 35, 300), E(f, 35, 300));
    case 1: return 1.0 + uni(g, -1e-3, 1e-3);
    case 2: return uni(g, 1e-3, 10);
    default: return mag(g, -5, 5);
    }
}

double g_nonneg(rng_t& g, bool f) { return g() % 12 == 0 ? 0.0 : g_pos(g, f); }

double g_gt_m1(rng_t& g, bool f)   // x > -1
{
    switch (g() % 5) {
    case 0: return uni(g, -0.999, 3);
    case 1: return -1.0 + mag(g, -6, 0);
    case 2: return smag(g, -E(f, 30, 300), -1) ;
    case 3: return mag(g, -5, E(f, 30, 300));
    default: return uni(g, -1e-3, 1e-3);
    }
}

double g_trig(rng_t& g, bool)   // sin / cos: includes critical points
{
    constexpr double pi = 3.141592653589793238462643383279502884;
    switch (g() % 5) {
    case 0: return uni(g, -10, 10);
    case 1: return uni(g, -1000, 1000);
    case 2: return std::round(uni(g, -20, 20)) * (pi / 2) + uni(g, -1e-3, 1e-3);
    case 3: return std::round(uni(g, -20, 20)) * (pi / 2);   // nearest to a critical point
    default: return smag(g, -10, 0);
    }
}

double g_tan(rng_t& g, bool)   // includes the neighbourhood of poles
{
    constexpr double pi = 3.141592653589793238462643383279502884;
    switch (g() % 4) {
    case 0: return uni(g, -1.5, 1.5);
    case 1: return uni(g, -10, 10);
    case 2: return (std::round(uni(g, -6, 6)) + 0.5) * pi + uni(g, -1e-3, 1e-3);
    default: return smag(g, -10, 0);
    }
}

double g_unit(rng_t& g, bool)   // [-1, 1]
{
    switch (g() % 5) {
    case 0: return uni(g, -1, 1);
    case 1: return sgn(g) * (1.0 - mag(g, -6, -1));
    case 2: return g() & 1 ? 1.0 : -1.0;
    case 3: return 0.0;
    default: return uni(g, -0.5, 0.5);
    }
}

double g_atan(rng_t& g, bool f)
{
    switch (g() % 3) {
    case 0: return uni(g, -10, 10);
    case 1: return smag(g, -E(f, 30, 300), E(f, 30, 300));
    default: return smag(g, -10, 0);
    }
}

double g_exp(rng_t& g, bool f)   // reaches overflow and (gradual) underflow
{
    switch (g() % 4) {
    case 0: return uni(g, -10, 10);
    case 1: return f ? uni(g, -105, 95) : uni(g, -745.5, 711);
    case 2: return smag(g, -10, 0);
    default: return f ? uni(g, 85, 90) : uni(g, 705, 712);
    }
}

double g_hyp(rng_t& g, bool f)   // sinh / cosh: reaches overflow in both directions
{
    switch (g() % 4) {
    case 0: return uni(g, -10, 10);
    case 1: return f ? uni(g, -95, 95) : uni(g, -715, 715);
    case 2: return smag(g, -10, 0);
    default: return sgn(g) * (f ? uni(g, 85, 90) : uni(g, 705, 712));
    }
}

double g_tanh(rng_t& g, bool)
{
    switch (g() % 3) {
    case 0: return uni(g, -20, 20);
    case 1: return smag(g, -10, 1);
    default: return uni(g, -1, 1);
    }
}

double g_cbrt(rng_t& g, bool f)
{
    switch (g() % 4) {
    case 0: return smag(g, -E(f, 35, 300), E(f, 35, 300));
    case 1: return uni(g, -10, 10);
    case 2: return g() % 10 == 0 ? 0.0 : smag(g, -3, 3);
    default: return uni(g, -1000, 1000);
    }
}

double g_binary(rng_t& g, bool f)
{
    switch (g() % 5) {
    case 0: return uni(g, -10, 10);
    case 1: return smag(g, -E(f, 12, 60), E(f, 12, 60));
    case 2: return uni(g, -1, 1);
    case 3: return g() % 10 == 0 ? 0.0 : uni(g, -1e3, 1e3);
    default: return smag(g, -2, 2);
    }
}

// ---------------------------------------------------------------------------
// Unary operation table
// ---------------------------------------------------------------------------

enum class U { sqrt, abs, sqr, powi, powr, inv, exp, log, log1p, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, cbrt, floor, ceil };

struct UOp {
    char const* name;
    U id;
    int p;      // integer exponent (powi)
    double pr;  // real exponent (powr)
    double (*gen)(rng_t&, bool);
    bool (*dom)(double);                 // is x in the domain?
    bool (*overlap)(double, double);     // does [lo, hi] meet the domain?
    void (*ref)(mpfr_ptr, mpfr_srcptr, UOp const&);
};

bool dom_all(double) { return true; }
bool dom_nonneg(double x) { return x >= 0; }
bool dom_pos(double x) { return x > 0; }
bool dom_gt_m1(double x) { return x > -1; }
bool dom_unit(double x) { return x >= -1 && x <= 1; }
bool dom_nonzero(double x) { return x != 0; }

bool ov_all(double, double) { return true; }
bool ov_nonneg(double, double hi) { return hi >= 0; }
bool ov_pos(double, double hi) { return hi > 0; }
bool ov_gt_m1(double, double hi) { return hi > -1; }
bool ov_unit(double lo, double hi) { return hi >= -1 && lo <= 1; }
bool ov_nonzero(double lo, double hi) { return !(lo == 0 && hi == 0); }

std::vector<UOp> const& unary_ops()
{
    static std::vector<UOp> const ops = {
        { "sqrt",  U::sqrt,  0, 0, g_nonneg,   dom_nonneg,  ov_nonneg,  [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_sqrt(r, x, MPFR_RNDN); } },
        { "abs",   U::abs,   0, 0, g_real,     dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_abs(r, x, MPFR_RNDN); } },
        { "square", U::sqr,  0, 0, g_real,     dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_sqr(r, x, MPFR_RNDN); } },
        { "pow(2)", U::powi, 2, 0, g_real,     dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const& o) { mpfr_pow_si(r, x, o.p, MPFR_RNDN); } },
        { "pow(3)", U::powi, 3, 0, g_small,    dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const& o) { mpfr_pow_si(r, x, o.p, MPFR_RNDN); } },
        { "pow(4)", U::powi, 4, 0, g_small,    dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const& o) { mpfr_pow_si(r, x, o.p, MPFR_RNDN); } },
        { "pow(5)", U::powi, 5, 0, g_small,    dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const& o) { mpfr_pow_si(r, x, o.p, MPFR_RNDN); } },
        { "pow(7)", U::powi, 7, 0, g_small,    dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const& o) { mpfr_pow_si(r, x, o.p, MPFR_RNDN); } },
        { "pow(-1)", U::powi, -1, 0, g_nonzero, dom_nonzero, ov_nonzero, [](mpfr_ptr r, mpfr_srcptr x, UOp const& o) { mpfr_pow_si(r, x, o.p, MPFR_RNDN); } },
        { "pow(-2)", U::powi, -2, 0, g_nonzero, dom_nonzero, ov_nonzero, [](mpfr_ptr r, mpfr_srcptr x, UOp const& o) { mpfr_pow_si(r, x, o.p, MPFR_RNDN); } },
        { "pow(-3)", U::powi, -3, 0, g_nonzero, dom_nonzero, ov_nonzero, [](mpfr_ptr r, mpfr_srcptr x, UOp const& o) { mpfr_pow_si(r, x, o.p, MPFR_RNDN); } },
        { "pow(-5)", U::powi, -5, 0, g_nonzero, dom_nonzero, ov_nonzero, [](mpfr_ptr r, mpfr_srcptr x, UOp const& o) { mpfr_pow_si(r, x, o.p, MPFR_RNDN); } },
        { "inv",   U::inv,   0, 0, g_nonzero,  dom_nonzero, ov_nonzero, [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_ui_div(r, 1, x, MPFR_RNDN); } },
        { "exp",   U::exp,   0, 0, g_exp,      dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_exp(r, x, MPFR_RNDN); } },
        { "log",   U::log,   0, 0, g_pos,      dom_pos,     ov_pos,     [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_log(r, x, MPFR_RNDN); } },
        { "log1p", U::log1p, 0, 0, g_gt_m1,    dom_gt_m1,   ov_gt_m1,   [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_log1p(r, x, MPFR_RNDN); } },
        { "sin",   U::sin,   0, 0, g_trig,     dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_sin(r, x, MPFR_RNDN); } },
        { "cos",   U::cos,   0, 0, g_trig,     dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_cos(r, x, MPFR_RNDN); } },
        { "tan",   U::tan,   0, 0, g_tan,      dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_tan(r, x, MPFR_RNDN); } },
        { "asin",  U::asin,  0, 0, g_unit,     dom_unit,    ov_unit,    [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_asin(r, x, MPFR_RNDN); } },
        { "acos",  U::acos,  0, 0, g_unit,     dom_unit,    ov_unit,    [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_acos(r, x, MPFR_RNDN); } },
        { "atan",  U::atan,  0, 0, g_atan,     dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_atan(r, x, MPFR_RNDN); } },
        { "sinh",  U::sinh,  0, 0, g_hyp,      dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_sinh(r, x, MPFR_RNDN); } },
        { "cosh",  U::cosh,  0, 0, g_hyp,      dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_cosh(r, x, MPFR_RNDN); } },
        { "tanh",  U::tanh,  0, 0, g_tanh,     dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_tanh(r, x, MPFR_RNDN); } },
        { "cbrt",  U::cbrt,  0, 0, g_cbrt,     dom_all,     ov_all,     [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_cbrt(r, x, MPFR_RNDN); } },
        { "floor", U::floor, 0, 0, g_real, dom_all, ov_all, [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_floor(r, x); } },
        { "ceil",  U::ceil,  0, 0, g_real, dom_all, ov_all, [](mpfr_ptr r, mpfr_srcptr x, UOp const&) { mpfr_ceil(r, x); } },
        // pow(real p): x > 0
        { "pow(2.5)",  U::powr, 0, 2.5,  g_pos, dom_pos, ov_pos, [](mpfr_ptr r, mpfr_srcptr x, UOp const& o) { mpfr_set_d(S().y, o.pr, MPFR_RNDN); mpfr_pow(r, x, S().y, MPFR_RNDN); } },
        { "pow(-1.5)", U::powr, 0, -1.5, g_pos, dom_pos, ov_pos, [](mpfr_ptr r, mpfr_srcptr x, UOp const& o) { mpfr_set_d(S().y, o.pr, MPFR_RNDN); mpfr_pow(r, x, S().y, MPFR_RNDN); } },
        { "pow(0.375)", U::powr, 0, 0.375, g_pos, dom_pos, ov_pos, [](mpfr_ptr r, mpfr_srcptr x, UOp const& o) { mpfr_set_d(S().y, o.pr, MPFR_RNDN); mpfr_pow(r, x, S().y, MPFR_RNDN); } },
        { "pow(0.5)",  U::powr, 0, 0.5,  g_pos, dom_pos, ov_pos, [](mpfr_ptr r, mpfr_srcptr x, UOp const& o) { mpfr_set_d(S().y, o.pr, MPFR_RNDN); mpfr_pow(r, x, S().y, MPFR_RNDN); } },
    };
    return ops;
}

// pow(x, p) with a large |p log x| overflows; keep the real-exponent bases
// moderate so the interesting failures are rounding, not range, errors.
double clamp_pow_base(UOp const& op, double x, bool f)
{
    if (op.id != U::powr) return x;
    double const lim = f ? 1e6 : 1e30;
    return std::clamp(x, 1.0 / lim, lim);
}

template<class I>
I apply(UOp const& op, I const& v)
{
    using E = std::remove_cvref_t<decltype(v.inf())>;
    switch (op.id) {
    case U::sqrt:  return v.sqrt();
    case U::abs:   return v.abs();
    case U::sqr:   return v.square();
    case U::powi:  return v.pow(op.p);
    case U::powr:  return v.pow(E(static_cast<eve::element_type_t<E>>(op.pr)));
    case U::inv:   return v.inv();
    case U::exp:   return v.exp();
    case U::log:   return v.log();
    case U::log1p: return v.log1p();
    case U::sin:   return v.sin();
    case U::cos:   return v.cos();
    case U::tan:   return v.tan();
    case U::asin:  return v.asin();
    case U::acos:  return v.acos();
    case U::atan:  return v.atan();
    case U::sinh:  return v.sinh();
    case U::cosh:  return v.cosh();
    case U::tanh:  return v.tanh();
    case U::cbrt:  return v.cbrt();
    case U::floor: return v.floor();
    case U::ceil:  return v.ceil();
    }
    return v;
}

// ---------------------------------------------------------------------------
// Interval shapes and sample points
// ---------------------------------------------------------------------------

template<class T>
T clampT(double v, T lo, T hi)
{
    T t = static_cast<T>(v);
    return std::clamp(t, lo, hi);
}

template<class T>
std::pair<T, T> make_bounds(int shape, rng_t& g, auto&& gen)
{
    constexpr T inf = std::numeric_limits<T>::infinity();
    T x = static_cast<T>(gen());
    switch (shape) {
    case POINT: return { x, x };
    case ULP1:
        return (g() & 1) ? std::pair<T, T>{ x, std::nextafter(x, inf) }
                         : std::pair<T, T>{ std::nextafter(x, -inf), x };
    case ULP3: {
        T a = std::nextafter(x, -inf), b = std::nextafter(std::nextafter(x, inf), inf);
        return { a, b };
    }
    case NARROW: {
        double const r = uni(g, 1e-4, 1e-2) * (x == 0 ? 0.1 : std::fabs(static_cast<double>(x)));
        return { static_cast<T>(static_cast<double>(x) - r), static_cast<T>(static_cast<double>(x) + r) };
    }
    default: {
        T y = static_cast<T>(gen());
        return { std::min(x, y), std::max(x, y) };
    }
    }
}

template<class T>
void sample_points(T lo, T hi, rng_t& g, std::vector<T>& out)
{
    out.clear();
    out.push_back(lo);
    out.push_back(hi);
    if (lo == hi) return;
    double const a = static_cast<double>(lo), b = static_cast<double>(hi);
    out.push_back(clampT<T>(0.5 * a + 0.5 * b, lo, hi));
    for (int k = 0; k < 3; ++k) out.push_back(clampT<T>(uni(g, a, b), lo, hi));
    if ((a > 0 && b > 0) || (a < 0 && b < 0))   // geometric mean: interior points of wide positive/negative ranges
        out.push_back(clampT<T>((a < 0 ? -1.0 : 1.0) * std::sqrt(a * b), lo, hi));
}

// ---------------------------------------------------------------------------
// Result checks
// ---------------------------------------------------------------------------

template<class T>
// strict_domain: a domain error must give the empty interval. SIMD lanes are
// checked non-strictly: a lane may answer a domain error with an enclosure
// (e.g. x / [0,0] -> [-inf, inf]), which is sound; only inputs inside the
// domain are required to produce a non-empty enclosing result.
void check_unary(Tally& t, UOp const& op, T lo, T hi, T rlo, T rhi, rng_t& g, bool strict_domain = true)
{
    ++t.checks;
    double const dlo = lo, dhi = hi;
    bool const empty = !(rlo <= rhi);
    auto const where = [&] { return std::string(op.name) + " in=" + iv_str(dlo, dhi) + " out=" + iv_str(rlo, rhi); };

    if (!op.overlap(dlo, dhi)) {
        if (!empty && strict_domain) t.fail(where() + ": domain error must be empty");
        return;
    }
    if (empty) {
        t.fail(where() + ": empty result for an input meeting the domain");
        return;
    }
    std::vector<T> pts;
    sample_points(lo, hi, g, pts);
    for (T p : pts) {
        double const x = p;
        if (!op.dom(x)) continue;
        mpfr_set_d(S().x, x, MPFR_RNDN);
        op.ref(S().r, S().x, op);
        if (!encloses(rlo, rhi, S().r))
            t.fail(where() + " at x=" + hex(x) + " exact=" + hex(mpfr_get_d(S().r, MPFR_RNDN)) + " not enclosed");
    }
}

template<class T>
void run_unary(Tally& t, UOp const& op, int shape, rng_t& g, int n)
{
    using I = pappus::interval<T>;
    for (int i = 0; i < n; ++i) {
        auto [lo, hi] = make_bounds<T>(shape, g, [&] { return clamp_pow_base(op, op.gen(g, std::is_same_v<T, float>), std::is_same_v<T, float>); });
        // a 1-ULP / narrow shape around a boundary point may leave the finite range
        if (!std::isfinite(lo) || !std::isfinite(hi)) continue;
        I const res = apply(op, I(lo, hi));
        check_unary<T>(t, op, lo, hi, res.inf(), res.sup(), g);
    }
}

template<class W>
void run_unary_wide(Tally& t, UOp const& op, int shape, rng_t& g, int n)
{
    using T = eve::element_type_t<W>;
    constexpr std::size_t L = W::size();
    using IW = pappus::interval<W>;
    for (int i = 0; i < n; i += static_cast<int>(L)) {
        std::array<T, L> lo{}, hi{};
        W wlo{}, whi{};
        for (std::size_t k = 0; k < L; ++k) {
            auto [a, b] = make_bounds<T>(shape, g, [&] { return clamp_pow_base(op, op.gen(g, std::is_same_v<T, float>), std::is_same_v<T, float>); });
            if (!std::isfinite(a) || !std::isfinite(b)) { a = b = T(1); }
            lo[k] = a; hi[k] = b;
            wlo.set(k, a); whi.set(k, b);
        }
        IW const res = apply(op, IW(wlo, whi));
        for (std::size_t k = 0; k < L; ++k) {
            T const rlo = res.inf().get(k), rhi = res.sup().get(k);
            check_unary<T>(t, op, lo[k], hi[k], rlo, rhi, g, /*strict_domain=*/false);
        }
    }
}

// ---------------------------------------------------------------------------
// Binary operations
// ---------------------------------------------------------------------------

enum class B { add, sub, mul, div, add_s, sub_s, mul_s, div_s, powi_iv };
struct BOp {
    char const* name;
    B id;
    void (*ref)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr);
};

std::vector<BOp> const& binary_ops()
{
    static std::vector<BOp> const ops = {
        { "x + y",   B::add,   [](mpfr_ptr r, mpfr_srcptr x, mpfr_srcptr y) { mpfr_add(r, x, y, MPFR_RNDN); } },
        { "x - y",   B::sub,   [](mpfr_ptr r, mpfr_srcptr x, mpfr_srcptr y) { mpfr_sub(r, x, y, MPFR_RNDN); } },
        { "x * y",   B::mul,   [](mpfr_ptr r, mpfr_srcptr x, mpfr_srcptr y) { mpfr_mul(r, x, y, MPFR_RNDN); } },
        { "x / y",   B::div,   [](mpfr_ptr r, mpfr_srcptr x, mpfr_srcptr y) { mpfr_div(r, x, y, MPFR_RNDN); } },
        { "x + T",   B::add_s, [](mpfr_ptr r, mpfr_srcptr x, mpfr_srcptr y) { mpfr_add(r, x, y, MPFR_RNDN); } },
        { "x - T",   B::sub_s, [](mpfr_ptr r, mpfr_srcptr x, mpfr_srcptr y) { mpfr_sub(r, x, y, MPFR_RNDN); } },
        { "x * T",   B::mul_s, [](mpfr_ptr r, mpfr_srcptr x, mpfr_srcptr y) { mpfr_mul(r, x, y, MPFR_RNDN); } },
        { "x / T",   B::div_s, [](mpfr_ptr r, mpfr_srcptr x, mpfr_srcptr y) { mpfr_div(r, x, y, MPFR_RNDN); } },
        { "x.pow(y)", B::powi_iv, [](mpfr_ptr r, mpfr_srcptr x, mpfr_srcptr y) { mpfr_pow(r, x, y, MPFR_RNDN); } },
    };
    return ops;
}

template<class I>
I apply2(B id, I const& x, I const& y)
{
    switch (id) {
    case B::add:   return x + y;
    case B::sub:   return x - y;
    case B::mul:   return x * y;
    case B::div:   return x / y;
    case B::add_s: return x + y.inf();
    case B::sub_s: return x - y.inf();
    case B::mul_s: return x * y.inf();
    case B::div_s: return x / y.inf();
    case B::powi_iv: return x.pow(y);
    }
    return x;
}

bool is_scalar_variant(B id) { return id == B::add_s || id == B::sub_s || id == B::mul_s || id == B::div_s; }
bool is_div(B id) { return id == B::div || id == B::div_s; }

template<class T>
void check_binary(Tally& t, BOp const& op, T xlo, T xhi, T ylo, T yhi, T rlo, T rhi, rng_t& g, bool strict_domain = true)
{
    ++t.checks;
    bool const empty = !(rlo <= rhi);
    auto const where = [&] {
        return std::string(op.name) + " x=" + iv_str(xlo, xhi) + " y=" + iv_str(ylo, yhi) + " out=" + iv_str(rlo, rhi);
    };
    bool const y_is_zero = ylo == 0 && yhi == 0;
    bool const x_nonpos_only = op.id == B::powi_iv && !(xhi > 0);
    if ((is_div(op.id) && y_is_zero) || x_nonpos_only) {
        if (!empty && strict_domain) t.fail(where() + ": domain error must be empty");
        return;
    }
    if (empty) {
        t.fail(where() + ": empty result for operands meeting the domain");
        return;
    }
    std::vector<T> xs, ys;
    sample_points(xlo, xhi, g, xs);
    sample_points(ylo, yhi, g, ys);
    for (T x : xs) {
        for (T y : ys) {
            if (is_div(op.id) && y == 0) continue;
            if (op.id == B::powi_iv && !(x > 0)) continue;
            mpfr_set_d(S().x, x, MPFR_RNDN);
            mpfr_set_d(S().y, y, MPFR_RNDN);
            op.ref(S().r, S().x, S().y);
            if (!encloses(rlo, rhi, S().r))
                t.fail(where() + " at (" + hex(x) + ", " + hex(y) + ") exact=" + hex(mpfr_get_d(S().r, MPFR_RNDN)) + " not enclosed");
        }
    }
}

template<class T>
std::pair<T, T> bin_bounds(BOp const& op, int shape, bool second, rng_t& g)
{
    bool const f = std::is_same_v<T, float>;
    if (op.id == B::powi_iv) {
        if (!second) return make_bounds<T>(shape, g, [&] { return uni(g, 0.05, 12.0); });
        return make_bounds<T>(shape, g, [&] { return uni(g, -2.5, 3.5); });
    }
    return make_bounds<T>(shape, g, [&] { return g_binary(g, f); });
}

template<class T>
void run_binary(Tally& t, BOp const& op, int shape, rng_t& g, int n)
{
    using I = pappus::interval<T>;
    for (int i = 0; i < n; ++i) {
        auto [xlo, xhi] = bin_bounds<T>(op, shape, false, g);
        int const yshape = is_scalar_variant(op.id) ? POINT : static_cast<int>(g() % N_SHAPES);
        auto [ylo, yhi] = bin_bounds<T>(op, yshape, true, g);
        if (!std::isfinite(xlo) || !std::isfinite(xhi) || !std::isfinite(ylo) || !std::isfinite(yhi)) continue;
        I const res = apply2(op.id, I(xlo, xhi), I(ylo, yhi));
        check_binary<T>(t, op, xlo, xhi, ylo, yhi, res.inf(), res.sup(), g);
    }
}

template<class W>
void run_binary_wide(Tally& t, BOp const& op, int shape, rng_t& g, int n)
{
    using T = eve::element_type_t<W>;
    constexpr std::size_t L = W::size();
    using IW = pappus::interval<W>;
    for (int i = 0; i < n; i += static_cast<int>(L)) {
        std::array<T, L> xlo{}, xhi{}, ylo{}, yhi{};
        W wxlo{}, wxhi{}, wylo{}, wyhi{};
        for (std::size_t k = 0; k < L; ++k) {
            auto [a, b] = bin_bounds<T>(op, shape, false, g);
            int const yshape = is_scalar_variant(op.id) ? POINT : static_cast<int>(g() % N_SHAPES);
            auto [c, d] = bin_bounds<T>(op, yshape, true, g);
            if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c) || !std::isfinite(d)) { a = b = c = d = T(1); }
            xlo[k] = a; xhi[k] = b; ylo[k] = c; yhi[k] = d;
            wxlo.set(k, a); wxhi.set(k, b); wylo.set(k, c); wyhi.set(k, d);
        }
        IW const res = [&] {
            if (is_scalar_variant(op.id)) {
                // the scalar-operand overloads take a single value; for SIMD that is a lane vector
                switch (op.id) {
                case B::add_s: return IW(wxlo, wxhi) + wylo;
                case B::sub_s: return IW(wxlo, wxhi) - wylo;
                case B::mul_s: return IW(wxlo, wxhi) * wylo;
                default:       return IW(wxlo, wxhi) / wylo;
                }
            }
            return apply2(op.id, IW(wxlo, wxhi), IW(wylo, wyhi));
        }();
        for (std::size_t k = 0; k < L; ++k) {
            T const rlo = res.inf().get(k), rhi = res.sup().get(k);
            check_binary<T>(t, op, xlo[k], xhi[k], ylo[k], yhi[k], rlo, rhi, g, /*strict_domain=*/false);
        }
    }
}

// ULP distance between two floats of the same type in the ordered representation.
template<class T>
std::int64_t ordered(T v)
{
    using Int = std::conditional_t<std::is_same_v<T, float>, std::int32_t, std::int64_t>;
    Int const i = std::bit_cast<Int>(v);
    return i < 0 ? -static_cast<std::int64_t>(i & std::numeric_limits<Int>::max()) : static_cast<std::int64_t>(i);
}

template<class T>
std::int64_t width_ulps(T lo, T hi) { return ordered(hi) - ordered(lo); }

} // namespace

// ===========================================================================
// Soundness: scalar float / double
// ===========================================================================

TEMPLATE_TEST_CASE("soundness: unary operations enclose the exact result", "[soundness]", float, double)
{
    rng_t g = make_rng(0x5eed0001);
    for (auto const& op : unary_ops()) {
        for (int shape = 0; shape < N_SHAPES; ++shape) {
            Tally t;
            run_unary<TestType>(t, op, shape, g, per_shape());
            INFO(op.name << " / " << shape_name[shape] << ": " << t.first);
            CHECK(t.bad == 0);
        }
    }
}

TEMPLATE_TEST_CASE("soundness: binary operations enclose the exact result", "[soundness]", float, double)
{
    rng_t g = make_rng(0x5eed0002);
    for (auto const& op : binary_ops()) {
        for (int shape = 0; shape < N_SHAPES; ++shape) {
            Tally t;
            run_binary<TestType>(t, op, shape, g, per_shape());
            INFO(op.name << " / " << shape_name[shape] << ": " << t.first);
            CHECK(t.bad == 0);
        }
    }
}

// SIMD lanes go through the same ropd/ropu and must be just as sound, and
// agree with the scalar implementation on which lanes are empty.
TEMPLATE_TEST_CASE("soundness: SIMD lanes enclose the exact result", "[soundness][simd]",
                   eve::wide<float>, eve::wide<double>)
{
    rng_t g = make_rng(0x5eed0003);
    for (auto const& op : unary_ops()) {
        for (int shape = 0; shape < N_SHAPES; ++shape) {
            Tally t;
            run_unary_wide<TestType>(t, op, shape, g, per_shape());
            INFO("wide " << op.name << " / " << shape_name[shape] << ": " << t.first);
            CHECK(t.bad == 0);
        }
    }
    for (auto const& op : binary_ops()) {
        if (op.id == B::powi_iv) continue;
        for (int shape = 0; shape < N_SHAPES; ++shape) {
            Tally t;
            run_binary_wide<TestType>(t, op, shape, g, per_shape());
            INFO("wide " << op.name << " / " << shape_name[shape] << ": " << t.first);
            CHECK(t.bad == 0);
        }
    }
}

// ===========================================================================
// Overflow / underflow: bounds must stay bounds when the exact result is
// outside the representable range.
// ===========================================================================

TEMPLATE_TEST_CASE("soundness: range edges (overflow, underflow)", "[soundness]", float, double)
{
    using T = TestType;
    using I = pappus::interval<T>;
    constexpr T big = std::numeric_limits<T>::max();
    constexpr T tiny = std::numeric_limits<T>::denorm_min();
    constexpr T mn = std::numeric_limits<T>::min();
    constexpr T h = std::is_same_v<T, float> ? T(1e20) : T(1e160);
    constexpr T hs = std::is_same_v<T, float> ? T(1e-20) : T(1e-160);

    struct Case { char const* name; I got; std::function<void(mpfr_ptr)> ref; };
    auto bin = [&](char const* nm, I got, T a, T b, int (*f)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t)) {
        return Case{ nm, got, [a, b, f](mpfr_ptr r) { mpfr_set_d(S().x, a, MPFR_RNDN); mpfr_set_d(S().y, b, MPFR_RNDN); f(r, S().x, S().y, MPFR_RNDN); } };
    };
    auto un = [&](char const* nm, I got, T a, int (*f)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)) {
        return Case{ nm, got, [a, f](mpfr_ptr r) { mpfr_set_d(S().x, a, MPFR_RNDN); f(r, S().x, MPFR_RNDN); } };
    };
    T const ex = std::is_same_v<T, float> ? T(89.5) : T(710.5);

    std::vector<Case> cases;
    cases.push_back(bin("big + big", I(big) + I(big), big, big, mpfr_add));
    cases.push_back(bin("-big - big", I(-big) - I(big), -big, big, mpfr_sub));
    cases.push_back(bin("h * h", I(h) * I(h), h, h, mpfr_mul));
    cases.push_back(bin("-h * h", I(-h) * I(h), -h, h, mpfr_mul));
    cases.push_back(bin("hs * hs (underflow)", I(hs) * I(hs), hs, hs, mpfr_mul));
    cases.push_back(bin("-hs * hs (underflow)", I(-hs) * I(hs), -hs, hs, mpfr_mul));
    cases.push_back(bin("h / hs", I(h) / I(hs), h, hs, mpfr_div));
    cases.push_back(bin("hs / h (underflow)", I(hs) / I(h), hs, h, mpfr_div));
    cases.push_back(bin("tiny / 3 (underflow)", I(tiny) / I(T(3)), tiny, T(3), mpfr_div));
    cases.push_back(bin("mn / 3 (subnormal)", I(mn) / I(T(3)), mn, T(3), mpfr_div));
    cases.push_back(bin("big / 0.5", I(big) / I(T(0.5)), big, T(0.5), mpfr_div));
    cases.push_back(un("exp(overflow)", I(ex).exp(), ex, mpfr_exp));
    cases.push_back(un("exp(-overflow)", I(-ex).exp(), -ex, mpfr_exp));
    cases.push_back(un("exp(-745)", I(T(-745)).exp(), T(-745), mpfr_exp));
    cases.push_back(un("sinh(overflow)", I(ex).sinh(), ex, mpfr_sinh));
    cases.push_back(un("sinh(-overflow)", I(-ex).sinh(), -ex, mpfr_sinh));
    cases.push_back(un("cosh(overflow)", I(ex).cosh(), ex, mpfr_cosh));
    cases.push_back(un("cbrt(big)", I(big).cbrt(), big, mpfr_cbrt));
    cases.push_back(un("cbrt(-big)", I(-big).cbrt(), -big, mpfr_cbrt));
    cases.push_back(un("cbrt(denorm_min)", I(tiny).cbrt(), tiny, mpfr_cbrt));
    cases.push_back(un("sqrt(denorm_min)", I(tiny).sqrt(), tiny, mpfr_sqrt));
    cases.push_back(un("atan(big)", I(big).atan(), big, mpfr_atan));
    cases.push_back(un("log(denorm_min)", I(tiny).log(), tiny, mpfr_log));
    cases.push_back(un("log(big)", I(big).log(), big, mpfr_log));
    cases.push_back(un("tanh(big)", I(big).tanh(), big, mpfr_tanh));
    {
        I const p = I(h).pow(2);
        cases.push_back(Case{ "pow(h, 2)", p, [h](mpfr_ptr r) { mpfr_set_d(S().x, h, MPFR_RNDN); mpfr_sqr(r, S().x, MPFR_RNDN); } });
        I const q = I(h).pow(-3);
        cases.push_back(Case{ "pow(h, -3) (underflow)", q, [h](mpfr_ptr r) { mpfr_set_d(S().x, h, MPFR_RNDN); mpfr_pow_si(r, S().x, -3, MPFR_RNDN); } });
        I const w = I(hs).pow(-2);
        cases.push_back(Case{ "pow(hs, -2) (overflow)", w, [hs](mpfr_ptr r) { mpfr_set_d(S().x, hs, MPFR_RNDN); mpfr_pow_si(r, S().x, -2, MPFR_RNDN); } });
    }

    for (auto const& c : cases) {
        mpfr_t r;
        mpfr_init2(r, 256);
        c.ref(r);
        INFO(c.name << " -> " << iv_str(c.got.inf(), c.got.sup()) << ", exact " << hex(mpfr_get_d(r, MPFR_RNDN)));
        CHECK_FALSE(c.got.is_empty());
        CHECK(encloses(c.got.inf(), c.got.sup(), r));
        mpfr_clear(r);
    }
}

// ===========================================================================
// ULP width: a point interval must not be wider than the directed-rounding
// margin; exact operations are at most 1 ULP wide.
// ===========================================================================

TEMPLATE_TEST_CASE("soundness: point intervals have bounded ULP width", "[soundness][width]", float, double)
{
    using T = TestType;
    using I = pappus::interval<T>;
    rng_t g = make_rng(0x5eed0004);
    bool const f = std::is_same_v<T, float>;
    // Transcendentals are widened by k_eve_ulps steps in the evaluation
    // grid; a float narrowing adds at most one more step per endpoint.
    std::int64_t const trans_max = f ? 2 : 2 * pappus::fp::detail::k_eve_ulps + 2;

    for (auto const& op : unary_ops()) {
        std::int64_t limit;
        switch (op.id) {
        case U::abs: case U::floor: case U::ceil: limit = 0; break;
        case U::sqrt: case U::sqr: case U::inv: limit = 1; break;
        case U::powi: limit = op.p == 2 ? 1 : (op.p == -1 ? 1 : 2 * std::abs(op.p) + 2); break;   // chains of directed multiplications
        // exp(p * log(x)) combines two 16-step transcendental margins; the
        // multiplicative scale p and binade crossings magnify that width.
        case U::powr: limit = 800; break;
        default: limit = trans_max; break;
        }
        std::int64_t worst = 0;
        T worst_x = 0;
        for (int i = 0; i < 400; ++i) {
            T x = static_cast<T>(clamp_pow_base(op, op.gen(g, f), f));
            if (!std::isfinite(x) || !op.dom(x)) continue;
            // Out of scope for a *margin* regression test:
            //  - sqrt below the underflow threshold: bounds are deliberately
            //    moved one more ULP outward (fp::detail::tiny_threshold);
            //  - real powers of extreme bases: exp(p * log x) amplifies the
            //    log margin by |p log x|, so restrict the base to [1e-3, 1e3];
            if (op.id == U::sqrt && x < 1e-28) continue;
            if (op.id == U::powr) x = std::clamp<T>(x, T(1e-3), T(1e3));
            I const r = apply(op, I(x, x));
            if (r.is_empty() || !std::isfinite(r.inf()) || !std::isfinite(r.sup())) continue;
            // results in the subnormal range: ULP width is not meaningful there (and for exp
            // eve flushes to 0, so the upper bound is a deliberately coarse 2^-1021)
            if (std::max(std::fabs(r.inf()), std::fabs(r.sup())) < 4 * std::numeric_limits<T>::min()) continue;
            auto const w = width_ulps(r.inf(), r.sup());
            if (w > worst) { worst = w; worst_x = x; }
        }
        INFO(op.name << ": worst width " << worst << " ulp at x=" << hex(worst_x) << ", limit " << limit);
        CHECK(worst <= limit);
    }
    // binary: +, -, *, / of two points
    for (int i = 0; i < 2000; ++i) {
        T a = static_cast<T>(g_binary(g, f)), b = static_cast<T>(g_binary(g, f));
        for (auto const& [nm, r] : std::array<std::pair<char const*, I>, 4>{{
                 { "+", I(a) + I(b) }, { "-", I(a) - I(b) }, { "*", I(a) * I(b) }, { "/", b == 0 ? I(T(1)) : I(a) / I(b) } }}) {
            INFO(nm << " a=" << hex(a) << " b=" << hex(b));
            CHECK(width_ulps(r.inf(), r.sup()) <= 1);
        }
    }
}

// ===========================================================================
// Named regressions
// ===========================================================================

// Captured from operon: eve::cbrt is not monotone, so with a 1-ULP margin the
// bounds of cbrt([x, nextafter(x)]) inverted (lo > hi) and
// EXPECT(!(lo > hi)) aborted the process even in NDEBUG builds.
TEST_CASE("regression: cbrt of a 1-ULP float interval does not invert", "[soundness][regression]")
{
    using I = pappus::interval<float>;
    float const lo = 0x1.407b58p+2f;
    float const hi = std::nextafter(lo, std::numeric_limits<float>::infinity());
    I const r = I(lo, hi).cbrt();   // must not abort
    CHECK_FALSE(r.is_empty());
    CHECK(r.inf() <= r.sup());
    for (float x : { lo, hi }) {
        mpfr_set_d(S().x, x, MPFR_RNDN);
        mpfr_cbrt(S().r, S().x, MPFR_RNDN);
        INFO("cbrt(" << hex(x) << ") not in " << iv_str(r.inf(), r.sup()));
        CHECK(encloses(r.inf(), r.sup(), S().r));
    }
}

TEMPLATE_TEST_CASE("regression: cbrt of thin intervals is never inverted or empty", "[soundness][regression]", float, double)
{
    using T = TestType;
    using I = pappus::interval<T>;
    constexpr T inf = std::numeric_limits<T>::infinity();
    rng_t g = make_rng(0x5eed0005);
    std::size_t bad = 0, n = 0;
    for (int i = 0; i < 100000; ++i) {
        T x = static_cast<T>(smag(g, -10, 10));
        for (int k : { 1, 2, 3, 5 }) {
            T y = x;
            for (int j = 0; j < k; ++j) y = std::nextafter(y, inf);
            I const r = I(x, y).cbrt();
            ++n;
            if (r.is_empty() || !(r.inf() <= r.sup())) ++bad;
        }
    }
    CHECK(n > 0);
    CHECK(bad == 0);
}

// eve::div[lower] returns a value one ULP too high for ~45% of the inexact
// quotients (is_ltz(e) where is_gtz(e) is needed, upstream eve div.hpp).
// ropd<op_div> must be a true lower bound regardless.
TEST_CASE("regression: float division lower bound (eve::div[lower] workaround)", "[soundness][regression]")
{
    namespace fp = pappus::fp;
    float const a = 0x1.2ce1c2p+0f, b = -0x1.856c22p+2f;
    float const lo = fp::ropd<fp::op_div>(a, b), hi = fp::ropu<fp::op_div>(a, b);
    mpfr_set_d(S().x, a, MPFR_RNDN);
    mpfr_set_d(S().y, b, MPFR_RNDN);
    mpfr_div(S().r, S().x, S().y, MPFR_RNDN);
    INFO("a=" << hex(a) << " b=" << hex(b) << " RD=" << hex(lo) << " RU=" << hex(hi));
    CHECK(encloses(lo, hi, S().r));
    CHECK(lo < hi);   // inexact quotient: exactly one ULP wide
    CHECK(width_ulps(lo, hi) == 1);

    using I = pappus::interval<float>;
    I const q = I(a, a) / I(b, b);
    CHECK(encloses(q.inf(), q.sup(), S().r));
    I const inv = I(b, b).inv();
    mpfr_ui_div(S().r, 1, S().y, MPFR_RNDN);
    CHECK(encloses(inv.inf(), inv.sup(), S().r));
}

TEMPLATE_TEST_CASE("regression: Eve lower division counterexamples are repaired", "[soundness][regression]", float, double)
{
    namespace fp = pappus::fp;
    using T = TestType;
    T const a = std::is_same_v<T, float> ? T(0x1.2fac24p+2f) : T(0x1.ce090a1a1b2bp+0);
    T const b = std::is_same_v<T, float> ? T(-0x1.8bfebcp+0f) : T(-0x1.894babf1d82e8p+1);
    mpfr_set_d(S().x, a, MPFR_RNDN);
    mpfr_set_d(S().y, b, MPFR_RNDN);
    mpfr_div(S().r, S().x, S().y, MPFR_RNDN);
    T const eve_lower = eve::div[eve::lower](a, b);
    mpfr_set_d(S().a, eve_lower, MPFR_RNDN);
    CHECK(mpfr_greater_p(S().a, S().r));
    CHECK(encloses(fp::ropd<fp::op_div>(a, b), fp::ropu<fp::op_div>(a, b), S().r));
}

TEMPLATE_TEST_CASE("regression: nonzero subnormal unary results are enclosed", "[soundness][regression]", float, double)
{
    using T = TestType;
    using I = pappus::interval<T>;
    T const x = std::numeric_limits<T>::denorm_min();
    struct Case { char const* name; I result; int (*ref)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t); };
    std::array<Case, 7> const cases {{
        { "sin", I(x).sin(), mpfr_sin }, { "tan", I(x).tan(), mpfr_tan },
        { "asin", I(x).asin(), mpfr_asin }, { "atan", I(x).atan(), mpfr_atan },
        { "sinh", I(x).sinh(), mpfr_sinh }, { "tanh", I(x).tanh(), mpfr_tanh },
        { "log1p", I(x).log1p(), mpfr_log1p },
    }};
    mpfr_set_d(S().x, x, MPFR_RNDN);
    for (auto const& c : cases) {
        c.ref(S().r, S().x, MPFR_RNDN);
        INFO(c.name << "(" << hex(x) << ") -> " << iv_str(c.result.inf(), c.result.sup()));
        CHECK(encloses(c.result.inf(), c.result.sup(), S().r));
    }
}

TEST_CASE("regression: exp of an empty interval remains canonical empty", "[soundness][regression]")
{
    CHECK(pappus::interval<float>::empty().exp() == pappus::interval<float>::empty());
    CHECK(pappus::interval<double>::empty().exp() == pappus::interval<double>::empty());
}

TEMPLATE_TEST_CASE("regression: SIMD exp preserves empty lanes", "[soundness][regression][simd]",
                   eve::wide<float>, eve::wide<double>)
{
    using W = TestType;
    using T = eve::element_type_t<W>;
    W lo { [](auto i, auto) { return (i % 2) == 0 ? T(2) : T(-1); } };
    W hi { [](auto i, auto) { return (i % 2) == 0 ? T(1) : T(1); } };
    pappus::interval<W> const r(lo, hi);
    auto const e = r.exp();
    for (std::size_t i = 0; i < W::size(); ++i) {
        CHECK(static_cast<bool>(e.is_empty().get(i)) == ((i % 2) == 0));
    }
}

TEMPLATE_TEST_CASE("regression: ropd/ropu<op_div> are the exact directed quotients", "[soundness][regression]", float, double)
{
    namespace fp = pappus::fp;
    using T = TestType;
    rng_t g = make_rng(0x5eed0006);
    std::size_t bad = 0, inexact = 0;
    std::string first;
    for (int i = 0; i < 200000; ++i) {
        T a = static_cast<T>(uni(g, -10, 10)), b = static_cast<T>(uni(g, -10, 10));
        if (b == 0) continue;
        T const lo = fp::ropd<fp::op_div>(a, b), hi = fp::ropu<fp::op_div>(a, b);
        mpfr_set_d(S().x, a, MPFR_RNDN);
        mpfr_set_d(S().y, b, MPFR_RNDN);
        mpfr_div(S().r, S().x, S().y, MPFR_RNDN);
        bool const exact = mpfr_cmp_d(S().r, static_cast<double>(lo)) == 0 && lo == hi;
        // RD(a/b) is the largest T <= a/b, RU(a/b) the smallest T >= a/b:
        // enclose, and be at most 1 step apart (== 1 step iff the quotient is inexact).
        bool ok = encloses(lo, hi, S().r) && width_ulps(lo, hi) <= 1;
        if (!exact) { ++inexact; ok = ok && width_ulps(lo, hi) == 1; }
        if (!ok && bad++ == 0) first = "a=" + hex(a) + " b=" + hex(b) + " RD=" + hex(lo) + " RU=" + hex(hi);
    }
    INFO(first);
    CHECK(bad == 0);
    CHECK(inexact > 100000);
}
// If this reports eve::div[lower] as fixed, fp::detail::div_rd can use
// eve::div[eve::lower] directly. Never fails: it only documents the upstream bug.
TEMPLATE_TEST_CASE("eve::div[lower] upstream bug sentinel", "[soundness][eve]", float, double)
{
    using T = TestType;
    rng_t g = make_rng(0x5eed0007);
    std::size_t miss = 0, n = 0;
    for (int i = 0; i < 100000; ++i) {
        T a = static_cast<T>(uni(g, -10, 10)), b = static_cast<T>(uni(g, -10, 10));
        if (b == 0) continue;
        T const lo = eve::div[eve::lower](a, b);
        mpfr_set_d(S().x, a, MPFR_RNDN);
        mpfr_set_d(S().y, b, MPFR_RNDN);
        mpfr_div(S().r, S().x, S().y, MPFR_RNDN);
        mpfr_set_d(S().a, static_cast<double>(lo), MPFR_RNDN);
        ++n;
        if (mpfr_greater_p(S().a, S().r)) ++miss;
    }
    if (miss == 0)
        WARN("eve::div[lower] returned a valid lower bound for all " << n << " samples: the upstream bug is fixed, simplify fp::detail::div_rd");
    else
        INFO("eve::div[lower] is not a lower bound for " << miss << " / " << n << " samples (workaround in fp::detail::div_rd is needed)");
    SUCCEED();
}

// ===========================================================================
// Non-aborting bound checks
// ===========================================================================

TEST_CASE("check_bounds: inverted bounds give the empty interval, not an abort", "[soundness][contracts]")
{
#ifdef NDEBUG
    WARN("check_bounds contract test: NDEBUG build (the old EXPECT aborted here too: libassert ignores NDEBUG)");
#else
    WARN("check_bounds contract test: assertions-enabled build");
#endif
    using IF = pappus::interval<float>;
    using ID = pappus::interval<double>;

    CHECK(IF(2.0f, 1.0f).is_empty());
    CHECK(ID(2.0, 1.0).is_empty());
    CHECK(IF(2.0f, 1.0f) == IF::empty());
    CHECK(ID(std::nextafter(1.0, 2.0), 1.0) == ID::empty());   // 1-ULP inversion

    IF a;
    a.set_bounds(3.0f, -3.0f);
    CHECK(a.is_empty());
    a.set_bounds(std::pair<float, float>{ 5.0f, 4.0f });
    CHECK(a.is_empty());

    // valid and degenerate pairs are unchanged (signed-zero convention kept)
    IF ok(1.0f, 2.0f);
    CHECK(!ok.is_empty());
    CHECK(ok.inf() == 1.0f);
    CHECK(ok.sup() == 2.0f);
    IF z(-0.0f, 0.0f);
    CHECK(!z.is_empty());
    CHECK(std::signbit(z.sup()) == true);
    CHECK(std::signbit(z.inf()) == false);

    // empty stays empty, and consumers of is_empty() see it
    CHECK((IF(2.0f, 1.0f) + IF(1.0f, 2.0f)).is_empty());
    CHECK(IF(2.0f, 1.0f).exp().is_empty());
}

TEMPLATE_TEST_CASE("check_bounds: inverted SIMD lanes become empty lanes", "[soundness][contracts][simd]",
                   eve::wide<float>, eve::wide<double>)
{
    using W = TestType;
    using T = eve::element_type_t<W>;
    constexpr std::size_t L = W::size();
    W lo{}, hi{};
    for (std::size_t k = 0; k < L; ++k) {
        bool const inverted = (k % 2) == 1;
        lo.set(k, static_cast<T>(k) + (inverted ? T(2) : T(0)));
        hi.set(k, static_cast<T>(k) + T(1));
        if (inverted) { lo.set(k, T(10)); hi.set(k, T(5)); }
    }
    pappus::interval<W> const iv(lo, hi);   // must not abort
    for (std::size_t k = 0; k < L; ++k) {
        bool const inverted = (k % 2) == 1;
        INFO("lane " << k);
        CHECK(static_cast<bool>(iv.is_empty().get(k)) == inverted);
        if (!inverted) {
            CHECK(iv.inf().get(k) == static_cast<T>(k));
            CHECK(iv.sup().get(k) == static_cast<T>(k) + T(1));
        }
    }
}

// ===========================================================================
// Margin guards: the ULP margins in fp/math.hpp must stay above the measured
// kernel error. These fail if an eve (or libm) update makes a function worse.
// ===========================================================================

namespace {

double ulp_error(double got, mpfr_srcptr r)
{
    mpfr_set_d(S().a, got, MPFR_RNDN);
    mpfr_sub(S().b, S().a, r, MPFR_RNDN);
    mpfr_abs(S().b, S().b, MPFR_RNDN);
    // ulp of the exact value in double precision
    if (mpfr_zero_p(r)) return 0.0;
    mpfr_exp_t const e = mpfr_get_exp(r);   // r = m * 2^e, 0.5 <= |m| < 1
    mpfr_div_2si(S().b, S().b, e - 53, MPFR_RNDN);
    return mpfr_get_d(S().b, MPFR_RNDN);
}

struct KernelCase {
    char const* name;
    double (*eve_fn)(double);
    double (*libm_fn)(double);
    int (*ref)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
    double lo, hi;
    bool logmode;
};

} // namespace

TEST_CASE("eve double kernels and libm stay within 3/4 of the ULP margins in fp/math.hpp", "[soundness][margin]")
{
    std::vector<KernelCase> const cases = {
        { "exp",   [](double x) { return eve::exp(x); },   [](double x) { return std::exp(x); },   mpfr_exp,   -5, 5, false },
        { "log",   [](double x) { return eve::log(x); },   [](double x) { return std::log(x); },   mpfr_log,   -3, 3, true },
        { "log1p", [](double x) { return eve::log1p(x); }, [](double x) { return std::log1p(x); }, mpfr_log1p, -0.9, 10, false },
        { "sin",   [](double x) { return eve::sin(x); },   [](double x) { return std::sin(x); },   mpfr_sin,   -3, 3, false },
        { "cos",   [](double x) { return eve::cos(x); },   [](double x) { return std::cos(x); },   mpfr_cos,   -3, 3, false },
        { "tan",   [](double x) { return eve::tan(x); },   [](double x) { return std::tan(x); },   mpfr_tan,   -1.5, 1.5, false },
        { "asin",  [](double x) { return eve::asin(x); },  [](double x) { return std::asin(x); },  mpfr_asin,  -1, 1, false },
        { "acos",  [](double x) { return eve::acos(x); },  [](double x) { return std::acos(x); },  mpfr_acos,  -1, 1, false },
        { "atan",  [](double x) { return eve::atan(x); },  [](double x) { return std::atan(x); },  mpfr_atan,  -10, 10, false },
        { "sinh",  [](double x) { return eve::sinh(x); },  [](double x) { return std::sinh(x); },  mpfr_sinh,  -5, 5, false },
        { "cosh",  [](double x) { return eve::cosh(x); },  [](double x) { return std::cosh(x); },  mpfr_cosh,  -5, 5, false },
        { "tanh",  [](double x) { return eve::tanh(x); },  [](double x) { return std::tanh(x); },  mpfr_tanh,  -5, 5, false },
        { "cbrt",  [](double x) { return eve::cbrt(x); },  [](double x) { return std::cbrt(x); },  mpfr_cbrt,  -10, 10, false },
    };
    rng_t g = make_rng(0x5eed0008);
    // The margin is expressed in the kernel-result grid, which can be half
    // the exact-result ULP at a binade boundary.
    double const eve_limit = 0.5 * pappus::fp::detail::k_eve_ulps;
    double const libm_limit = 0.75 * pappus::fp::detail::k_libm_ulps;
    for (auto const& c : cases) {
        double worst_eve = 0, worst_libm = 0;
        for (int i = 0; i < 20000; ++i) {
            double const x = c.logmode ? mag(g, c.lo, c.hi) : uni(g, c.lo, c.hi);
            mpfr_set_d(S().x, x, MPFR_RNDN);
            c.ref(S().r, S().x, MPFR_RNDN);
            worst_eve = std::max(worst_eve, ulp_error(c.eve_fn(x), S().r));
            worst_libm = std::max(worst_libm, ulp_error(c.libm_fn(x), S().r));
        }
        INFO(c.name << ": eve worst " << worst_eve << " ulp (limit " << eve_limit << "), libm worst " << worst_libm
                    << " ulp (limit " << libm_limit << ")");
        CHECK(worst_eve <= eve_limit);
        // widen_lo/widen_hi (affine_form) use k_libm_ulps
        CHECK(worst_libm <= libm_limit);
    }
}

// fp::widen_lo / widen_hi (used by affine_form on libm values) must enclose
// the exact value of std:: results, in float and double.
TEMPLATE_TEST_CASE("fp::widen_lo/widen_hi enclose libm results", "[soundness][margin]", float, double)
{
    using T = TestType;
    namespace fp = pappus::fp;
    struct LCase {
        char const* name;
        T (*fn)(T);
        int (*ref)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
        double lo, hi;
    };
    std::vector<LCase> const cases = {
        { "exp", [](T x) { return std::exp(x); }, mpfr_exp, -20, 20 },
        { "log", [](T x) { return std::log(x); }, mpfr_log, 1e-3, 1e3 },
        { "log1p", [](T x) { return std::log1p(x); }, mpfr_log1p, -0.9, 10 },
        { "sin", [](T x) { return std::sin(x); }, mpfr_sin, -10, 10 },
        { "cos", [](T x) { return std::cos(x); }, mpfr_cos, -10, 10 },
        { "tan", [](T x) { return std::tan(x); }, mpfr_tan, -1.5, 1.5 },
        { "asin", [](T x) { return std::asin(x); }, mpfr_asin, -1, 1 },
        { "acos", [](T x) { return std::acos(x); }, mpfr_acos, -1, 1 },
        { "atan", [](T x) { return std::atan(x); }, mpfr_atan, -10, 10 },
        { "sinh", [](T x) { return std::sinh(x); }, mpfr_sinh, -8, 8 },
        { "cosh", [](T x) { return std::cosh(x); }, mpfr_cosh, -8, 8 },
        { "tanh", [](T x) { return std::tanh(x); }, mpfr_tanh, -4, 4 },
        { "cbrt", [](T x) { return std::cbrt(x); }, mpfr_cbrt, -10, 10 },
    };
    rng_t g = make_rng(0x5eed0009);
    for (auto const& c : cases) {
        std::size_t bad = 0;
        double bad_x = 0;
        for (int i = 0; i < 20000; ++i) {
            T const x = static_cast<T>(uni(g, c.lo, c.hi));
            mpfr_set_d(S().x, x, MPFR_RNDN);
            c.ref(S().r, S().x, MPFR_RNDN);
            T const y = c.fn(x);
            if (!encloses(fp::widen_lo(y), fp::widen_hi(y), S().r) && bad++ == 0) bad_x = x;
        }
        INFO(c.name << " first miss at x=" << hex(bad_x));
        CHECK(bad == 0);
    }
}
