//
// Soundness fuzz tests for interval and affine arithmetic.
//
// For each operation we:
//   1. generate random input intervals (varied magnitudes, signs, zero-crossing)
//   2. compute the interval/affine result
//   3. sample scalar points from the input interval(s)
//   4. compute the scalar reference result for each sample
//   5. verify every scalar result is contained in the output interval
//
// A violation means the library produced an unsound (too-narrow) enclosure.
//
// Transcendental operations use eve-based approximations that may differ from
// libm by 1-2 ULP, so a small tolerance (in ULPs) is used for those.
// Basic arithmetic uses strict (0-ULP) checking.
//

#include <catch2/catch_all.hpp>

#include <cmath>
#include <cstdio>
#include <functional>
#include <random>
#include <string>
#include <vector>

#include "pappus/pappus.hpp"

using T = double;
using I = pappus::interval<T>;
using AF = pappus::affine_form<T>;

namespace {

constexpr int N_ITER = 200;    // random inputs per section
constexpr int N_SAMPLE = 100;  // scalar evaluations per unary check
constexpr int N_GRID = 14;     // samples per dimension for binary checks

// Deterministic seed for reproducibility.
std::mt19937_64 rng(987654321ULL);

T urand(T lo, T hi)
{
    std::uniform_real_distribution<T> d(lo, hi);
    return d(rng);
}

int irand(int lo, int hi)
{
    std::uniform_int_distribution<int> d(lo, hi);
    return d(rng);
}

// ---------------------------------------------------------------------------
// Interval generators — produce diverse finite intervals
// ---------------------------------------------------------------------------

I random_interval()
{
    switch (irand(0, 11)) {
    case 0: case 1: { // normal range [-10,10]
        T a = urand(-10, 10), b = urand(-10, 10);
        return a <= b ? I(a, b) : I(b, a);
    }
    case 2: { // straddling zero
        return I(urand(-10, -0.01), urand(0.01, 10));
    }
    case 3: { // positive
        return I(urand(0.01, 2), urand(2, 10));
    }
    case 4: { // negative
        return I(urand(-10, -2), urand(-2, -0.01));
    }
    case 5: { // large magnitude
        T a = urand(1e2, 1e6), b = urand(1e2, 1e6);
        if (a > b) std::swap(a, b);
        return irand(0, 1) ? I(a, b) : I(-b, -a);
    }
    case 6: { // very large magnitude
        T a = urand(1e8, 1e15), b = urand(1e8, 1e15);
        if (a > b) std::swap(a, b);
        return irand(0, 1) ? I(a, b) : I(-b, -a);
    }
    case 7: { // tiny around a point
        T mid = urand(-1, 1), r = urand(1e-12, 1e-4);
        return I(mid - r, mid + r);
    }
    case 8: { // point interval
        T v = urand(-5, 5);
        return I(v, v);
    }
    case 9: { // symmetric
        T m = urand(0.01, 8);
        return I(-m, m);
    }
    case 10: { // [0, b]
        return I(T(0), urand(0.1, 10));
    }
    case 11: { // [a, 0]
        return I(urand(-10, -0.1), T(0));
    }
    default:
        return I(-1, 1);
    }
}

// Bounded interval: |endpoints| <= limit, for overflow-prone ops.
I bounded_interval(T limit = 50)
{
    switch (irand(0, 5)) {
    case 0: case 1: {
        T a = urand(-limit, limit), b = urand(-limit, limit);
        return a <= b ? I(a, b) : I(b, a);
    }
    case 2: { return I(urand(-limit, -0.01), urand(0.01, limit)); }
    case 3: { return I(urand(0.01, limit * 0.3), urand(limit * 0.3, limit)); }
    case 4: { T v = urand(-limit, limit); return I(v, v); }
    case 5: { T m = urand(0.01, limit); return I(-m, m); }
    default: return I(-1, 1);
    }
}

// Domain-constrained generators.

I bounded_iv_50()     { return bounded_interval(50); }
I pos_iv() { // inf > 0 — for log (unbounded version)
    I iv = random_interval();
    if (iv.inf() <= 0)
        return I(urand(1e-6, 0.5), urand(0.5, 10));
    return iv;
}

I above_neg1_iv() { // inf > -1 — for log1p
    I iv = random_interval();
    if (iv.inf() <= -1)
        return I(urand(-0.99, 0), urand(0, 5));
    return iv;
}

I unit_iv() { // subset of [-1,1] — for asin/acos
    return I(urand(-1, 0), urand(0, 1));
}

I nneg_iv() { // inf >= 0 — for sqrt
    I iv = random_interval();
    if (iv.inf() < 0)
        return I(T(0), std::max(iv.sup(), T(1)));
    return iv;
}

I nz_iv() { // does not contain zero — for inv
    I iv = random_interval();
    if (iv.inf() <= 0 && iv.sup() >= 0)
        return irand(0, 1) ? I(urand(0.01, 5), urand(5, 10))
                           : I(urand(-10, -5), urand(-5, -0.01));
    return iv;
}

I tan_iv() { // within one branch of tan, away from asymptotes
    T c = urand(-1.2, 1.2), r = urand(0.01, 0.3);
    return I(c - r, c + r);
}

I nz_iv_bounded() {
    I iv = bounded_iv_50();
    if (iv.inf() <= 0 && iv.sup() >= 0)
        return irand(0, 1) ? I(urand(0.01, 5), urand(5, 10))
                           : I(urand(-10, -5), urand(-5, -0.01));
    return iv;
}

std::vector<T> sample_points(I const& iv, int n)
{
    std::vector<T> pts;
    pts.reserve(n + 3);
    pts.push_back(iv.inf());
    pts.push_back(iv.sup());
    if (iv.inf() != iv.sup())
        pts.push_back(iv.mid());
    while ((int)pts.size() < n)
        pts.push_back(urand(iv.inf(), iv.sup()));
    return pts;
}

// Widen a value by tol_ulps ULPs toward dir (+1 = +inf, -1 = -inf).
T widen(T v, int tol_ulps, int dir)
{
    if (!std::isfinite(v)) return v;
    for (int i = 0; i < tol_ulps; ++i)
        v = std::nextafter(v, dir > 0 ? std::numeric_limits<T>::infinity()
                                      : -std::numeric_limits<T>::infinity());
    return v;
}

// Widen interval bounds proportionally to |center| to account for affine
// construction rounding, plus tol_ulps for eve/Chebyshev approximation error.
struct sound_bounds {
    T lo, hi;
    sound_bounds(I const& iv, int tol_ulps)
    {
        lo = iv.inf();
        hi = iv.sup();
        for (int i = 0; i < tol_ulps; ++i) {
            lo = std::nextafter(lo, -std::numeric_limits<T>::infinity());
            hi = std::nextafter(hi, std::numeric_limits<T>::infinity());
        }
    }
    bool contains(T v) const { return !std::isnan(v) && lo <= v && v <= hi; }
};

// Unary: verify f(x) ∈ result (with tol_ulps widening) for all sampled x.
template<typename F>
void check_unary_sound(I const& input, I const& result, F f,
                       const char* op_name, int tol_ulps = 0)
{
    sound_bounds b(result, tol_ulps);
    for (T x : sample_points(input, N_SAMPLE)) {
        T r = f(x);
        if (std::isnan(r)) continue;
        if (!b.contains(r)) {
            FAIL_CHECK(op_name << ": f(" << x << ")=" << r
                       << " not in [" << b.lo << ", " << b.hi << "]"
                       << "  input=[" << input.inf() << ", " << input.sup() << "]");
            return;
        }
    }
}

// Binary: verify f(x,y) ∈ result for all sampled (x,y).
template<typename F>
void check_binary_sound(I const& xiv, I const& yiv, I const& result, F f,
                        const char* op_name, int tol_ulps = 0)
{
    sound_bounds b(result, tol_ulps);
    auto xs = sample_points(xiv, N_GRID);
    auto ys = sample_points(yiv, N_GRID);
    for (T x : xs) {
        for (T y : ys) {
            T r = f(x, y);
            if (std::isnan(r)) continue;
            if (!b.contains(r)) {
                FAIL_CHECK(op_name << ": f(" << x << ", " << y << ")=" << r
                           << " not in [" << b.lo << ", " << b.hi << "]"
                           << "  x=[" << xiv.inf() << ", " << xiv.sup() << "]"
                           << "  y=[" << yiv.inf() << ", " << yiv.sup() << "]");
                return;
            }
        }
    }
}

// Affine-aware version: widens proportionally to center magnitude.
template<typename F>
void check_aa_unary_sound(I const& input, AF const& result_af, F f,
                          const char* op_name, int extra_ulps = 0)
{
    I result = result_af.to_interval();
    // Proportional widening: center was rounded by up to 0.5 ULP, which
    // translates to many ULPs of the bounds when |center| >> |bound|.
    auto mag = std::fabs(result_af.center());
    auto w = std::ldexp(mag, 2 - std::numeric_limits<T>::digits);
    T lo = result.inf(), hi = result.sup();
    if (std::isfinite(w) && w > T(0)) {
        lo -= w;
        hi += w;
    }
    for (int i = 0; i < extra_ulps; ++i) {
        lo = std::nextafter(lo, -std::numeric_limits<T>::infinity());
        hi = std::nextafter(hi, std::numeric_limits<T>::infinity());
    }
    for (T x : sample_points(input, N_SAMPLE)) {
        T r = f(x);
        if (std::isnan(r)) continue;
        if (!(lo <= r && r <= hi)) {
            FAIL_CHECK(op_name << ": f(" << x << ")=" << r
                       << " not in [" << lo << ", " << hi << "]"
                       << "  input=[" << input.inf() << ", " << input.sup() << "]");
            return;
        }
    }
}

template<typename F>
void check_aa_binary_sound(I const& xiv, I const& yiv, AF const& result_af, F f,
                           const char* op_name, int extra_ulps = 0)
{
    I result = result_af.to_interval();
    auto mag = std::fabs(result_af.center());
    auto w = std::ldexp(mag, 2 - std::numeric_limits<T>::digits);
    T lo = result.inf(), hi = result.sup();
    if (std::isfinite(w) && w > T(0)) {
        lo -= w;
        hi += w;
    }
    for (int i = 0; i < extra_ulps; ++i) {
        lo = std::nextafter(lo, -std::numeric_limits<T>::infinity());
        hi = std::nextafter(hi, std::numeric_limits<T>::infinity());
    }
    auto xs = sample_points(xiv, N_GRID);
    auto ys = sample_points(yiv, N_GRID);
    for (T x : xs) {
        for (T y : ys) {
            T r = f(x, y);
            if (std::isnan(r)) continue;
            if (!(lo <= r && r <= hi)) {
                FAIL_CHECK(op_name << ": f(" << x << ", " << y << ")=" << r
                           << " not in [" << lo << ", " << hi << "]"
                           << "  x=[" << xiv.inf() << ", " << xiv.sup() << "]"
                           << "  y=[" << yiv.inf() << ", " << yiv.sup() << "]");
                return;
            }
        }
    }
}

// Tolerances: 0 for exact arithmetic, 2 for eve-based interval transcendentals,
// 8 for affine ops (accumulated center/radius + Chebyshev rounding).
constexpr int TOL_EXACT = 0;
constexpr int TOL_IV_TRAN = 2;
constexpr int TOL_AA = 8;

} // namespace

// ===========================================================================
// Interval arithmetic soundness
// ===========================================================================

TEST_CASE("fuzz: interval binary ops", "[fuzz]")
{
    SECTION("add") {
        for (int i = 0; i < N_ITER; ++i) {
            I x = random_interval(), y = random_interval();
            check_binary_sound(x, y, x + y, [](T a, T b) { return a + b; }, "add", TOL_EXACT);
        }
    }
    SECTION("sub") {
        for (int i = 0; i < N_ITER; ++i) {
            I x = random_interval(), y = random_interval();
            check_binary_sound(x, y, x - y, [](T a, T b) { return a - b; }, "sub", TOL_EXACT);
        }
    }
    SECTION("mul") {
        for (int i = 0; i < N_ITER; ++i) {
            I x = random_interval(), y = random_interval();
            check_binary_sound(x, y, x * y, [](T a, T b) { return a * b; }, "mul", TOL_EXACT);
        }
    }
    SECTION("div") {
        for (int i = 0; i < N_ITER; ++i) {
            I x = random_interval(), y = nz_iv();
            check_binary_sound(x, y, x / y, [](T a, T b) { return a / b; }, "div", TOL_EXACT);
        }
    }
}

TEST_CASE("fuzz: interval unary ops", "[fuzz]")
{
#define FUZZ_IV(sec, gen, method, sfn, tol)                                   \
    SECTION(sec) {                                                            \
        for (int i = 0; i < N_ITER; ++i) {                                    \
            I x = gen();                                                      \
            check_unary_sound(x, x.method(), [](T v) { return sfn; }, sec, tol); \
        }                                                                     \
    }

    FUZZ_IV("neg",     random_interval, operator-,  -v,             TOL_EXACT)
    FUZZ_IV("square",  random_interval, square,     v* v,           TOL_EXACT)
    FUZZ_IV("abs",     random_interval, abs,        std::fabs(v),   TOL_EXACT)
    FUZZ_IV("floor",   random_interval, floor,      std::floor(v),  TOL_EXACT)
    FUZZ_IV("ceil",    random_interval, ceil,       std::ceil(v),   TOL_EXACT)
    FUZZ_IV("exp",     bounded_iv_50,   exp,        std::exp(v),     TOL_IV_TRAN)
    FUZZ_IV("log",     pos_iv,          log,        std::log(v),     TOL_IV_TRAN)
    FUZZ_IV("log1p",   above_neg1_iv,   log1p,      std::log1p(v),   TOL_IV_TRAN)
    FUZZ_IV("sqrt",    nneg_iv,         sqrt,       std::sqrt(v),    TOL_IV_TRAN)
    FUZZ_IV("inv",     nz_iv,           inv,        T(1) / v,        TOL_EXACT)
    FUZZ_IV("cbrt",    random_interval, cbrt,       std::cbrt(v),    TOL_IV_TRAN)
    FUZZ_IV("sin",     random_interval, sin,        std::sin(v),     TOL_IV_TRAN)
    FUZZ_IV("cos",     random_interval, cos,        std::cos(v),     TOL_IV_TRAN)
    FUZZ_IV("tan",     tan_iv,          tan,        std::tan(v),     TOL_IV_TRAN)
    FUZZ_IV("asin",    unit_iv,         asin,       std::asin(v),   TOL_IV_TRAN)
    FUZZ_IV("acos",    unit_iv,         acos,       std::acos(v),    TOL_IV_TRAN)
    FUZZ_IV("atan",    random_interval, atan,       std::atan(v),    TOL_IV_TRAN)
    FUZZ_IV("sinh",    bounded_iv_50,   sinh,       std::sinh(v),    TOL_IV_TRAN)
    FUZZ_IV("cosh",    bounded_iv_50,   cosh,       std::cosh(v),    TOL_IV_TRAN)
    FUZZ_IV("tanh",    random_interval, tanh,       std::tanh(v),    TOL_IV_TRAN)

#undef FUZZ_IV
}

TEST_CASE("fuzz: interval power", "[fuzz]")
{
    SECTION("pow(int)") {
        for (int i = 0; i < N_ITER; ++i) {
            int p = irand(-5, 10);
            I x = (p < 0) ? nz_iv() : random_interval();
            check_unary_sound(x, x.pow(p), [p](T v) { return std::pow(v, p); }, "pow(int)", TOL_EXACT);
        }
    }
    SECTION("pow(T) fractional, positive base") {
        for (int i = 0; i < N_ITER; ++i) {
            T p = urand(0.1, 5.0);
            I x = pos_iv();
            check_unary_sound(x, x.pow(p), [p](T v) { return std::pow(v, p); }, "pow(T)", TOL_IV_TRAN);
        }
    }
    SECTION("pow(T) negative exponent, positive base") {
        for (int i = 0; i < N_ITER; ++i) {
            T p = urand(-4.0, -0.1);
            I x = pos_iv();
            check_unary_sound(x, x.pow(p), [p](T v) { return std::pow(v, p); }, "pow(T) neg", TOL_IV_TRAN);
        }
    }
}

TEST_CASE("fuzz: interval composite ops", "[fuzz]")
{
    SECTION("aq") {
        for (int i = 0; i < N_ITER; ++i) {
            I x = random_interval(), y = random_interval();
            check_binary_sound(x, y, pappus::ops::aq(x, y),
                [](T a, T b) { return a / std::sqrt(1 + b * b); }, "aq", TOL_IV_TRAN);
        }
    }
    SECTION("sqrtabs") {
        for (int i = 0; i < N_ITER; ++i) {
            I x = random_interval();
            check_unary_sound(x, pappus::ops::sqrtabs(x),
                [](T v) { return std::sqrt(std::fabs(v)); }, "sqrtabs", TOL_IV_TRAN);
        }
    }
    SECTION("logabs") {
        for (int i = 0; i < N_ITER; ++i) {
            I x = nz_iv();
            check_unary_sound(x, pappus::ops::logabs(x),
                [](T v) { return std::log(std::fabs(v)); }, "logabs", TOL_IV_TRAN);
        }
    }
}

// ===========================================================================
// Affine arithmetic soundness
// ===========================================================================

TEST_CASE("fuzz: affine binary ops (independent)", "[fuzz]")
{
    pappus::affine_context ctx;

    SECTION("add") {
        for (int i = 0; i < N_ITER; ++i) {
            I xi = random_interval(), yi = random_interval();
            AF x(ctx, xi), y(ctx, yi);
            check_aa_binary_sound(xi, yi, x + y,
                [](T a, T b) { return a + b; }, "aa:add", TOL_AA);
        }
    }
    SECTION("sub") {
        for (int i = 0; i < N_ITER; ++i) {
            I xi = random_interval(), yi = random_interval();
            AF x(ctx, xi), y(ctx, yi);
            check_aa_binary_sound(xi, yi, x - y,
                [](T a, T b) { return a - b; }, "aa:sub", TOL_AA);
        }
    }
    SECTION("mul") {
        for (int i = 0; i < N_ITER; ++i) {
            I xi = random_interval(), yi = random_interval();
            AF x(ctx, xi), y(ctx, yi);
            check_aa_binary_sound(xi, yi, x * y,
                [](T a, T b) { return a * b; }, "aa:mul", TOL_AA);
        }
    }
    SECTION("div") {
        for (int i = 0; i < N_ITER; ++i) {
            I xi = random_interval(), yi = nz_iv();
            AF x(ctx, xi), y(ctx, yi);
            check_aa_binary_sound(xi, yi, x / y,
                [](T a, T b) { return a / b; }, "aa:div", TOL_AA);
        }
    }
}

TEST_CASE("fuzz: affine self-composition (dependency)", "[fuzz]")
{
    pappus::affine_context ctx;

    SECTION("x + x") {
        for (int i = 0; i < N_ITER; ++i) {
            I xi = random_interval();
            AF x(ctx, xi);
            I result = (x + x).to_interval();
            for (T v : sample_points(xi, N_SAMPLE))
                if (!result.contains(v + v))
                    FAIL_CHECK("aa x+x: " << v << "+" << v << "=" << (v+v)
                               << " not in [" << result.inf() << "," << result.sup() << "]");
        }
    }
    SECTION("x - x") {
        for (int i = 0; i < N_ITER; ++i) {
            I xi = random_interval();
            AF x(ctx, xi);
            I result = (x - x).to_interval();
            for (T v : sample_points(xi, N_SAMPLE))
                if (!result.contains(T(0)))
                    FAIL_CHECK("aa x-x: 0 not in [" << result.inf() << "," << result.sup() << "]");
        }
    }
    SECTION("x * x") {
        for (int i = 0; i < N_ITER; ++i) {
            I xi = random_interval();
            AF x(ctx, xi);
            I result = (x * x).to_interval();
            for (T v : sample_points(xi, N_SAMPLE))
                if (!result.contains(v * v))
                    FAIL_CHECK("aa x*x: " << v << "^2=" << (v*v)
                               << " not in [" << result.inf() << "," << result.sup() << "]");
        }
    }
    SECTION("x / x") {
        for (int i = 0; i < N_ITER; ++i) {
            I xi = nz_iv();
            AF x(ctx, xi);
            I result = (x / x).to_interval();
            for (T v : sample_points(xi, N_SAMPLE))
                if (!result.contains(T(1)))
                    FAIL_CHECK("aa x/x: 1 not in [" << result.inf() << "," << result.sup() << "]");
        }
    }
}

TEST_CASE("fuzz: affine unary ops", "[fuzz]")
{
    pappus::affine_context ctx;

#define FUZZ_AA(sec, gen, method, sfn, tol)                                   \
    SECTION(sec) {                                                            \
        for (int i = 0; i < N_ITER; ++i) {                                    \
            I xi = gen();                                                     \
            AF x(ctx, xi);                                                    \
            check_aa_unary_sound(xi, x.method(),                              \
                [](T v) { return sfn; }, "aa:" sec, tol);                     \
        }                                                                     \
    }

    FUZZ_AA("neg",     random_interval, operator-,  -v,               TOL_EXACT)
    FUZZ_AA("abs",     random_interval, abs,        std::fabs(v),     TOL_AA)
    FUZZ_AA("floor",   random_interval, floor,      std::floor(v),    TOL_AA)
    FUZZ_AA("ceil",    random_interval, ceil,       std::ceil(v),     TOL_AA)
    FUZZ_AA("cbrt",    random_interval, cbrt,       std::cbrt(v),     TOL_AA)
    FUZZ_AA("exp",     bounded_iv_50,   exp,        std::exp(v),      TOL_AA)
    FUZZ_AA("log",     pos_iv,          log,        std::log(v),      TOL_AA)
    FUZZ_AA("log1p",   above_neg1_iv,   log1p,      std::log1p(v),   TOL_AA)
    FUZZ_AA("sqrt",    nneg_iv,         sqrt,       std::sqrt(v),     TOL_AA)
    FUZZ_AA("isqrt",   pos_iv,          isqrt,      T(1)/std::sqrt(v), TOL_AA)
    FUZZ_AA("inv",     nz_iv,           inv,        T(1) / v,         TOL_AA)
    FUZZ_AA("sin",     random_interval, sin,        std::sin(v),      TOL_AA)
    FUZZ_AA("cos",     random_interval, cos,        std::cos(v),      TOL_AA)
    FUZZ_AA("tan",     tan_iv,          tan,        std::tan(v),      TOL_AA)
    FUZZ_AA("asin",    unit_iv,         asin,       std::asin(v),     TOL_AA)
    FUZZ_AA("acos",    unit_iv,         acos,       std::acos(v),     TOL_AA)
    FUZZ_AA("atan",    random_interval, atan,       std::atan(v),     TOL_AA)
    FUZZ_AA("sinh",    bounded_iv_50,   sinh,       std::sinh(v),     TOL_AA)
    FUZZ_AA("cosh",    bounded_iv_50,   cosh,       std::cosh(v),     TOL_AA)
    FUZZ_AA("tanh",    random_interval, tanh,       std::tanh(v),     TOL_AA)

#undef FUZZ_AA
}

TEST_CASE("fuzz: affine power", "[fuzz]")
{
    pappus::affine_context ctx;

    SECTION("pow(int)") {
        for (int i = 0; i < N_ITER; ++i) {
            int p = irand(0, 8);
            I xi = random_interval();
            if (p == 0) continue;
            AF x(ctx, xi);
            check_aa_unary_sound(xi, x.pow(p),
                [p](T v) { return std::pow(v, p); }, "aa:pow(int)", TOL_AA);
        }
    }
    SECTION("pow(T) positive base") {
        for (int i = 0; i < N_ITER; ++i) {
            T p = urand(0.1, 5.0);
            I xi = pos_iv();
            AF x(ctx, xi);
            check_aa_unary_sound(xi, x.pow(p),
                [p](T v) { return std::pow(v, p); }, "aa:pow(T)", TOL_AA);
        }
    }
}

TEST_CASE("fuzz: affine composite ops", "[fuzz]")
{
    pappus::affine_context ctx;

    SECTION("aq") {
        for (int i = 0; i < N_ITER; ++i) {
            I xi = random_interval(), yi = random_interval();
            AF x(ctx, xi), y(ctx, yi);
            check_aa_binary_sound(xi, yi, pappus::ops::aq(x, y),
                [](T a, T b) { return a / std::sqrt(1 + b * b); }, "aa:aq", TOL_AA);
        }
    }
    SECTION("sqrtabs") {
        for (int i = 0; i < N_ITER; ++i) {
            I xi = random_interval();
            AF x(ctx, xi);
            check_aa_unary_sound(xi, pappus::ops::sqrtabs(x),
                [](T v) { return std::sqrt(std::fabs(v)); }, "aa:sqrtabs", TOL_AA);
        }
    }
    SECTION("logabs") {
        for (int i = 0; i < N_ITER; ++i) {
            I xi = nz_iv();
            AF x(ctx, xi);
            check_aa_unary_sound(xi, pappus::ops::logabs(x),
                [](T v) { return std::log(std::fabs(v)); }, "aa:logabs", TOL_AA);
        }
    }
}
