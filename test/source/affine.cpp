#include <catch2/catch_all.hpp>

#include "aa.h"
#include "pappus/pappus.hpp"


using af = pappus::affine_form<double>;
using ai = pappus::interval<double>;

constexpr double eps = std::numeric_limits<double>::epsilon();

// Relative tolerance: directed rounding in pappus min()/max() can shift
// bounds by 1 ULP, which exceeds eps for values with magnitude > 1.
bool operator==(ai const& lhs, AAInterval const& rhs)
{
    auto tol = [](double a, double b) {
        return eps * 32.0 * std::max({1.0, std::fabs(a), std::fabs(b)});
    };
    return std::fabs(lhs.inf() - rhs.getlo()) < tol(lhs.inf(), rhs.getlo())
        && std::fabs(lhs.sup() - rhs.gethi()) < tol(lhs.sup(), rhs.gethi());
}

bool operator==(AAInterval const& lhs, ai const& rhs)
{
    return rhs == lhs;
}

bool operator==(af const& lhs, AAF const& rhs)
{
    if (lhs.length() > rhs.getlength()) {
        return false;
    }

    for (size_t i = 0; i < lhs.length(); ++i) {
        if (lhs[i] != rhs[i]) {
            return false;
        }
    }
    return true;
}

// Domain errors return af::invalid() (NaN-centered) rather than throwing --
// see affine.hpp's affine_form::invalid() for why.
bool is_invalid(af const& x)
{
    auto iv = x.to_interval();
    return !std::isfinite(iv.inf()) || !std::isfinite(iv.sup());
}

bool operator==(AAF const& lhs, af const& rhs)
{
    return rhs == lhs;
}

// pappus deliberately widens its unary-op bounds by a small (a few ULP)
// safety margin to guarantee soundness against RN-computed endpoint/
// critical-point approximations (see affine_form::apply_unary's doc
// comment) -- aaflib carries no equivalent margin, so pappus's bound is
// expected to be a (very slightly) wider superset of aaflib's, not
// bit-identical to it. The right cross-check for those cases is
// containment, not equality; `slop` absorbs aaflib's own last-bit
// rounding, not pappus's deliberate margin (which is allowed to be
// arbitrarily larger on either side).
bool SoundlyContains(ai const& bound, AAInterval const& ref)
{
    auto slop = eps * 32.0 * std::max({1.0, std::fabs(ref.getlo()), std::fabs(ref.gethi())});
    return bound.inf() <= ref.getlo() + slop && ref.gethi() - slop <= bound.sup();
}


TEST_CASE("affine_form::operator+(af)")
{
    // pappus
    ai u1(-2, 3);
    ai v1(-1, 1);

    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1(ctx, v1);

    af z1 = x1 + y1;

    // aaflib
    AAInterval u2(u1.inf(), u1.sup());
    AAInterval v2(v1.inf(), v1.sup());

    AAF x2(u2);
    AAF y2(v2);
    AAF z2 = x2 + y2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("affine_form::operator+(double, af)")
{
    ai u1(-2, 3);
    double v = 2;

    pappus::affine_context ctx;
    af x1(ctx, u1);

    af z1 = v + x1;

    AAInterval u2(u1.inf(), u1.sup());
    AAF x2(u2);
    AAF z2 = v + x2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("affine_form::operator+=(af)")
{
    // pappus
    ai u1(-2, 3);
    ai v1(-1, 1);

    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1(ctx, v1);

    af z1 = x1 + y1;
    x1 += y1;

    CHECK(z1.to_interval() == x1.to_interval());
}

TEST_CASE("affine_form::operator-(af)")
{
    // pappus
    ai u1(-2, 3);
    ai v1(-1, 1);

    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1(ctx, v1);

    af z1 = x1 - y1;

    // aaflib
    AAInterval u2(u1.inf(), u1.sup());
    AAInterval v2(v1.inf(), v1.sup());

    AAF x2(u2);
    AAF y2(v2);
    AAF z2 = x2 - y2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("affine_form::operator-(double, af)")
{
    ai u1(-2, 3);
    double v = 2;

    pappus::affine_context ctx;
    af x1(ctx, u1);

    af z1 = v - x1;

    AAInterval u2(u1.inf(), u1.sup());
    AAF x2(u2);
    AAF z2 = v - x2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("affine_form::operator-=(af)")
{
    // pappus
    ai u1(-2, 3);
    ai v1(-1, 1);

    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1(ctx, v1);

    af z1 = x1 - y1;
    x1 -= y1;

    CHECK(z1.to_interval() == x1.to_interval());
}

TEST_CASE("affine_form::operator*(af)")
{
    ai u1(-2, 3);
    ai v1(-1, 1);

    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1(ctx, v1);

    af z1 = x1 * y1;

    AAInterval u2(u1.inf(), u1.sup());
    AAInterval v2(v1.inf(), v1.sup());

    AAF x2(u2);
    AAF y2(v2);
    AAF z2 = x2 * y2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("affine_form::operator*(double, af)")
{
    ai u1(-2, 3);
    double v = 2;

    pappus::affine_context ctx;
    af x1(ctx, u1);

    af z1 = v * x1;

    AAInterval u2(u1.inf(), u1.sup());
    AAF x2(u2);
    AAF z2 = v * x2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("affine_form::operator*=(af)")
{
    // pappus
    ai u1(-2, 3);
    ai v1(-1, 1);

    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1(ctx, v1);

    af z1 = x1 * y1;
    x1 *= y1;

    CHECK(z1.to_interval() == x1.to_interval());
}

TEST_CASE("affine_form::operator/(af)")
{
    ai u1(2, 3);
    ai v1(1.5, 2);

    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1(ctx, v1);

    af z1 = x1 / y1;

    AAInterval u2(u1.inf(), u1.sup());
    AAInterval v2(v1.inf(), v1.sup());

    AAF x2(u2);
    AAF y2(v2);
    AAF z2 = x2 / y2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("affine_form::operator/(double, af)")
{
    ai u1(2, 3);
    double v = 2;

    pappus::affine_context ctx;
    af x1(ctx, u1);

    af z1 = v / x1;

    AAInterval u2(u1.inf(), u1.sup());
    AAF x2(u2);
    AAF z2 = v / x2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("affine_form::inv()")
{
    ai u1(1, 4);
    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1 = x1.inv();

    AAInterval u2(u1.inf(), u1.sup());
    AAF x2(u2);
    AAF y2 = inv(x2);

    CHECK(y1.to_interval() == y2.convert());
}

TEST_CASE("affine_form::inv() rejects intervals containing zero")
{
    pappus::affine_context ctx;
    CHECK(is_invalid(af(ctx, ai(-1.0, 1.0)).inv()));
    CHECK(is_invalid(af(ctx, 0.0).inv()));
}

TEST_CASE("affine_form::pow(int)")
{
    ai u1(1, 4);
    int exponent = 3;
    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1 = x1.pow(exponent);

    AAInterval u2(u1.inf(), u1.sup());
    AAF x2(u2);
    AAF y2 = x2 ^ exponent;

    CHECK(SoundlyContains(y1.to_interval(), y2.convert()));
}

TEST_CASE("affine_form::pow(double)")
{
    ai u1(1, 4);
    pappus::affine_context ctx;
    af x1(ctx, u1);

    AAInterval u2(u1.inf(), u1.sup());
    AAF x2(u2);

    for (auto exponent = 0.5; exponent < 5; exponent += 0.5) {
        af y1 = x1.pow(exponent);
        AAF y2 = x2 ^ exponent;
        // Containment, not equality/structural-term comparison: pow(double)
        // routes through apply_unary_bounded, which now adds a deliberate
        // small safety-margin noise term aaflib has no equivalent of (see
        // SoundlyContains's doc comment) -- so both the interval bound and
        // the raw term count are expected to differ from aaflib's.
        CHECK(SoundlyContains(y1.to_interval(), y2.convert()));
    }

    af y1 = x1.pow(0.0);
    AAF y2 = x2 ^ 0.0;
    CHECK(SoundlyContains(y1.to_interval(), y2.convert()));
}

TEST_CASE("affine_form::pow(double) domain checks")
{
    pappus::affine_context ctx;

    CHECK(is_invalid(af(ctx, ai(-1.0, 4.0)).pow(0.5)));
    CHECK(is_invalid(af(ctx, ai(0.0, 4.0)).pow(-0.5)));
}

TEST_CASE("affine_form::last_index() on constant form")
{
    pappus::affine_context ctx;
    af x(ctx, 1.0);
    CHECK_THROWS_AS(x.last_index(), std::logic_error);
}

TEST_CASE("affine_form point interval stays constant")
{
    pappus::affine_context ctx;
    af x(ctx, ai(2.0, 2.0));
    CHECK(x.length() == 0);
    CHECK(x.radius() == 0.0);

    // pow(T) has no tagged directed-rounding primitive (there is no
    // fp::op_pow), so its constant path widens the RN result by 1 ULP
    // each direction via widen_scalar_result -- exactly like pappus's own
    // pre-existing interval<T>::cosh()/exp()/cos()/cbrt() already do
    // unconditionally for their own constant/point inputs (see
    // fp::detail::outward_lo/outward_hi: only the "_keep_zero" variants
    // preserve an exact result). 2.0^2.0 = 4.0 exactly in real arithmetic,
    // but the result is a thin (length 1) sound enclosure of 4.0, not
    // required to be bit-exact.
    auto y = x.pow(2.0);
    CHECK(y.length() <= 1);
    CHECK(y.to_interval().contains(4.0));

    auto z = af::pow(2.0, x);
    CHECK(z.length() <= 1);
    CHECK(z.to_interval().contains(4.0));
}


TEST_CASE("affine_form rejects mixed contexts")
{
    pappus::affine_context lhs_ctx;
    pappus::affine_context rhs_ctx;
    af x(lhs_ctx, ai(-1.0, 1.0));
    af y(rhs_ctx, ai(2.0, 3.0));

    CHECK_THROWS_AS(x + y, std::invalid_argument);
    CHECK_THROWS_AS(x - y, std::invalid_argument);
    CHECK_THROWS_AS(x * y, std::invalid_argument);
    CHECK_THROWS_AS(x / y, std::invalid_argument);
    CHECK_THROWS_AS(x.pow(y), std::invalid_argument);
}

TEST_CASE("affine_form keeps context state alive")
{
    auto make_form = []() {
        pappus::affine_context ctx;
        return af(ctx, ai(1.0, 2.0));
    };

    auto x = make_form();
    auto y = x + 1.0;
    CHECK(y.to_interval() == ai(2.0, 3.0));
}

TEST_CASE("affine_context approximation mode setter")
{
    pappus::affine_context ctx;
    ctx.set_approximation_mode(pappus::approximation_mode::SECANT);
    CHECK(ctx.approximation_mode() == pappus::approximation_mode::SECANT);
}

TEST_CASE("affine_form::pow(af)")
{
    auto check_exp = [](auto const& base, auto const& exponent) {
        pappus::affine_context ctx;
        af x1(ctx, base);
        af e1(ctx, exponent);
        af y1 = x1.pow(e1);

        AAF x2 = AAF(AAInterval(base.inf(), base.sup()));
        AAF e2 = AAF(AAInterval(exponent.inf(), exponent.sup()));
        AAF y2 = x2 ^ e2;

        auto w1 = y1.to_interval();
        auto w2 = y2.convert();
        CHECK(Catch::Approx(w1.inf()) == w2.getlo());
        CHECK(Catch::Approx(w1.sup()) == w2.gethi());
    };

    SECTION("[1, 4] ^ [0, 0.5]") { check_exp(ai(1.0, 4.0), ai(0.0, 0.5)); }
    SECTION("[0.1, 2] ^ [-1, 0]") { check_exp(ai(0.1, 2.0), ai(-1.0, 0.0)); }
}

TEST_CASE("affine_form::pow(double, af)")
{
    auto check_exp = [](double base, auto const& exponent) {
        pappus::affine_context ctx;
        af e1(ctx, exponent);
        af y1 = af::pow(base, e1);

        AAF e2 = AAF(AAInterval(exponent.inf(), exponent.sup()));
        AAF y2 = aaf_pow(base, e2);
        CHECK(y1.to_interval() == y2.convert());
        CHECK(y1 == y2);
    };

    check_exp(0.4, ai(1.0, 2.0));
    check_exp(3.14, ai(0.0, 0.5));
    check_exp(2, ai(2, 3));
}

TEST_CASE("affine_form::sqrt()")
{
    ai u1(1.0, 4.0);
    pappus::affine_context ctx;
    af x1(ctx, u1);
    auto y1 = x1.sqrt();

    AAInterval u2(u1.inf(), u1.sup());
    AAF x2(u2);
    auto y2 = sqrt(x2);

    CHECK(y1.to_interval() == y2.convert());
}

TEST_CASE("affine_form::isqrt()")
{
    ai u1(2.0, 7.0);
    pappus::affine_context ctx;
    af x1(ctx, u1);
    auto y1 = x1.isqrt();

    AAInterval u2(u1.inf(), u1.sup());
    AAF x2(u2);
    auto y2 = isqrt(x2);

    CHECK(y1.to_interval() == y2.convert());
}

/******************************************************
 * Transcendental function tests                      *
 *****************************************************/

// Returns true if result.to_interval() contains f(x) (to within 4 ULP) for 200 sampled x.
// The 4-ULP slack accounts for rounding in the Chebyshev coefficient computation — the
// approximation is mathematically exact at the endpoints but FP rounding can push the
// computed bound by 1-2 ULP relative to f(endpoint).
template<typename Scalar>
bool sound(pappus::affine_form<double> const& result, Scalar f, ai const& range, int n = 200)
{
    ai ri = result.to_interval();
    double lo = range.inf(), hi = range.sup(), step = (hi - lo) / n;
    constexpr double eps4 = 4 * std::numeric_limits<double>::epsilon();
    double slack = eps4 * std::max({1.0, std::fabs(ri.inf()), std::fabs(ri.sup())});
    for (int i = 0; i <= n; ++i) {
        double y = f(lo + i * step);
        if (y < ri.inf() - slack || y > ri.sup() + slack) return false;
    }
    return true;
}

TEST_CASE("affine_form::exp()")
{
    pappus::affine_context ctx;

    // constant form
    // op_exp has no "_keep_zero" variant -- unconditional 1-ULP nudge,
    // same as pappus's pre-existing interval<T>::exp().
    CHECK(af(ctx, 1.0).exp().to_interval().contains(std::exp(1.0)));

    ai u(0.5, 1.5);
    af x(ctx, u);
    auto y = x.exp();
    CHECK(sound(y, [](double x) { return std::exp(x); }, u));
    CHECK(y.to_interval().inf() <= std::exp(0.5));
    CHECK(y.to_interval().sup() >= std::exp(1.5));

    ai u2(-1.0, 1.0);
    af x2(ctx, u2);
    CHECK(sound(x2.exp(), [](double x) { return std::exp(x); }, u2));
}

TEST_CASE("affine_form::log()")
{
    pappus::affine_context ctx;

    CHECK(af(ctx, 1.0).log().to_interval() == ai(0.0));

    ai u(1.0, 4.0);
    af x(ctx, u);
    auto y = x.log();
    CHECK(sound(y, [](double x) { return std::log(x); }, u));
    CHECK(y.to_interval().sup() >= std::log(4.0));

    CHECK(is_invalid(af(ctx, ai(0.0, 1.0)).log()));
    CHECK(is_invalid(af(ctx, -1.0).log()));
}

TEST_CASE("affine_form::sin()")
{
    pappus::affine_context ctx;

    CHECK(af(ctx, 0.0).sin().to_interval() == ai(0.0));

    // Monotone piece: sin is decreasing on [π/2, π]
    ai u(pappus::fp::half_pi_v<double>, pappus::fp::pi_v<double>);
    af x(ctx, u);
    auto y = x.sin();
    CHECK(sound(y, [](double x) { return std::sin(x); }, u));

    // Non-monotone: sin on [0, π] has max at π/2
    ai u2(0.0, pappus::fp::pi_v<double>);
    af x2(ctx, u2);
    auto y2 = x2.sin();
    CHECK(sound(y2, [](double x) { return std::sin(x); }, u2));

    // Full period → returns [-1, 1]
    ai u3(0.0, pappus::fp::two_pi_v<double> + 0.1);
    af x3(ctx, u3);
    auto y3 = x3.sin();
    CHECK(y3.to_interval().inf() <= -1.0);
    CHECK(y3.to_interval().sup() >= 1.0);
}

TEST_CASE("affine_form::cos()")
{
    pappus::affine_context ctx;

    // op_cos has no "_keep_zero" variant -- unconditional 1-ULP nudge,
    // same as pappus's pre-existing interval<T>::cos().
    CHECK(af(ctx, 0.0).cos().to_interval().contains(1.0));

    // Monotone decreasing: cos on [0, π]
    ai u(0.0, pappus::fp::pi_v<double>);
    af x(ctx, u);
    auto y = x.cos();
    CHECK(sound(y, [](double x) { return std::cos(x); }, u));

    // Contains max at 0 and min at π
    CHECK(y.to_interval().sup() >= 1.0);
    CHECK(y.to_interval().inf() <= -1.0);

    ai u2(-0.5, 0.5);
    af x2(ctx, u2);
    CHECK(sound(x2.cos(), [](double x) { return std::cos(x); }, u2));
}

TEST_CASE("affine_form::tan()")
{
    pappus::affine_context ctx;

    CHECK(af(ctx, 0.0).tan().to_interval() == ai(0.0));

    ai u(-0.5, 0.5);
    af x(ctx, u);
    auto y = x.tan();
    CHECK(sound(y, [](double x) { return std::tan(x); }, u));

    // Crosses asymptote at π/2
    CHECK(is_invalid(af(ctx, ai(0.0, 2.0)).tan()));
}

TEST_CASE("affine_form::sinh()")
{
    pappus::affine_context ctx;

    CHECK(af(ctx, 0.0).sinh().to_interval() == ai(0.0));

    ai u(-1.0, 1.0);
    af x(ctx, u);
    auto y = x.sinh();
    CHECK(sound(y, [](double x) { return std::sinh(x); }, u));

    ai u2(0.5, 2.0);
    af x2(ctx, u2);
    CHECK(sound(x2.sinh(), [](double x) { return std::sinh(x); }, u2));

    ai u3(-2.0, -0.5);
    af x3(ctx, u3);
    CHECK(sound(x3.sinh(), [](double x) { return std::sinh(x); }, u3));
}

TEST_CASE("affine_form::cosh()")
{
    pappus::affine_context ctx;

    // op_cosh has no "_keep_zero" directed-rounding variant (matches
    // pappus's pre-existing interval<T>::cosh(), which has the same
    // unconditional 1-ULP nudge) -- cosh(0)=1 exactly in real arithmetic,
    // but the sound enclosure is thin, not required to be bit-exact [1,1].
    CHECK(af(ctx, 0.0).cosh().to_interval().contains(1.0));

    // Symmetric, minimum at 0
    ai u(-1.0, 1.0);
    af x(ctx, u);
    auto y = x.cosh();
    CHECK(sound(y, [](double x) { return std::cosh(x); }, u));

    // Positive-only
    ai u2(0.5, 2.0);
    af x2(ctx, u2);
    CHECK(sound(x2.cosh(), [](double x) { return std::cosh(x); }, u2));

    // Negative-only
    ai u3(-2.0, -0.5);
    af x3(ctx, u3);
    CHECK(sound(x3.cosh(), [](double x) { return std::cosh(x); }, u3));
}

TEST_CASE("affine_form::tanh()")
{
    pappus::affine_context ctx;

    CHECK(af(ctx, 0.0).tanh().to_interval() == ai(0.0));

    ai u(-1.0, 1.0);
    af x(ctx, u);
    auto y = x.tanh();
    CHECK(sound(y, [](double x) { return std::tanh(x); }, u));
    CHECK(y.to_interval().inf() >= -1.0);
    CHECK(y.to_interval().sup() <=  1.0);

    ai u2(0.5, 2.0);
    af x2(ctx, u2);
    CHECK(sound(x2.tanh(), [](double x) { return std::tanh(x); }, u2));

    ai u3(-2.0, -0.5);
    af x3(ctx, u3);
    CHECK(sound(x3.tanh(), [](double x) { return std::tanh(x); }, u3));
}

TEST_CASE("affine_form::asin()")
{
    pappus::affine_context ctx;

    CHECK(af(ctx, 0.0).asin().to_interval() == ai(0.0));
    // outward_lo/hi_keep_zero only preserves exactness when the RESULT is
    // exactly zero -- asin(1)=half_pi is nonzero, so it gets the ordinary
    // 1-ULP outward nudge each direction.
    CHECK(af(ctx, 1.0).asin().to_interval().contains(pappus::fp::half_pi_v<double>));

    ai u(-0.8, 0.8);
    af x(ctx, u);
    auto y = x.asin();
    CHECK(sound(y, [](double x) { return std::asin(x); }, u));

    CHECK(is_invalid(af(ctx, ai(-2.0, 0.0)).asin()));
}

TEST_CASE("affine_form::acos()")
{
    pappus::affine_context ctx;

    CHECK(af(ctx, 1.0).acos().to_interval() == ai(0.0));
    // outward_lo/hi_keep_zero only preserves exactness when the RESULT is
    // exactly zero (see fp/math.hpp) -- acos(0)=half_pi is a nonzero
    // result, so it gets the ordinary 1-ULP outward nudge each direction.
    CHECK(af(ctx, 0.0).acos().to_interval().contains(pappus::fp::half_pi_v<double>));

    ai u(-0.8, 0.8);
    af x(ctx, u);
    auto y = x.acos();
    CHECK(sound(y, [](double x) { return std::acos(x); }, u));

    CHECK(is_invalid(af(ctx, ai(0.0, 2.0)).acos()));
}

TEST_CASE("affine_form::atan()")
{
    pappus::affine_context ctx;

    CHECK(af(ctx, 0.0).atan().to_interval() == ai(0.0));

    ai u(-2.0, 2.0);
    af x(ctx, u);
    auto y = x.atan();
    CHECK(sound(y, [](double x) { return std::atan(x); }, u));
    CHECK(y.to_interval().inf() >= -pappus::fp::half_pi_v<double>);
    CHECK(y.to_interval().sup() <=  pappus::fp::half_pi_v<double>);

    ai u2(0.5, 3.0);
    af x2(ctx, u2);
    CHECK(sound(x2.atan(), [](double x) { return std::atan(x); }, u2));

    ai u3(-3.0, -0.5);
    af x3(ctx, u3);
    CHECK(sound(x3.atan(), [](double x) { return std::atan(x); }, u3));
}

/******************************************************
 * Arbitrary expressions tests                        *
 *****************************************************/
TEST_CASE("X^2 + X")
{
    SECTION("Dependency problem")
    {
        ai u1(-1, 1);
        pappus::affine_context ctx;

        af x1(ctx, u1);

        auto z1 = x1 * x1 + x1;

        AAInterval u2(u1.inf(), u1.sup());
        AAF x2(u2);
        AAF z2 = x2 * x2 + x2;

        CHECK(z1.to_interval() == z2.convert());
    }

    SECTION("Rewrite")
    {
        pappus::affine_context ctx;
        ai u1(-1, 1);
        af x1(ctx, u1);
        auto z1 = (x1 + 0.5).pow(2.0) - 0.25;

        AAInterval u2(u1.inf(), u1.sup());
        AAF x2(u2);
        AAF z2 = ((x2 + 0.5) ^ 2) - 0.25;

        CHECK(z1.to_interval() == z2.convert());
    }
}

/******************************************************
 * Domain predicates                                  *
 *****************************************************/

TEST_CASE("domain predicates: log")
{
    CHECK(pappus::log_domain_ok(ai(0.5, 2.0)));
    CHECK(pappus::log_domain_ok(ai(1.0, 1.0)));
    CHECK_FALSE(pappus::log_domain_ok(ai(0.0, 2.0)));   // inf == 0
    CHECK_FALSE(pappus::log_domain_ok(ai(-1.0, 2.0)));  // inf < 0
    CHECK_FALSE(pappus::log_domain_ok(ai(-2.0, -1.0))); // all negative
    CHECK_FALSE(pappus::log_domain_ok(ai::empty()));

    pappus::affine_context ctx;
    CHECK(pappus::log_domain_ok(af(ctx, ai(1.0, 4.0))));
    CHECK_FALSE(pappus::log_domain_ok(af(ctx, ai(0.0, 4.0))));
}

TEST_CASE("domain predicates: sqrt")
{
    CHECK(pappus::sqrt_domain_ok(ai(0.0, 2.0)));  // sqrt(0) valid
    CHECK(pappus::sqrt_domain_ok(ai(1.0, 4.0)));
    CHECK_FALSE(pappus::sqrt_domain_ok(ai(-0.5, 2.0)));
    CHECK_FALSE(pappus::sqrt_domain_ok(ai(-2.0, -1.0)));
    CHECK_FALSE(pappus::sqrt_domain_ok(ai::empty()));
}

TEST_CASE("domain predicates: isqrt")
{
    CHECK(pappus::isqrt_domain_ok(ai(1.0, 4.0)));
    CHECK_FALSE(pappus::isqrt_domain_ok(ai(0.0, 4.0)));  // isqrt(0) = inf
    CHECK_FALSE(pappus::isqrt_domain_ok(ai(-1.0, 4.0)));
    CHECK_FALSE(pappus::isqrt_domain_ok(ai::empty()));
}

TEST_CASE("domain predicates: inv")
{
    CHECK(pappus::inv_domain_ok(ai(1.0, 2.0)));
    CHECK(pappus::inv_domain_ok(ai(-2.0, -1.0)));
    CHECK_FALSE(pappus::inv_domain_ok(ai(-1.0, 1.0)));  // straddles 0
    CHECK_FALSE(pappus::inv_domain_ok(ai(0.0, 1.0)));   // inf == 0
    CHECK_FALSE(pappus::inv_domain_ok(ai(-1.0, 0.0)));  // sup == 0
    CHECK_FALSE(pappus::inv_domain_ok(ai::empty()));
}

TEST_CASE("domain predicates: tan")
{
    CHECK(pappus::tan_domain_ok(ai(0.0, 1.0)));
    CHECK(pappus::tan_domain_ok(ai(-0.5, 0.5)));
    CHECK(pappus::tan_domain_ok(ai(pappus::fp::pi_v<double> + 0.1, 1.5 * pappus::fp::pi_v<double> - 0.1)));
    CHECK_FALSE(pappus::tan_domain_ok(ai(0.0, 2.0)));   // crosses π/2
    CHECK_FALSE(pappus::tan_domain_ok(ai::empty()));
}

TEST_CASE("domain predicates: asin / acos")
{
    CHECK(pappus::asin_domain_ok(ai(-1.0, 1.0)));
    CHECK(pappus::asin_domain_ok(ai(-0.5, 0.5)));
    CHECK_FALSE(pappus::asin_domain_ok(ai(-1.5, 0.5)));
    CHECK_FALSE(pappus::asin_domain_ok(ai(-0.5, 1.5)));

    CHECK(pappus::acos_domain_ok(ai(-1.0, 1.0)));
    CHECK_FALSE(pappus::acos_domain_ok(ai(-2.0, 0.0)));
}

TEST_CASE("domain predicates: pow")
{
    CHECK(pappus::pow_domain_ok(ai(-2.0, 2.0), 2.0));   // integer exp: any base
    CHECK(pappus::pow_domain_ok(ai(-2.0, 2.0), 3.0));
    CHECK_FALSE(pappus::pow_domain_ok(ai(-2.0, 2.0), 0.5));  // fractional: needs inf >= 0
    CHECK(pappus::pow_domain_ok(ai(0.0, 2.0), 0.5));
    CHECK(pappus::pow_domain_ok(ai(1.0, 2.0), -1.0));   // negative exp: needs inf > 0
    CHECK_FALSE(pappus::pow_domain_ok(ai(0.0, 2.0), -1.0)); // inf == 0
    CHECK_FALSE(pappus::pow_domain_ok(ai(-1.0, 2.0), -1.0));
}

/******************************************************
 * try_* wrappers                                     *
 *****************************************************/

TEST_CASE("try_ wrappers: valid domain returns value matching throwing version")
{
    pappus::affine_context ctx;

    auto check = [&](auto try_fn, auto member_fn, ai const& u) {
        af x(ctx, u);
        auto opt = try_fn(x);
        REQUIRE(opt.has_value());
        CHECK(opt->to_interval() == (x.*member_fn)().to_interval());
    };

    check(pappus::try_log<double>,   &af::log,   ai(1.0, 4.0));
    check(pappus::try_log1p<double>, &af::log1p, ai(0.0, 3.0));
    check(pappus::try_sqrt<double>,  &af::sqrt,  ai(1.0, 4.0));
    check(pappus::try_isqrt<double>, &af::isqrt, ai(1.0, 4.0));
    check(pappus::try_inv<double>,   &af::inv,   ai(1.0, 4.0));
    check(pappus::try_tan<double>,   &af::tan,   ai(-0.5, 0.5));
    check(pappus::try_asin<double>,  &af::asin,  ai(-0.8, 0.8));
    check(pappus::try_acos<double>,  &af::acos,  ai(-0.8, 0.8));
}

TEST_CASE("try_ wrappers: invalid domain returns nullopt")
{
    pappus::affine_context ctx;

    CHECK_FALSE(pappus::try_log  (af(ctx, ai(-1.0, 2.0))).has_value()); // inf <= 0
    CHECK_FALSE(pappus::try_log  (af(ctx, ai(-2.0, -1.0))).has_value());
    CHECK_FALSE(pappus::try_log1p(af(ctx, ai(-2.0, 0.0))).has_value());  // inf <= -1
    CHECK_FALSE(pappus::try_sqrt (af(ctx, ai(-1.0, 2.0))).has_value()); // inf < 0
    CHECK_FALSE(pappus::try_isqrt(af(ctx, ai(0.0, 4.0))).has_value());  // inf == 0
    CHECK_FALSE(pappus::try_inv  (af(ctx, ai(-1.0, 1.0))).has_value()); // straddles 0
    CHECK_FALSE(pappus::try_tan  (af(ctx, ai(0.0, 2.0))).has_value());  // crosses π/2
    CHECK_FALSE(pappus::try_asin (af(ctx, ai(-2.0, 0.5))).has_value());
    CHECK_FALSE(pappus::try_pow  (af(ctx, ai(-1.0, 2.0)), 0.5).has_value());
}

/******************************************************
 * safe_* wrappers                                    *
 *****************************************************/

// Check that every sample point in the valid subdomain is enclosed.
template<typename Scalar>
bool safe_sound(pappus::affine_form<double> const& result, Scalar f,
                ai const& valid_range, int n = 200)
{
    ai ri = result.to_interval();
    double lo = valid_range.inf(), hi = valid_range.sup(), step = (hi - lo) / n;
    constexpr double eps4 = 4 * std::numeric_limits<double>::epsilon();
    double slack = eps4 * std::max({1.0, std::fabs(ri.inf()), std::fabs(ri.sup())});
    for (int i = 0; i <= n; ++i) {
        double y = f(lo + i * step);
        if (!std::isfinite(y)) continue;
        if (y < ri.inf() - slack || y > ri.sup() + slack) return false;
    }
    return true;
}

TEST_CASE("safe_log")
{
    pappus::affine_context ctx;

    SECTION("entirely valid — same as log()")
    {
        af x(ctx, ai(1.0, 4.0));
        auto y = pappus::safe_log(x);
        REQUIRE(y.has_value());
        CHECK(y->to_interval() == x.log().to_interval());
    }

    SECTION("partial overlap — clamps, result is sound on valid subdomain")
    {
        af x(ctx, ai(-1.0, 4.0));
        auto y = pappus::safe_log(x);
        REQUIRE(y.has_value());
        auto ri = y->to_interval();
        // Upper bound must reach log(4); lower bound must be finite (clamped lo > 0)
        CHECK(ri.sup() >= std::log(4.0));
        CHECK(std::isfinite(ri.inf()));
        // Sound on a conservative inner subdomain well within the clamped range
        CHECK(safe_sound(*y, [](double v) { return std::log(v); }, ai(1e-6, 4.0)));
    }

    SECTION("entirely invalid — nullopt")
    {
        CHECK_FALSE(pappus::safe_log(af(ctx, ai(-2.0, -1.0))).has_value());
        CHECK_FALSE(pappus::safe_log(af(ctx, ai(-1.0, 0.0))).has_value());
    }
}

TEST_CASE("safe_sqrt")
{
    pappus::affine_context ctx;

    SECTION("entirely valid")
    {
        af x(ctx, ai(0.0, 4.0));
        auto y = pappus::safe_sqrt(x);
        REQUIRE(y.has_value());
        CHECK(y->to_interval() == x.sqrt().to_interval());
    }

    SECTION("partial overlap")
    {
        af x(ctx, ai(-2.0, 4.0));
        auto y = pappus::safe_sqrt(x);
        REQUIRE(y.has_value());
        CHECK(safe_sound(*y, [](double v) { return std::sqrt(v); }, ai(0.0, 4.0)));
        CHECK(y->to_interval().sup() >= std::sqrt(4.0));
    }

    SECTION("entirely invalid")
    {
        CHECK_FALSE(pappus::safe_sqrt(af(ctx, ai(-4.0, -1.0))).has_value());
    }
}

TEST_CASE("safe_isqrt")
{
    pappus::affine_context ctx;

    SECTION("entirely valid")
    {
        af x(ctx, ai(1.0, 4.0));
        auto y = pappus::safe_isqrt(x);
        REQUIRE(y.has_value());
        CHECK(y->to_interval() == x.isqrt().to_interval());
    }

    SECTION("partial overlap (includes zero)")
    {
        af x(ctx, ai(0.0, 4.0));
        auto y = pappus::safe_isqrt(x);
        REQUIRE(y.has_value());
        double lo = std::numeric_limits<double>::min();
        CHECK(safe_sound(*y, [](double v) { return 1.0 / std::sqrt(v); }, ai(lo, 4.0)));
    }

    SECTION("entirely invalid")
    {
        CHECK_FALSE(pappus::safe_isqrt(af(ctx, ai(-2.0, 0.0))).has_value());
    }
}

TEST_CASE("safe_asin")
{
    pappus::affine_context ctx;

    SECTION("entirely valid")
    {
        af x(ctx, ai(-0.8, 0.8));
        auto y = pappus::safe_asin(x);
        REQUIRE(y.has_value());
        CHECK(y->to_interval() == x.asin().to_interval());
    }

    SECTION("clamps inf below -1")
    {
        af x(ctx, ai(-2.0, 0.5));
        auto y = pappus::safe_asin(x);
        REQUIRE(y.has_value());
        CHECK(safe_sound(*y, [](double v) { return std::asin(v); }, ai(-1.0, 0.5)));
    }

    SECTION("clamps sup above 1")
    {
        af x(ctx, ai(-0.5, 2.0));
        auto y = pappus::safe_asin(x);
        REQUIRE(y.has_value());
        CHECK(safe_sound(*y, [](double v) { return std::asin(v); }, ai(-0.5, 1.0)));
    }

    SECTION("entirely outside [-1, 1]")
    {
        CHECK_FALSE(pappus::safe_asin(af(ctx, ai(2.0, 3.0))).has_value());
        CHECK_FALSE(pappus::safe_asin(af(ctx, ai(-3.0, -2.0))).has_value());
    }
}

TEST_CASE("safe_acos")
{
    pappus::affine_context ctx;

    SECTION("entirely valid")
    {
        af x(ctx, ai(-0.8, 0.8));
        auto y = pappus::safe_acos(x);
        REQUIRE(y.has_value());
        CHECK(y->to_interval() == x.acos().to_interval());
    }

    SECTION("clamps to [-1, 1]")
    {
        af x(ctx, ai(-2.0, 2.0));
        auto y = pappus::safe_acos(x);
        REQUIRE(y.has_value());
        CHECK(safe_sound(*y, [](double v) { return std::acos(v); }, ai(-1.0, 1.0)));
    }

    SECTION("entirely outside [-1, 1]")
    {
        CHECK_FALSE(pappus::safe_acos(af(ctx, ai(2.0, 3.0))).has_value());
    }
}

TEST_CASE("safe_log1p")
{
    pappus::affine_context ctx;

    SECTION("entirely valid")
    {
        af x(ctx, ai(0.0, 3.0));
        auto y = pappus::safe_log1p(x);
        REQUIRE(y.has_value());
        CHECK(y->to_interval() == x.log1p().to_interval());
    }

    SECTION("partial overlap")
    {
        af x(ctx, ai(-2.0, 3.0));
        auto y = pappus::safe_log1p(x);
        REQUIRE(y.has_value());
        auto ri = y->to_interval();
        CHECK(ri.sup() >= std::log1p(3.0));
        CHECK(std::isfinite(ri.inf()));
        CHECK(safe_sound(*y, [](double v) { return std::log1p(v); }, ai(-0.9, 3.0)));
    }

    SECTION("entirely invalid")
    {
        CHECK_FALSE(pappus::safe_log1p(af(ctx, ai(-3.0, -1.0))).has_value());
    }
}

TEST_CASE("affine_form::log1p encloses its interior maximum residual")
{
    // On [0, 3], log1p's Chebyshev tangent point lies strictly inside the
    // domain. Endpoint-only checks miss a misplaced tangent residual.
    pappus::affine_context ctx;
    ai const input{0.0, 3.0};
    af const x{ctx, input};
    auto const result = x.log1p().to_interval();

    auto const alpha = std::log1p(input.sup()) / input.sup();
    auto const critical = 1.0 / alpha - 1.0;
    REQUIRE(critical > input.inf());
    REQUIRE(critical < input.sup());
    CHECK(result.contains(std::log1p(critical)));
}

// --- Regression tests: 2026-09 affine-arithmetic soundness review ---

TEST_CASE("affine_form::pow(int) rejects negative exponent crossing zero")
{
    pappus::affine_context ctx;
    // x in [-1, 1], x^-3 is unbounded/undefined at the interior point x=0
    // -- pow(int) used to have no guard at all for exponent < -1 (only
    // exponent == -1 delegated to inv(), which does check this).
    CHECK(is_invalid(af(ctx, ai(-1.0, 1.0)).pow(-3)));
    CHECK(is_invalid(af(ctx, ai(-1.0, 1.0)).pow(-2)));
    // Domain touching zero from one side is equally invalid (mirrors
    // pow(double)'s existing, unmodified guard).
    CHECK(is_invalid(af(ctx, ai(0.0, 2.0)).pow(-2)));
    // Purely negative or purely positive, not touching zero, must still work.
    CHECK_FALSE(is_invalid(af(ctx, ai(1.0, 2.0)).pow(-2)));
    CHECK_FALSE(is_invalid(af(ctx, ai(-2.0, -1.0)).pow(-3)));
}

TEST_CASE("affine_form::pow(int) result is sound against dense sampling")
{
    // pow(int) used to build its result by hand, bypassing the
    // self-correcting endpoint check every other unary op goes through.
    pappus::affine_context ctx;
    ai const input{1.0, 4.0};
    af const x{ctx, input};
    CHECK(sound(x.pow(3), [](double v) { return std::pow(v, 3); }, input));
    CHECK(sound(x.pow(-2), [](double v) { return std::pow(v, -2); }, input));
    CHECK(sound(x.pow(4), [](double v) { return std::pow(v, 4); }, input));
}

TEST_CASE("affine_form::operator*=(T) introduces no spurious slack for exact scaling")
{
    // Regression: an earlier version of this fix used an eps-proportional
    // upper bound on rounding error instead of an exact error-free
    // transform, which added nonzero slack even to EXACT multiplies (e.g.
    // by -1). That silently turned abs()'s exactly-zero-at-boundary case
    // into a small negative value, which then wrongly tripped sqrt()'s
    // `a < 0` domain guard on a mathematically valid (touching-zero) input.
    pappus::affine_context ctx;
    af const x{ctx, ai(-8.65874, -0.0)};
    auto const a = x.abs();
    CHECK(a.min() == 0.0); // exact: negating [-8.65874, -0] must hit 0 exactly, not ~-1.4e-6
    CHECK_FALSE(is_invalid(a.sqrt()));

    // Multiplying by other exact scale factors (powers of two) must also
    // introduce no term growth/slack.
    pappus::affine_context ctx2;
    af const y{ctx2, ai(1.0, 3.0)};
    auto const doubled = y * 2.0;
    CHECK(doubled.length() == y.length());
    CHECK(doubled.to_interval() == ai(2.0, 6.0));
}

TEST_CASE("affine_form::tan() rejects a pole even when it sits at the domain's own lower edge")
{
    // interval<T>::tan()'s old quadrant-comparison guard treated
    // [pi/2, pi] as pole-free; tan is actually unbounded immediately to
    // the right of pi/2. affine_form::tan() uses its own (now-shared)
    // predicate, exercised here on the same case for parity with
    // interval.cpp's equivalent regression test.
    pappus::affine_context ctx;
    constexpr double half_pi = pappus::fp::half_pi_v<double>;
    constexpr double pi = pappus::fp::pi_v<double>;
    CHECK(is_invalid(af(ctx, ai(half_pi, pi)).tan()));
    // A domain entirely within one branch, not touching any pole, must
    // still produce a normal finite result.
    CHECK_FALSE(is_invalid(af(ctx, ai(0.0, 1.0)).tan()));
}

TEST_CASE("affine_form::condense() error term is a sound upper bound, not just RN-close")
{
    // Build a form with several small noise terms of similar magnitude
    // (so the RN transform_reduce sum and the correctly-rounded outward
    // sum can plausibly diverge across many additions), condense it down
    // to fewer terms, and check the resulting bound still soundly encloses
    // the true combined range. NOTE: `after` is not required to be a
    // superset of `before` -- both are independently-computed, independently
    // outward-rounded upper bounds on the same true L1 magnitude (built via
    // different summation groupings), so floating-point summation order can
    // make one a few ULP tighter than the other despite both remaining
    // sound; the real contract is soundness against the true value, not
    // monotonicity between two different valid computations of it.
    pappus::affine_context ctx;
    af acc(ctx, 0.0);
    for (int i = 0; i < 20; ++i) {
        af xi(ctx, ai(-1.0, 1.0));
        acc = acc + xi * 0.1;
    }
    // True combined range: 20 independent [-0.1, 0.1] contributions, worst
    // case all aligned -> exactly [-2, 2].
    ai const truth{-2.0, 2.0};
    CHECK(acc.to_interval().contains(truth));
    acc.condense(4);
    CHECK(acc.length() <= 4);
    CHECK(acc.to_interval().contains(truth));
}


