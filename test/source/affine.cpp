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

bool operator==(AAF const& lhs, af const& rhs)
{
    return rhs == lhs;
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
    CHECK_THROWS_AS(af(ctx, ai(-1.0, 1.0)).inv(), std::invalid_argument);
    CHECK_THROWS_AS(af(ctx, 0.0).inv(), std::invalid_argument);
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

    CHECK(y1.to_interval() == y2.convert());
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
        CHECK(y1.to_interval() == y2.convert());
        CHECK(y1 == y2);
    }

    af y1 = x1.pow(0.0);
    AAF y2 = x2 ^ 0.0;
    CHECK(y1.to_interval() == y2.convert());
    CHECK(y1 == y2);
}

TEST_CASE("affine_form::pow(double) domain checks")
{
    pappus::affine_context ctx;

    CHECK_THROWS_AS(af(ctx, ai(-1.0, 4.0)).pow(0.5), std::invalid_argument);
    CHECK_THROWS_AS(af(ctx, ai(0.0, 4.0)).pow(-0.5), std::invalid_argument);
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

    auto y = x.pow(2.0);
    CHECK(y.length() == 0);
    CHECK(y.to_interval() == ai(4.0, 4.0));

    auto z = af::pow(2.0, x);
    CHECK(z.length() == 0);
    CHECK(z.to_interval() == ai(4.0, 4.0));
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
    CHECK(af(ctx, 1.0).exp().to_interval() == ai(std::exp(1.0)));

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

    CHECK_THROWS_AS(af(ctx, ai(0.0, 1.0)).log(), std::invalid_argument);
    CHECK_THROWS_AS(af(ctx, -1.0).log(),          std::invalid_argument);
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

    CHECK(af(ctx, 0.0).cos().to_interval() == ai(1.0));

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
    CHECK_THROWS_AS(af(ctx, ai(0.0, 2.0)).tan(), std::domain_error);
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

    CHECK(af(ctx, 0.0).cosh().to_interval() == ai(1.0));

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
    CHECK(af(ctx, 1.0).asin().to_interval() == ai(pappus::fp::half_pi_v<double>));

    ai u(-0.8, 0.8);
    af x(ctx, u);
    auto y = x.asin();
    CHECK(sound(y, [](double x) { return std::asin(x); }, u));

    CHECK_THROWS_AS(af(ctx, ai(-2.0, 0.0)).asin(), std::domain_error);
}

TEST_CASE("affine_form::acos()")
{
    pappus::affine_context ctx;

    CHECK(af(ctx, 1.0).acos().to_interval() == ai(0.0));
    CHECK(af(ctx, 0.0).acos().to_interval() == ai(pappus::fp::half_pi_v<double>));

    ai u(-0.8, 0.8);
    af x(ctx, u);
    auto y = x.acos();
    CHECK(sound(y, [](double x) { return std::acos(x); }, u));

    CHECK_THROWS_AS(af(ctx, ai(0.0, 2.0)).acos(), std::domain_error);
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
