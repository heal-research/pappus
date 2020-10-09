#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "aa.h"
#include "pappus.hpp"

#include <Eigen/Eigen>

using af = pappus::affine_form;
using ai = pappus::interval;

constexpr double eps = std::numeric_limits<double>::epsilon();

bool operator==(ai const& lhs, AAInterval const& rhs)
{
    auto l = std::fabs(lhs.inf() - rhs.getlo());
    auto u = std::fabs(lhs.sup() - rhs.gethi());
    return l < eps && u < eps;
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
        CHECK(doctest::Approx(w1.inf()) == w2.getlo());
        CHECK(doctest::Approx(w1.sup()) == w2.gethi());
    };

    SUBCASE("[1, 4] ^ [0, 0.5]") { check_exp(ai(1.0, 4.0), ai(0.0, 0.5)); }
    SUBCASE("[0.1, 2] ^ [-1, 0]") { check_exp(ai(0.1, 2.0), ai(-1.0, 0.0)); }
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
    pappus::interval u1(1, 4);
    pappus::affine_context ctx;
    af x1(ctx, u1);
    auto y1 = x1.sqrt();

    AAInterval u2(u1.inf(), u1.sup());
    AAF x2(u2);
    auto y2 = sqrt(x2);

    CHECK_EQ(y1.to_interval(), y2.convert());
}

TEST_CASE("affine_form::isqrt()")
{
    pappus::interval u1(2, 7);
    pappus::affine_context ctx;
    af x1(ctx, u1);
    auto y1 = x1.isqrt();

    AAInterval u2(u1.inf(), u1.sup());
    AAF x2(u2);
    auto y2 = isqrt(x2);

    CHECK_EQ(y1.to_interval(), y2.convert());
}

/******************************************************
 * Arbitrary expressions tests                        *
 *****************************************************/
TEST_CASE("X^2 + X")
{
    SUBCASE("Dependency problem")
    {
        ai u1(-1, 1);
        pappus::affine_context ctx;

        af x1(ctx, u1);

        auto z1 = x1 * x1 + x1;

        AAInterval u2(u1.inf(), u1.sup());
        AAF x2(u2);
        AAF z2 = x2 * x2 + x2;

        CHECK(z1.to_interval() == z2.convert());

        std::cout << z1.to_interval() << "\n";
        std::cout << z2.convert() << "\n";
    }

    SUBCASE("Rewrite")
    {
        pappus::affine_context ctx;
        ai u1(-1, 1);
        af x1(ctx, u1);
        std::cout << "x1:\n"
                  << x1 << "\n";
        std::cout << "x1 interval: " << x1.to_interval() << "\n";
        auto z1 = (x1 + 0.5).pow(2.0) - 0.25;
        std::cout << "z1:\n"
                  << z1 << "\n";
        std::cout << "z1 interval: " << z1.to_interval() << "\n";

        AAInterval u2(u1.inf(), u1.sup());
        AAF x2(u2);
        AAF z2 = ((x2 + 0.5) ^ 2) - 0.25;
        std::cout << "x2:\n"
                  << x2;
        std::cout << "x2 interval:\n"
                  << x2.convert() << "\n";

        std::cout << "z2:\n"
                  << z2;
        std::cout << "z2 interval:\n"
                  << z2.convert() << "\n";
    }
}
