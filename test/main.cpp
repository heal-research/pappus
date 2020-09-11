#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "aa.h"
#include "pappus.hpp"

#include <Eigen/Eigen>

using af = pappus::affine_form;
using ai = pappus::affine_interval;

constexpr double eps = std::numeric_limits<double>::epsilon();

bool operator==(ai const& lhs, AAInterval const& rhs)
{
    auto l = std::fabs(lhs.lower() - rhs.getlo());
    auto u = std::fabs(lhs.upper() - rhs.gethi());
    return l < eps && u < eps; 
}

bool operator==(AAInterval const& lhs, ai const& rhs)
{
    return rhs == lhs;
}

TEST_CASE("operator+(af)")
{
    // pappus
    ai u1(-2, 3);
    ai v1(-1, 1);

    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1(ctx, v1);

    af z1 = x1 + y1;

    // aaflib
    AAInterval u2(u1.lower(), u1.upper());
    AAInterval v2(v1.lower(), v1.upper());

    AAF x2(u2);
    AAF y2(v2);
    AAF z2 = x2 + y2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("operator+(double, af)")
{
    ai u1(-2, 3);
    double v = 2;

    pappus::affine_context ctx;
    af x1(ctx, u1);

    af z1 = v + x1;

    AAInterval u2(u1.lower(), u1.upper());
    AAF x2(u2);
    AAF z2 = v + x2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("operator+=(af)")
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

TEST_CASE("operator-(af)")
{
    // pappus
    ai u1(-2, 3);
    ai v1(-1, 1);

    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1(ctx, v1);

    af z1 = x1 - y1;

    // aaflib
    AAInterval u2(u1.lower(), u1.upper());
    AAInterval v2(v1.lower(), v1.upper());

    AAF x2(u2);
    AAF y2(v2);
    AAF z2 = x2 - y2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("operator-(double, af)")
{
    ai u1(-2, 3);
    double v = 2;

    pappus::affine_context ctx;
    af x1(ctx, u1);

    af z1 = v - x1;

    AAInterval u2(u1.lower(), u1.upper());
    AAF x2(u2);
    AAF z2 = v - x2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("operator-=(af)")
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

TEST_CASE("operator*(af)")
{
    ai u1(-2, 3);
    ai v1(-1, 1);

    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1(ctx, v1);

    af z1 = x1 * y1;

    AAInterval u2(u1.lower(), u1.upper());
    AAInterval v2(v1.lower(), v1.upper());

    AAF x2(u2);
    AAF y2(v2);
    AAF z2 = x2 * y2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("operator*(double, af)")
{
    ai u1(-2, 3);
    double v = 2;

    pappus::affine_context ctx;
    af x1(ctx, u1);

    af z1 = v * x1;

    AAInterval u2(u1.lower(), u1.upper());
    AAF x2(u2);
    AAF z2 = v * x2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("operator*=(af)")
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

TEST_CASE("operator/(af)")
{
    ai u1(2, 3);
    ai v1(1.5, 2);

    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1(ctx, v1);

    af z1 = x1 / y1;

    AAInterval u2(u1.lower(), u1.upper());
    AAInterval v2(v1.lower(), v1.upper());

    AAF x2(u2);
    AAF y2(v2);
    AAF z2 = x2 / y2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("operator/(double, af)")
{
    ai u1(2, 3);
    double v = 2;

    pappus::affine_context ctx;
    af x1(ctx, u1);

    af z1 = v / x1;

    AAInterval u2(u1.lower(), u1.upper());
    AAF x2(u2);
    AAF z2 = v / x2;

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("affine_form.inv()")
{
    ai u1(1, 4);
    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1 = x1.inv();

    AAInterval u2(u1.lower(), u1.upper());
    AAF x2(u2);
    AAF y2 = inv(x2);

    CHECK(y1.to_interval() == y2.convert());
}
