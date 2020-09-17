#include "interval.hpp"
#include <doctest/doctest.h>

namespace dt = doctest;

const auto inf = pappus::fp::inf;
const auto max_ = std::numeric_limits<double>::max();
const auto min_ = std::numeric_limits<double>::min();

using I = pappus::interval;
const auto F = I::infinite();
const auto Z = I::zero();
const auto E = I::empty();

TEST_CASE("extended real arithmetic rules" * dt::test_suite("IA"))
{
    const auto x = 1.2345;

    pappus::fp::op_add add {};
    pappus::fp::op_sub sub {};
    pappus::fp::op_div div {};
    pappus::fp::op_mul mul {};

    CHECK_EQ(add(inf, x), inf);
    CHECK_EQ(add(-inf, x), -inf);
    CHECK_EQ(add(-inf, -inf), -inf);
    CHECK_EQ(add(inf, inf), inf);
    CHECK(std::isnan(sub(inf, inf)));
    CHECK(std::isnan(add(-inf, inf)));

    CHECK_EQ(mul(-inf, inf), -inf);
    CHECK_EQ(mul(-inf, -inf), inf);
    CHECK_EQ(mul(0, inf), 0);
    CHECK_EQ(mul(inf, 0), 0);
    CHECK_EQ(mul(0, -inf), 0);
    CHECK_EQ(mul(-inf, 0), 0);

    CHECK_EQ(div(x, inf), 0);
    CHECK_EQ(div(x, -inf), 0);
    CHECK_EQ(div(x, 0), inf);
    CHECK_EQ(div(x, -0), inf);
}

TEST_CASE("addition" * dt::test_suite("IA"))
{
    CHECK_EQ(F + F, F);
    CHECK_EQ(E + E, E);
    CHECK_EQ(+I(3,4), I(3,4));
    CHECK_EQ(I(3,4) + I(-6,2), I(-3,6));
    CHECK_EQ(I(3,4) + 5, I(8,9));
    CHECK_EQ(5 + I(3,4), I(8,9));
}

TEST_CASE("subtraction" * dt::test_suite("IA"))
{
    CHECK_EQ(-I(4,5), I(-5,-4));
    CHECK_EQ(-I(-max_,inf), I(-inf,max_));
    CHECK(E.is_empty());
    CHECK((-E).is_empty());
    CHECK_EQ(F - F,F);
    CHECK_EQ(I(4,5) - I(3,7), I(-3,2));
    CHECK_EQ(I(-inf,-max_) - I(-inf,-max_),F);
    CHECK_EQ(I(3,4) - 5, I(-2,-1));
}

TEST_CASE("multiplication" * dt::test_suite("IA"))
{
    CHECK_EQ(E * E, E);
    CHECK_EQ(F * F, F);
    CHECK_EQ(E * I(4,5), E);
    CHECK_EQ(I(4,5) * E, E);
    CHECK_EQ(E * Z, E);
    CHECK_EQ(Z * E, E);
    CHECK_EQ(3 * E, E);
    CHECK_EQ(E * 3, E);
    CHECK_EQ(E * 0.0, E);
    CHECK_EQ(0.0 * E, E);
    CHECK_EQ(F * 0,Z);
    CHECK_EQ(I(2,4) * 5, I(10,20));
    CHECK_EQ(I(-inf,max_) * 0,Z);
    CHECK_EQ(I(max_,inf) * Z,Z);
    CHECK_EQ(Z * I(max_,inf),Z);
    CHECK_EQ(Z * I(-inf,max_),Z);

    CHECK_EQ(I(4,5) * I(2,3), I(8,15)); // P*P
    CHECK_EQ(I(4,5) * I(-2,3), I(-10,15)); // P*M
    CHECK_EQ(I(4,5) * I(-3,-2), I(-15,-8)); // P*N
    CHECK_EQ(I(-3,-2) * I(4,5), I(-15,-8)); // N*P
    CHECK_EQ(I(-3,2) * I(-4,5), I(-15,12)); // M*M
    CHECK_EQ(I(-3,2) * I(-5,-4), I(-10,15)); // M*N
    CHECK_EQ(I(-3,-2) * I(4,5), I(-15,-8)); // N*P
    CHECK_EQ(I(-3,-2) * I(-4,5), I(-15,12)); // N*M
    CHECK_EQ(I(-3,-2) * I(-5,-4), I(8,15)); // N*N
}

TEST_CASE("division" * dt::test_suite("IA"))
{
    CHECK_EQ(E / E, E);
    CHECK_EQ(E / F, E);
    CHECK_EQ(F / E, E);
    CHECK_EQ(I(3,4) / 2, I(1.5,2));
    CHECK_EQ(I(3,4) / -2, I(-2,-1.5));
    CHECK_EQ(I(-4,-3) / 2, I(-2,-1.5));
    CHECK_EQ(I(-4,-3) / -2, I(1.5,2));
    CHECK_EQ(I(3,4) / 0, E);
    CHECK_EQ(F / 0, E);

    CHECK_EQ(I(-3,-2) / I(-5,-4), I("0.4","0.75")); // N1/N1
    CHECK_EQ(I(-3,-2) / Z, E); // N1/Z
    CHECK_EQ(I(-3,-2) / I(-4,0), I(0.5,inf)); // N1/N0
    CHECK_EQ(I(-3,-2) / I(-4, 5), F); // N1/M
    CHECK_EQ(I(-3,-2) / I(0, 4), I(-inf, -0.5)); // N1/P0
    CHECK_EQ(I(-3,-2) / I(4, 5), I("-0.75", "-0.4")); // N1/P1
    CHECK_EQ(Z / Z, E); // Z/Z
    CHECK_EQ(Z / I(-3,-2), Z); // Z/N
    CHECK_EQ(Z / I(-3,2), Z); // Z/M
    CHECK_EQ(Z / I(2,3), Z); // Z/P
    CHECK_EQ(I(-2,0) / I(-5,-4), I(0, 0.5)); // N0/N1
    CHECK_EQ(I(-2,0) / Z, E); // N0/Z
    CHECK_EQ(I(-2,0) / I(-4,0), I(0, inf)); // N0/N0
    CHECK_EQ(I(-3,0) / I(-4,2), I::infinite()); // N0/M
    CHECK_EQ(I(-2,0) / I(0,4), I(-inf, 0)); // N0/P0
    CHECK_EQ(I(-2,0) / I(4,5), I(-0.5, 0)); // N0/P1
    CHECK_EQ(I(-3,2) / I(-5,-4), I(-0.5, 0.75)); // M/N1
    CHECK_EQ(I(-3,2) / Z, E); // M/Z
    CHECK_EQ(I(-3,2) / I(-4,0), F); // M/N0
    CHECK_EQ(I(-3,2) / I(-4,5), F); // M/M
    CHECK_EQ(I(-3,2) / I(0,5), F); // M/P0
    CHECK_EQ(I(-3,2) / I(4,5), I(-0.75,0.5)); // M/P1
    CHECK_EQ(I(0,2) / I(-5,-4), I(-0.5,0.0)); // P0/N1
    CHECK_EQ(I(0,2) / Z, E); // P0/Z
    CHECK_EQ(I(0,2) / I(-4,0), I(-inf,0)); // P0/N0
    CHECK_EQ(I(0,2)/I(-4,5), F); // P0/M
    CHECK_EQ(I(0,2) / I(0,5), I(0,inf)); // P0/P0
    CHECK_EQ(I(0,2) / I(4,5), I(0,0.5)); // P0/P1
    CHECK_EQ(I(2,3) / I(-5,-4), I("-0.75","-0.4")); // P1/N1
    CHECK_EQ(I(2,3) / Z, E); // P1/Z
    CHECK_EQ(I(2,3) / I(-4,0), I(-inf, -0.5)); // P1/N0
    CHECK_EQ(I(2,3)/I(-4,5), F); // P1/M
    CHECK_EQ(I(2,3) / I(0,4), I(0.5, inf)); // P1/P0
    CHECK_EQ(I(2,3) / I(4,5), I("0.4","0.75")); // P1/P1
}
