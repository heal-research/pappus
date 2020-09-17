#include <doctest/doctest.h>

#include "interval.hpp"

using I = pappus::interval;

TEST_CASE("extended real arithmetic rules")
{
    const auto inf = pappus::fp::inf;
    const auto x = 1.2345;

    pappus::fp::op_add add{};
    pappus::fp::op_sub sub{};
    pappus::fp::op_div div{};
    pappus::fp::op_mul mul{};

    CHECK(add(inf, x) == inf);
    CHECK(add(-inf, x) == -inf);
    CHECK(add(-inf, -inf) == -inf);
    CHECK(add(inf, inf) == inf);
    CHECK(std::isnan(sub(inf, inf)));
    CHECK(std::isnan(add(-inf, inf)));

    CHECK(mul(-inf, inf) == -inf);
    CHECK(mul(-inf, -inf) == inf);
    CHECK(mul(0, inf) == 0);
    CHECK(mul(inf, 0) == 0);
    CHECK(mul(0, -inf) == 0);
    CHECK(mul(-inf, 0) == 0);

    CHECK(div(x, inf) == 0);
    CHECK(div(x, -inf) == 0);
    CHECK(div(x, 0) == inf);
    CHECK(div(x, -0) == inf);
}

TEST_CASE("interval::operator+")
{
    SUBCASE("[a1, b1] + [a2, b2]")
    {
        auto a = I(2, 3);
        auto b = I(2, 5);
        auto c = a + b;
        CHECK(c == I(4, 8));
    }

    SUBCASE("[a1, b1] + [inf, b2]")
    {
        auto a = I(2, 3);
        auto b = I(-INFINITY, 5);
        auto c = a + b;
        std::cout << c << "\n";
        CHECK(c == I(-INFINITY, 8));
    }

    SUBCASE("[inf, b1] + [a2, b2]")
    {
        auto a = I(-INFINITY, 3);
        auto b = I(2, 5);
        auto c = a + b;
        CHECK(c == I(-INFINITY, 8));
    }

    SUBCASE("[inf, b1] + [a2, inf]")
    {
        auto a = I(-INFINITY, 3);
        auto b = I(2, +INFINITY);
        auto c = a + b;
        CHECK(c == I(-INFINITY, INFINITY));
    }

    SUBCASE("[inf, inf] + [inf, inf]")
    {
        auto a = I(-INFINITY, INFINITY);
        auto b = I(-INFINITY, INFINITY);
        auto c = a + b;
        CHECK(c == I(-INFINITY, INFINITY));
    }
}
