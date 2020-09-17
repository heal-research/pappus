#include <doctest/doctest.h>

#include "interval.hpp"

using I = pappus::interval;

TEST_CASE("extended real arithmetic rules")
{
    const auto inf = std::numeric_limits<double>::infinity();
    const auto x = 1.2345;

    CHECK(inf + x == inf);
    CHECK(-inf + x == -inf);
    CHECK(-inf + (-inf) == -inf);
    CHECK(inf + inf == inf);
    CHECK(inf - inf == 0);
    CHECK(-inf + inf == 0);

    CHECK(-inf * inf == -inf);
    CHECK(-inf * -inf == inf);
    CHECK(0 * inf == 0);
    CHECK(inf * 0 == 0);
    CHECK(0 * -inf == 0);
    CHECK(-inf * 0 == 0);

    CHECK(x / inf == 0);
    CHECK(x / -inf == 0);
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
