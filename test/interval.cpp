#include "interval.hpp"
#include <iomanip>
#include <doctest/doctest.h>

namespace dt = doctest;

namespace fp = pappus::fp;

const auto inf = fp::inf;
const auto max = std::numeric_limits<double>::max();
const auto min = std::numeric_limits<double>::min();

using I = pappus::interval;
const auto F = I::infinite();
const auto Z = I::zero();
const auto E = I::empty();

TEST_CASE("addition" * dt::test_suite("IA"))
{
    // operator+
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
    CHECK_EQ(-I(-max,inf), I(-inf,max));
    CHECK(E.is_empty());
    CHECK((-E).is_empty());
    CHECK_EQ(F - F,F);
    CHECK_EQ(I(4,5) - I(3,7), I(-3,2));
    CHECK_EQ(I(-inf,-max) - I(-inf,-max),F);
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
    CHECK_EQ(5 * I(2,4), I(10,20));
    CHECK_EQ(I(-inf,max) * 0,Z);
    CHECK_EQ(I(max,inf) * Z,Z);
    CHECK_EQ(Z * I(max,inf),Z);
    CHECK_EQ(Z * I(-inf,max),Z);

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
    CHECK_EQ(I(-3,0) / I(-4,2), F); // N0/M
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
    CHECK_EQ(I(0,2) / I(-4,5), F); // P0/M
    CHECK_EQ(I(0,2) / I(0,5), I(0,inf)); // P0/P0
    CHECK_EQ(I(0,2) / I(4,5), I(0,0.5)); // P0/P1
    CHECK_EQ(I(2,3) / I(-5,-4), I("-0.75","-0.4")); // P1/N1
    CHECK_EQ(I(2,3) / Z, E); // P1/Z
    CHECK_EQ(I(2,3) / I(-4,0), I(-inf, -0.5)); // P1/N0
    CHECK_EQ(I(2,3) / I(-4,5), F); // P1/M
    CHECK_EQ(I(2,3) / I(0,4), I(0.5, inf)); // P1/P0
    CHECK_EQ(I(2,3) / I(4,5), I("0.4","0.75")); // P1/P1
    CHECK_EQ(1 / I(10, 10), I("0.1"));
    CHECK_EQ(1 / I(10, 10), I(10, 10).inv());
}

TEST_CASE("inverse" * dt::test_suite("IA"))
{
    CHECK_EQ(E.inv(), E);
    CHECK_EQ(Z.inv(), E);
    CHECK_EQ(F.inv(), F);
    CHECK_EQ(I(10,10).inv(), I("0.1"));
    CHECK_EQ(I(0.5, 1.5).inv(), 1 / I(0.5, 1.5));
    CHECK_EQ(I(-4,-2).inv(), I(-0.5, -0.25));
    CHECK_EQ(I(-2,0).inv(), I(-inf, -0.5));
    CHECK_EQ(I(-4,2).inv(), F);
    CHECK_EQ(I(-2,4).inv(), F);
    CHECK_EQ(I(0,2).inv(), I(0.5, inf));
    CHECK_EQ(I(2,4).inv(), I(0.25, 0.5));
}

TEST_CASE("trigonometric functions" * dt::test_suite("IA"))
{
    auto print = [](auto v) { std::cout << std::setprecision(60) << std::fixed << v << "\n"; };

    SUBCASE("sin")
    {
        CHECK_EQ(I(0, fp::pi).sin(), I(0, 1));
        CHECK_EQ(F.sin(), I(-1, 1));
    }

    SUBCASE("cos")
    {
        CHECK_EQ(E.cos(), E);
        CHECK_EQ(I(-1.0,-0.5).cos(), I("0.540302305868139765010482733487151563167572021484375", 
                                       "0.8775825618903727587394314468838274478912353515625"));
        CHECK_EQ(I(-10,10).cos(), I(-1,1));
        CHECK_EQ(F.cos(), I(-1, 1));
        CHECK_EQ(I(fp::pi, fp::pi).cos(), I(-1, -1));
        CHECK_EQ(I(3,3.5).cos(), I("-1","-0.936456687290796341294196736271260306239128112792969")); 
        CHECK_EQ(I(-3.5,-3).cos(), I("-1","-0.936456687290796341294196736271260306239128112792969"));
        CHECK_EQ(I(-3.5,3).cos(), I(-1,1));

        CHECK_EQ(I(10,12).cos(), I("-0.83907152907645243811174395887064747512340545654296875",
                                   "0.84385395873249213760658449245966039597988128662109375"));

        CHECK_EQ(I(13,14).cos(), I("0.1367372182078336051436195930364192463457584381103515625",
                                   "0.90744678145019619375233332903007976710796356201171875"));

        CHECK_EQ(I(10,14).cos(), I("-0.83907152907645243811174395887064747512340545654296875","1"));
        CHECK_EQ(I(14,16).cos(), I("-1","0.1367372182078336051436195930364192463457584381103515625"));
        print(I(14,16).cos().upper());
        //TEST_EQ(cos(interval(14,16)),interval("[-1,0.136737218207833]"));
        //TEST_EQ(cos(interval(-11,-10)),interval("[-0.839071529076452,0.004425697988051]"));
        //TEST_EQ(cos(interval(-14,-13)),interval("[0.136737218207833,0.907446781450197]"));
        //TEST_EQ(cos(interval(-16,-14)),interval("[-1,0.136737218207833]"));
        //TEST_EQ(cos(interval(-102,-100)),interval("[0.101585703696621,1]"));
        //TEST_EQ(cos(interval(4.6e15,4.7e15)),interval(-1,1));
        //TEST_EQ(cos(interval(4503599627370495,4503599627370496)),interval("[-0.48553486774222065, 0.4732928859543091]"));
    }
}
