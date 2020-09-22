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

TEST_CASE("interval")
{
    CHECK_EQ(I(0, fp::pi).diameter(), fp::pi);
}

TEST_CASE("addition" * dt::test_suite("IA"))
{
    // operator+
    CHECK_EQ(F + F, F);
    CHECK_EQ(E + E, E);
    CHECK_EQ(+I(3, 4), I(3, 4));
    CHECK_EQ(I(3, 4) + I(-6, 2), I(-3, 6));
    CHECK_EQ(I(3, 4) + 5, I(8, 9));
    CHECK_EQ(5 + I(3, 4), I(8, 9));
}

TEST_CASE("subtraction" * dt::test_suite("IA"))
{
    CHECK_EQ(-I(4, 5), I(-5, -4));
    CHECK_EQ(-I(-max, inf), I(-inf, max));
    CHECK(E.is_empty());
    CHECK((-E).is_empty());
    CHECK_EQ(F - F, F);
    CHECK_EQ(I(4, 5) - I(3, 7), I(-3, 2));
    CHECK_EQ(I(-inf, -max) - I(-inf, -max), F);
    CHECK_EQ(I(3, 4) - 5, I(-2, -1));
}

TEST_CASE("multiplication" * dt::test_suite("IA"))
{
    CHECK_EQ(E * E, E);
    CHECK_EQ(F * F, F);
    CHECK_EQ(E * I(4, 5), E);
    CHECK_EQ(I(4, 5) * E, E);
    CHECK_EQ(E * Z, E);
    CHECK_EQ(Z * E, E);
    CHECK_EQ(3 * E, E);
    CHECK_EQ(E * 3, E);
    CHECK_EQ(E * 0.0, E);
    CHECK_EQ(0.0 * E, E);
    CHECK_EQ(F * 0, Z);
    CHECK_EQ(I(2, 4) * 5, I(10, 20));
    CHECK_EQ(5 * I(2, 4), I(10, 20));
    CHECK_EQ(I(-inf, max) * 0, Z);
    CHECK_EQ(I(max, inf) * Z, Z);
    CHECK_EQ(Z * I(max, inf), Z);
    CHECK_EQ(Z * I(-inf, max), Z);

    CHECK_EQ(I(4, 5) * I(2, 3), I(8, 15)); // P*P
    CHECK_EQ(I(4, 5) * I(-2, 3), I(-10, 15)); // P*M
    CHECK_EQ(I(4, 5) * I(-3, -2), I(-15, -8)); // P*N
    CHECK_EQ(I(-3, -2) * I(4, 5), I(-15, -8)); // N*P
    CHECK_EQ(I(-3, 2) * I(-4, 5), I(-15, 12)); // M*M
    CHECK_EQ(I(-3, 2) * I(-5, -4), I(-10, 15)); // M*N
    CHECK_EQ(I(-3, -2) * I(4, 5), I(-15, -8)); // N*P
    CHECK_EQ(I(-3, -2) * I(-4, 5), I(-15, 12)); // N*M
    CHECK_EQ(I(-3, -2) * I(-5, -4), I(8, 15)); // N*N
}

TEST_CASE("division" * dt::test_suite("IA"))
{
    CHECK_EQ(E / E, E);
    CHECK_EQ(E / F, E);
    CHECK_EQ(F / E, E);
    CHECK_EQ(I(3, 4) / 2, I(1.5, 2));
    CHECK_EQ(I(3, 4) / -2, I(-2, -1.5));
    CHECK_EQ(I(-4, -3) / 2, I(-2, -1.5));
    CHECK_EQ(I(-4, -3) / -2, I(1.5, 2));
    CHECK_EQ(I(3, 4) / 0, E);
    CHECK_EQ(F / 0, E);

    CHECK_EQ(I(-3, -2) / I(-5, -4), I("0.4", "0.75")); // N1/N1
    CHECK_EQ(I(-3, -2) / Z, E); // N1/Z
    CHECK_EQ(I(-3, -2) / I(-4, 0), I(0.5, inf)); // N1/N0
    CHECK_EQ(I(-3, -2) / I(-4, 5), F); // N1/M
    CHECK_EQ(I(-3, -2) / I(0, 4), I(-inf, -0.5)); // N1/P0
    CHECK_EQ(I(-3, -2) / I(4, 5), I("-0.75", "-0.4")); // N1/P1
    CHECK_EQ(Z / Z, E); // Z/Z
    CHECK_EQ(Z / I(-3, -2), Z); // Z/N
    CHECK_EQ(Z / I(-3, 2), Z); // Z/M
    CHECK_EQ(Z / I(2, 3), Z); // Z/P
    CHECK_EQ(I(-2, 0) / I(-5, -4), I(0, 0.5)); // N0/N1
    CHECK_EQ(I(-2, 0) / Z, E); // N0/Z
    CHECK_EQ(I(-2, 0) / I(-4, 0), I(0, inf)); // N0/N0
    CHECK_EQ(I(-3, 0) / I(-4, 2), F); // N0/M
    CHECK_EQ(I(-2, 0) / I(0, 4), I(-inf, 0)); // N0/P0
    CHECK_EQ(I(-2, 0) / I(4, 5), I(-0.5, 0)); // N0/P1
    CHECK_EQ(I(-3, 2) / I(-5, -4), I(-0.5, 0.75)); // M/N1
    CHECK_EQ(I(-3, 2) / Z, E); // M/Z
    CHECK_EQ(I(-3, 2) / I(-4, 0), F); // M/N0
    CHECK_EQ(I(-3, 2) / I(-4, 5), F); // M/M
    CHECK_EQ(I(-3, 2) / I(0, 5), F); // M/P0
    CHECK_EQ(I(-3, 2) / I(4, 5), I(-0.75, 0.5)); // M/P1
    CHECK_EQ(I(0, 2) / I(-5, -4), I(-0.5, 0.0)); // P0/N1
    CHECK_EQ(I(0, 2) / Z, E); // P0/Z
    CHECK_EQ(I(0, 2) / I(-4, 0), I(-inf, 0)); // P0/N0
    CHECK_EQ(I(0, 2) / I(-4, 5), F); // P0/M
    CHECK_EQ(I(0, 2) / I(0, 5), I(0, inf)); // P0/P0
    CHECK_EQ(I(0, 2) / I(4, 5), I(0, 0.5)); // P0/P1
    CHECK_EQ(I(2, 3) / I(-5, -4), I("-0.75", "-0.4")); // P1/N1
    CHECK_EQ(I(2, 3) / Z, E); // P1/Z
    CHECK_EQ(I(2, 3) / I(-4, 0), I(-inf, -0.5)); // P1/N0
    CHECK_EQ(I(2, 3) / I(-4, 5), F); // P1/M
    CHECK_EQ(I(2, 3) / I(0, 4), I(0.5, inf)); // P1/P0
    CHECK_EQ(I(2, 3) / I(4, 5), I("0.4", "0.75")); // P1/P1
    CHECK_EQ(1 / I(10, 10), I("0.1"));
    CHECK_EQ(1 / I(10, 10), I(10, 10).inv());
}

TEST_CASE("inverse" * dt::test_suite("IA"))
{
    CHECK_EQ(E.inv(), E);
    CHECK_EQ(Z.inv(), E);
    CHECK_EQ(F.inv(), F);
    CHECK_EQ(I(10, 10).inv(), I("0.1"));
    CHECK_EQ(I(0.5, 1.5).inv(), 1 / I(0.5, 1.5));
    CHECK_EQ(I(-4, -2).inv(), I(-0.5, -0.25));
    CHECK_EQ(I(-2, 0).inv(), I(-inf, -0.5));
    CHECK_EQ(I(-4, 2).inv(), F);
    CHECK_EQ(I(-2, 4).inv(), F);
    CHECK_EQ(I(0, 2).inv(), I(0.5, inf));
    CHECK_EQ(I(2, 4).inv(), I(0.25, 0.5));
}

TEST_CASE("trigonometric functions" * dt::test_suite("IA"))
{
    SUBCASE("sin")
    {
        CHECK_EQ(F.sin(), I(-1, 1));
        CHECK_EQ(E.sin(), E);
        CHECK_EQ(Z.sin(), Z);
        CHECK_EQ(I(0, fp::pi).sin(), I(0, 1));
        CHECK_EQ(I(fp::pi, fp::pi).sin(), I(1.2246467991473529607e-16, 1.2246467991473532072e-16));
        CHECK_EQ(I(0.5, 0.5).sin(), I(0.47942553860420295, 0.47942553860420301));
        CHECK_EQ(I(0.5, 1.67).sin(), I(0.47942553860420295, 1.0));
        CHECK_EQ(I(0.5, 8.5).sin(), I(-1, 1));
        CHECK_EQ(I(1.67, 3.2).sin(), I(-5.8374143427580093e-02, 9.9508334981018021e-01));
        CHECK_EQ(I(1.67,3.2).sin(), I(-5.8374143427580093e-02, 9.9508334981018021e-01));
        CHECK_EQ(I(2.1, 5.6).sin(), I(-1.0, 0.8632093666488738));
        CHECK_EQ(I(-4.5, 0.1).sin(), I(-1.0, 0.9775301176650971));
        CHECK_EQ(I(1.3, 6.3).sin(), I(-1.0, 1.0));
        CHECK_EQ(I(-1.0, -0.5).sin(), I(-0.8414709848078966, -0.47942553860420295));
        CHECK_EQ(I(-10, 10).sin(), I(-1, 1));
        CHECK_EQ(I(3, 3.5).sin(), I(-0.3507832276896199, 0.14112000805986724)); 
        CHECK_EQ(I(-3.5, -3).sin(), I(-0.14112000805986724, 0.3507832276896199));
        CHECK_EQ(I(-3.5, 3).sin(), I(-1, 1));
        CHECK_EQ(I(10, 12).sin(), I(-1, -0.5365729180004349));
        CHECK_EQ(I(13, 14).sin(), I(0.4201670368266409, 0.9906073556948704));
        CHECK_EQ(I(10, 14).sin(), I(-1, 0.9906073556948704));
        CHECK_EQ(I(14, 16).sin(), I(-0.2879033166650653, 1));
        CHECK_EQ(I(-11, -10).sin(), I(0.5440211108893698, 1));
        CHECK_EQ(I(-14, -13).sin(), I(-0.9906073556948704, -0.4201670368266409));
        CHECK_EQ(I(-16, -14).sin(), I(-1, 0.2879033166650653));
        CHECK_EQ(I(-102, -100).sin(), I(-0.9948267913584065, 0.5063656411097589));
        CHECK_EQ(I(4.6e15, 4.7e15).sin(), I(-1, 1));
        CHECK_EQ(I(4503599627370495, 4503599627370496).sin(), I(0.8742173026236351, 1.0));
    }

    SUBCASE("cos")
    {
        CHECK_EQ(F.cos(), I(-1, 1));
        CHECK_EQ(E.cos(), E);
        CHECK_EQ(Z.cos(), I(1, 1));
        CHECK_EQ(I(0.5).cos(), I(0.87758256189037265, 0.87758256189037276));
        CHECK_EQ(I(0.5, 1.67).cos(), I(-0.09904103659872802, 0.8775825618903728));
        CHECK_EQ(I(2.1, 5.6).cos(), I(-1.0, 0.7755658785102496));
        CHECK_EQ(I(0.5, 8.5).cos(), I(-1.0, 1.0));
        CHECK_EQ(I(1.67, 3.2).cos(), I(-1.0, -0.09904103659872801));
        CHECK_EQ(I(-1.0, -0.5).cos(), I(0.5403023058681397, 0.8775825618903728));
        CHECK_EQ(I(-10, 10).cos(), I(-1, 1));
        CHECK_EQ(I(fp::pi, fp::pi).cos(), I(-1, -0.9999999999999999));
        CHECK_EQ(I(3, 3.5).cos(), I(-1, -0.9364566872907962)); 
        CHECK_EQ(I(-3.5, -3).cos(), I(-1, -0.9364566872907962));
        CHECK_EQ(I(-3.5, 3).cos(), I(-1, 1));
        CHECK_EQ(I(10, 12).cos(), I(-0.8390715290764525, 0.8438539587324921));
        CHECK_EQ(I(13, 14).cos(), I(0.13673721820783358, 0.9074467814501963));
        CHECK_EQ(I(10, 14).cos(), I(-0.8390715290764525, 1));
        CHECK_EQ(I(14, 16).cos(), I(-1, 0.1367372182078336));
        CHECK_EQ(I(-11, -10).cos(), I(-0.8390715290764525, 0.004425697988050786));
        CHECK_EQ(I(-14, -13).cos(), I(0.13673721820783358, 0.9074467814501963));
        CHECK_EQ(I(-16, -14).cos(), I(-1, 0.1367372182078336));
        CHECK_EQ(I(-102, -100).cos(), I(0.10158570369662133, 1));
        CHECK_EQ(I(4.6e15, 4.7e15).cos(), I(-1, 1));
        CHECK_EQ(I(4503599627370495, 4503599627370496).cos(), I(-0.48553486774222065, 0.4732928859543091));
    }

    SUBCASE("tan")
    {
        CHECK_EQ(I(0.5).tan(), I(0.5463024898437905, 0.5463024898437906));
        CHECK_EQ(I(0.5, 1.67).tan(), F);
        CHECK_EQ(I(1.67, 3.2).tan(), I(-10.047182299210307, 0.05847385445957865));
        CHECK_EQ(I(6.638314112824137, 8.38263151220128).tan(), F);
    }
}

