#include "pappus/interval/interval.hpp"
#include "pappus/interval/box.hpp"
#include "pappus/interval/iterator.hpp"
#include <iomanip>
#include <catch2/catch_all.hpp>
#include <functional>
#include <random>
#include <vector>

namespace fp = pappus::fp;

const auto inf = fp::inf;
const auto max = std::numeric_limits<double>::max();
const auto min = std::numeric_limits<double>::min();

using I = pappus::interval<double>;
const auto F = I::infinite();
const auto Z = I::zero();
const auto E = I::empty();

namespace {

bool interval_approx_equal(I const& lhs, I const& rhs, double scale = 32.0)
{
    if (lhs.is_empty() || rhs.is_empty()) return lhs == rhs;

    auto close = [scale](double a, double b) {
        if (a == b) return true;
        if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && std::isnan(b);
        if (std::isinf(a) || std::isinf(b)) return a == b;

        auto tol = scale * std::numeric_limits<double>::epsilon() * std::max({1.0, std::fabs(a), std::fabs(b)});
        tol = std::max(tol, scale * std::numeric_limits<double>::denorm_min());
        return std::fabs(a - b) <= tol;
    };

    return close(lhs.inf(), rhs.inf()) && close(lhs.sup(), rhs.sup());
}

} // namespace

#define CHECK_IA_APPROX(actual, expected) CHECK(interval_approx_equal((actual), (expected)))

TEST_CASE("interval")
{
    CHECK(I(0, fp::pi).diameter() == fp::pi);
    CHECK(I::force_interval(3, 2) == I(2, 3));

    SECTION("subdivision_iterator")
    {
        I iv(0, 1);
        auto subdiv = iv.split(5);
        std::vector<I> parts;
        for (auto vv : subdiv)
            parts.push_back(vv);

        REQUIRE(parts.size() == 5);
        for (size_t i = 0; i < parts.size(); ++i)
            CHECK(parts[i] == iv.segment(i, parts.size()));

        I shifted(-2, 3);
        std::vector<I> shifted_parts;
        for (auto vv : shifted.split(5))
            shifted_parts.push_back(vv);

        REQUIRE(shifted_parts.size() == 5);
        for (size_t i = 0; i < shifted_parts.size(); ++i)
            CHECK(shifted_parts[i] == shifted.segment(i, shifted_parts.size()));
    }
}

TEST_CASE("addition", "[IA]")
{
    // operator+
    CHECK(F + F == F);
    CHECK(E + E == E);
    CHECK(+I(3, 4) == I(3, 4));
    CHECK(I(3, 4) + I(-6, 2) == I(-3, 6));
    CHECK(I(3, 4) + 5 == I(8, 9));
    CHECK(5 + I(3, 4) == I(8, 9));
}

TEST_CASE("subtraction", "[IA]")
{
    CHECK(-I(4, 5) == I(-5, -4));
    CHECK(-I(-max, inf) == I(-inf, max));
    CHECK(E.is_empty());
    CHECK((-E).is_empty());
    CHECK(F - F == F);
    CHECK(I(4, 5) - I(3, 7) == I(-3, 2));
    CHECK(I(-inf, -max) - I(-inf, -max) == F);
    CHECK(I(3, 4) - 5 == I(-2, -1));
}

TEST_CASE("multiplication", "[IA]")
{
    CHECK(E * E == E);
    CHECK(F * F == F);
    CHECK(E * I(4, 5) == E);
    CHECK(I(4, 5) * E == E);
    CHECK(E * Z == E);
    CHECK(Z * E == E);
    CHECK(3 * E == E);
    CHECK(E * 3 == E);
    CHECK(E * 0.0 == E);
    CHECK(0.0 * E == E);
    CHECK(F * 0 == Z);
    CHECK(I(2, 4) * 5 == I(10, 20));
    CHECK(5 * I(2, 4) == I(10, 20));
    CHECK(I(-inf, max) * 0 == Z);
    CHECK(I(max, inf) * Z == Z);
    CHECK(Z * I(max, inf) == Z);
    CHECK(Z * I(-inf, max) == Z);

    CHECK(I(4, 5) * I(2, 3) == I(8, 15)); // P*P
    CHECK(I(4, 5) * I(-2, 3) == I(-10, 15)); // P*M
    CHECK(I(4, 5) * I(-3, -2) == I(-15, -8)); // P*N
    CHECK(I(-3, -2) * I(4, 5) == I(-15, -8)); // N*P
    CHECK(I(-3, 2) * I(-4, 5) == I(-15, 12)); // M*M
    CHECK(I(-3, 2) * I(-5, -4) == I(-10, 15)); // M*N
    CHECK(I(-3, -2) * I(4, 5) == I(-15, -8)); // N*P
    CHECK(I(-3, -2) * I(-4, 5) == I(-15, 12)); // N*M
    CHECK(I(-3, -2) * I(-5, -4) == I(8, 15)); // N*N
}

TEST_CASE("division", "[IA]")
{
    CHECK(E / E == E);
    CHECK(E / F == E);
    CHECK(F / E == E);
    CHECK(I(3, 4) / 2 == I(1.5, 2));
    CHECK(I(3, 4) / -2 == I(-2, -1.5));
    CHECK(I(-4, -3) / 2 == I(-2, -1.5));
    CHECK(I(-4, -3) / -2 == I(1.5, 2));
    CHECK(I(3, 4) / 0 == E);
    CHECK(F / 0 == E);

    CHECK_IA_APPROX(I(-3, -2) / I(-5, -4), I("0.4", "0.75")); // N1/N1
    CHECK(I(-3, -2) / Z == E); // N1/Z
    CHECK(I(-3, -2) / I(-4, 0) == I(0.5, inf)); // N1/N0
    CHECK(I(-3, -2) / I(-4, 5) == F); // N1/M
    CHECK(I(-3, -2) / I(0, 4) == I(-inf, -0.5)); // N1/P0
    CHECK(I(-3, -2) / I(4, 5) == I("-0.75", "-0.4")); // N1/P1
    CHECK(Z / Z == E); // Z/Z
    CHECK(Z / I(-3, -2) == Z); // Z/N
    CHECK(Z / I(-3, 2) == Z); // Z/M
    CHECK(Z / I(2, 3) == Z); // Z/P
    CHECK(I(-2, 0) / I(-5, -4) == I(0, 0.5)); // N0/N1
    CHECK(I(-2, 0) / Z == E); // N0/Z
    CHECK(I(-2, 0) / I(-4, 0) == I(0, inf)); // N0/N0
    CHECK(I(-3, 0) / I(-4, 2) == F); // N0/M
    CHECK(I(-2, 0) / I(0, 4) == I(-inf, 0)); // N0/P0
    CHECK(I(-2, 0) / I(4, 5) == I(-0.5, 0)); // N0/P1
    CHECK(I(-3, 2) / I(-5, -4) == I(-0.5, 0.75)); // M/N1
    CHECK(I(-3, 2) / Z == E); // M/Z
    CHECK(I(-3, 2) / I(-4, 0) == F); // M/N0
    CHECK(I(-3, 2) / I(-4, 5) == F); // M/M
    CHECK(I(-3, 2) / I(0, 5) == F); // M/P0
    CHECK(I(-3, 2) / I(4, 5) == I(-0.75, 0.5)); // M/P1
    CHECK(I(0, 2) / I(-5, -4) == I(-0.5, 0.0)); // P0/N1
    CHECK(I(0, 2) / Z == E); // P0/Z
    CHECK(I(0, 2) / I(-4, 0) == I(-inf, 0)); // P0/N0
    CHECK(I(0, 2) / I(-4, 5) == F); // P0/M
    CHECK(I(0, 2) / I(0, 5) == I(0, inf)); // P0/P0
    CHECK(I(0, 2) / I(4, 5) == I(0, 0.5)); // P0/P1
    CHECK(I(2, 3) / I(-5, -4) == I("-0.75", "-0.4")); // P1/N1
    CHECK(I(2, 3) / Z == E); // P1/Z
    CHECK(I(2, 3) / I(-4, 0) == I(-inf, -0.5)); // P1/N0
    CHECK(I(2, 3) / I(-4, 5) == F); // P1/M
    CHECK(I(2, 3) / I(0, 4) == I(0.5, inf)); // P1/P0
    CHECK_IA_APPROX(I(2, 3) / I(4, 5), I("0.4", "0.75")); // P1/P1
    CHECK_IA_APPROX(1 / I(10, 10), I("0.1"));
    CHECK(1 / I(10, 10) == I(10, 10).inv());

    CHECK(I(-30, -15) / I(-5, -3) == I(3, 10));
    CHECK_IA_APPROX(I(0.9, 2.0) / I(0.1, 1.1), I(0.8181818181818181, 20));
    CHECK(I(0.1, 1.1) / I(0.25, 4.0) == I(0.025, 4.4));
    CHECK(I(0.25, 4) / 4 == I(0.0625, 1));
    CHECK(I(0.25, 4) / Z == E);
    CHECK(I(0, 1) / I(0, 1) == I(0, fp::inf));
    CHECK(I(-1, 1) / I(0, 1) == F);
    CHECK(I(-1, 1) / I(-1, 1) == F);
}

TEST_CASE("inverse", "[IA]")
{
    CHECK(E.inv() == E);
    CHECK(Z.inv() == E);
    CHECK(F.inv() == F);
    CHECK_IA_APPROX(I(10, 10).inv(), I("0.1"));
    CHECK(I(0.5, 1.5).inv() == 1 / I(0.5, 1.5));
    CHECK(I(-4, -2).inv() == I(-0.5, -0.25));
    CHECK(I(-2, 0).inv() == I(-inf, -0.5));
    CHECK(I(-4, 2).inv() == F);
    CHECK(I(-2, 4).inv() == F);
    CHECK(I(0, 2).inv() == I(0.5, inf));
    CHECK(I(2, 4).inv() == I(0.25, 0.5));
}

TEST_CASE("optimize_bounds", "[IA]")
{
    pappus::box<double> x { I(0.0, 1.0), I(0.0, 2.0) };
    auto f = [](auto const& bx) { return bx[0] + bx[1]; };

    auto result = pappus::optimize_bounds(f, x, 1e-6, 64);
    CHECK(result.contains(I(0.0, 3.0)));
}

TEST_CASE("scalar base pow", "[IA]")
{
    CHECK(pow(2.0, I(2.0, 3.0)) == I(4.0, 8.0));
    CHECK(pow(4.0, I(-1.0, 1.0)) == I(0.25, 4.0));
    CHECK(pow(-2.0, I(0.5, 1.5)) == E);
}

TEST_CASE("trigonometric functions", "[IA]")
{
    SECTION("sin")
    {
        CHECK(F.sin() == I(-1, 1));
        CHECK(E.sin() == E);
        CHECK(Z.sin() == Z);
        CHECK(I(0, fp::pi).sin() == I(0, 1));
        CHECK_IA_APPROX(I(fp::pi, fp::pi).sin(), I(1.224646799147353e-16, 1.2246467991473535e-16));
        CHECK_IA_APPROX(I(0.5, 0.5).sin(), I(0.47942553860420295, 0.47942553860420306));
        CHECK(I(0.5, 1.67).sin() == I(0.47942553860420295, 1.0));
        CHECK(I(0.5, 8.5).sin() == I(-1, 1));
        CHECK(I(1.67, 3.2).sin() == I(-0.058374143427580093, 0.99508334981018032));
        CHECK(I(1.67,3.2).sin() == I(-0.058374143427580093, 0.99508334981018032));
        CHECK_IA_APPROX(I(2.1, 5.6).sin(), I(-1.0, 0.86320936664887382));
        CHECK(I(-4.5, 0.1).sin() == I(-1.0, 0.97753011766509712));
        CHECK(I(1.3, 6.3).sin() == I(-1.0, 1.0));
        CHECK(I(-1.0, -0.5).sin() == I(-0.84147098480789662, -0.47942553860420295));
        CHECK(I(-10, 10).sin() == I(-1, 1));
        CHECK(I(3, 3.5).sin() == I(-0.35078322768961989, 0.14112000805986724));
        CHECK(I(-3.5, -3).sin() == I(-0.14112000805986724, 0.35078322768961989));
        CHECK(I(-3.5, 3).sin() == I(-1, 1));
        CHECK(I(10, 12).sin() == I(-1, -0.53657291800043483));
        CHECK(I(13, 14).sin() == I(0.42016703682664086, 0.99060735569487046));
        CHECK(I(10, 14).sin() == I(-1, 0.99060735569487046));
        CHECK(I(14, 16).sin() == I(-0.28790331666506536, 1));
        CHECK_IA_APPROX(I(-11, -10).sin(), I(0.54402111088936966, 1));
        CHECK(I(-14, -13).sin() == I(-0.99060735569487046, -0.42016703682664086));
        CHECK(I(-16, -14).sin() == I(-1, 0.28790331666506536));
        CHECK(I(-102, -100).sin() == I(-0.99482679135840646, 0.5063656411097589));
        CHECK(I(4.6e15, 4.7e15).sin() == I(-1, 1));
        CHECK(I(4503599627370495, 4503599627370496).sin() == I(0.87421730262363495, 1.0));
    }

    SECTION("cos")
    {
        CHECK(F.cos() == I(-1, 1));
        CHECK(E.cos() == E);
        CHECK_IA_APPROX(Z.cos(), I(1, 1));
        CHECK_IA_APPROX(I(0.5).cos(), I(0.87758256189037265, 0.87758256189037287));
        CHECK_IA_APPROX(I(fp::pi, fp::pi).cos(), I(-1, -1));
        CHECK_IA_APPROX(I(0.5, 1.67).cos(), I(-0.099041036598728024, 0.87758256189037287));
        CHECK_IA_APPROX(I(2.1, 5.6).cos(), I(-1.0, 0.77556587851024972));
        CHECK(I(0.5, 8.5).cos() == I(-1.0, 1.0));
        CHECK_IA_APPROX(I(1.67, 3.2).cos(), I(-1.0, -0.099041036598727997));
        CHECK(I(-1.0, -0.5).cos() == I(0.54030230586813965, 0.87758256189037287));
        CHECK(I(-10, 10).cos() == I(-1, 1));
        CHECK(I(3, 3.5).cos() == I(-1, -0.93645668729079623));
        CHECK(I(-3.5, -3).cos() == I(-1, -0.93645668729079623));
        CHECK(I(-3.5, 3).cos() == I(-1, 1));
        CHECK_IA_APPROX(I(10, 12).cos(), I(-0.83907152907645255, 0.84385395873249225));
        CHECK(I(13, 14).cos() == I(0.13673721820783358, 0.9074467814501963));
        CHECK(I(10, 14).cos() == I(-0.83907152907645255, 1));
        CHECK(I(14, 16).cos() == I(-1, 0.13673721820783363));
        CHECK_IA_APPROX(I(-11, -10).cos(), I(-0.83907152907645255, 0.0044256979880507863));
        CHECK(I(-14, -13).cos() == I(0.13673721820783358, 0.9074467814501963));
        CHECK(I(-16, -14).cos() == I(-1, 0.13673721820783363));
        CHECK(I(-102, -100).cos() == I(0.10158570369662133, 1));
        CHECK(I(4.6e15, 4.7e15).cos() == I(-1, 1));
        CHECK(I(4503599627370495, 4503599627370496).cos() == I(-0.48553486774222065, 0.47329288595430913));
    }

    SECTION("tan")
    {
        CHECK_IA_APPROX(I(0.5).tan(), I(0.54630248984379037, 0.5463024898437906));
        CHECK(I(0.5, 1.67).tan() == F);
        CHECK(I(1.67, 3.2).tan() == I(-10.047182299210307, 0.058473854459578652));
        CHECK(I(6.638314112824137, 8.38263151220128).tan() == F);
        CHECK_IA_APPROX(I(0.0, 1.0).tan(), I(0.0, 1.5574077246549025));
        CHECK_IA_APPROX(I(-1.0, 0.0).tan(), I(-1.5574077246549025, 0.0));
        CHECK(I(-2.0, -1.0).tan() == F);
        CHECK(I(202, 203).tan() == F);
    }

    SECTION("asin")
    {
        CHECK(E.asin() == E);
        CHECK(I(-6, -3).asin() == E);
        CHECK(I(0.9, 2).asin() == I(0.9, 1).asin());
        CHECK(I(3, 4).asin() == E);
        CHECK(I(3, 5).asin() == E);
        CHECK(I(-1.5, -0.5).asin() == I(-1.5707963267948968, -0.52359877559829882));
        CHECK(I(-0.5, 0.5).asin() == I(-0.52359877559829904, 0.52359877559829904));
        // I("0.1") spans [nextafter(0.1,-inf), nextafter(0.1,+inf)] due to directed rounding of "0.1"
        CHECK(I("0.1").asin() == I(0x1.9a49276037882p-4, 0x1.9a49276037885p-4));
    }

    SECTION("acos")
    {
        CHECK(E.acos() == E);
        CHECK(I(-1.5, -0.5).acos() == I(2.0943951023931953, 3.1415926535897936));
        CHECK_IA_APPROX(I(-0.5,0.5).acos(), I(1.0471975511965976, 2.0943951023931962));
        CHECK(I("0.1").acos() == I(1.4706289056333366, 1.470628905633337));
    }

    SECTION("atan")
    {
        CHECK(E.atan() == E);
        // atan(0) = 0 exactly; the zero guard in cr_floor/cr_ceil preserves it
        CHECK(Z.atan() == Z);
        // atan(±∞) = ±π/2; with nextafter the result is 1 ULP beyond ±π/2
        CHECK(!F.atan().is_infinite());
        CHECK_IA_APPROX(I(1.0).atan(), I(0x1.921fb54442d17p-1, 0x1.921fb54442d19p-1));
        CHECK_IA_APPROX(I(-1.0).atan(), I(-0x1.921fb54442d19p-1, -0x1.921fb54442d17p-1));
        CHECK_IA_APPROX(I(0.5).atan(), I(0x1.dac670561bb4ep-2, 0x1.dac670561bb5p-2));
        CHECK_IA_APPROX(I(-0.5).atan(), I(-0x1.dac670561bb5p-2, -0x1.dac670561bb4ep-2));
        CHECK_IA_APPROX(I(2.0).atan(), I(0x1.1b6e192ebbe43p+0, 0x1.1b6e192ebbe45p+0));
        CHECK_IA_APPROX(I(-2.0).atan(), I(-0x1.1b6e192ebbe45p+0, -0x1.1b6e192ebbe43p+0));
        CHECK_IA_APPROX(I(0.5, 1.0).atan(), I(0x1.dac670561bb4ep-2, 0x1.921fb54442d19p-1));
        CHECK(I(-1.0, 1.0).atan() == I(-0x1.921fb54442d19p-1, 0x1.921fb54442d19p-1));
        CHECK_IA_APPROX(I(-10.0, 10.0).atan(), I(-0x1.789bd2c160055p+0, 0x1.789bd2c160055p+0));
        // atan is monotone — result must be contained in (-π/2, π/2)
        CHECK(I(0.0, 1.0).atan().inf() == 0.0);
        CHECK(I(-1.0, 0.0).atan().sup() == 0.0);
    }
}

TEST_CASE("hyperbolic functions", "[IA]")
{
    SECTION("sinh")
    {
        CHECK(E.sinh() == E);
        CHECK_IA_APPROX(I(0.5).sinh(), I(0.5210953054937473, 0.5210953054937474));
        CHECK_IA_APPROX(I(0.5, 1.67).sinh(), I(0.5210953054937473, 2.5619603657712102));
        CHECK(I(-4.5, 0.1).sinh() == I(-45.00301115199179, 0.10016675001984404));
        // negative-only intervals: verifies both lo and hi are sound
        CHECK(I(-1.0, -0.5).sinh() == I(-1.1752011936438016004, -0.52109530549374727393));
        CHECK(I(-2.0, -0.5).sinh() == I(-3.6268604078470194629, -0.52109530549374727393));
    }

    SECTION("cosh")
    {
        CHECK(E.cosh() == E);
        CHECK_IA_APPROX(I(0.5).cosh(), I(1.1276259652063807, 1.127625965206381));
        CHECK_IA_APPROX(I(0.5, 1.67).cosh(), I(1.1276259652063807, 2.750207431409957));
        CHECK_IA_APPROX(I(-4.5, 0.1).cosh(), I(1.0, 45.01412014853004));
    }

    SECTION("tanh")
    {
        // glibc's tanh ignores fesetround entirely; we use nextafter of nearest.
        CHECK(E.tanh() == E);
        CHECK(Z.tanh() == Z);  // tanh(0) = 0 exactly, zero guard preserves it
        CHECK_IA_APPROX(I(0.5, 1.67).tanh(), I(0x1.d9353d7568af2p-2, 0x1.dcf457a7e97c4p-1));
        CHECK_IA_APPROX(I(-4.5, 0.1).tanh(), I(-0x1.ffdfa72153984p-1, 0x1.983d7795f413bp-4));
        // negative-only intervals: previously broken by fmin/fmax workaround
        CHECK_IA_APPROX(I(-1.0, -0.5).tanh(), I(-0x1.85efab514f395p-1, -0x1.d9353d7568af2p-2));
        CHECK_IA_APPROX(I(-2.0, -0.5).tanh(), I(-0x1.ed9505e1bc3d5p-1, -0x1.d9353d7568af2p-2));
    }
}

TEST_CASE("trancendental functions", "[IA]")
{
    SECTION("exp")
    {
        CHECK(E.exp() == E);
        CHECK_IA_APPROX(Z.exp(), I(1));

        CHECK_IA_APPROX(I(1/2.).exp(), I(1.648721270700128, 1.6487212707001282));
        CHECK_IA_APPROX(I(0.1).exp(), I(1.1051709180756475e+00, 1.1051709180756477e+00));
        CHECK_IA_APPROX(I(1).exp(), I(2.718281828459045, 2.7182818284590455));
        CHECK_IA_APPROX(I(999).exp(), I(max, fp::inf));
        CHECK_IA_APPROX(I(-999).exp(), I(0, 4.94066e-324)); 
        CHECK_IA_APPROX(I(-10, -5).exp(), I(4.539992976248485e-5, 0.006737946999085468));
        CHECK_IA_APPROX(I(-5, 9).exp(), I(0.006737946999085467, 8103.083927575384));
        CHECK(I(9.0,11.0).exp() == I(8103.083927575383, 59874.141715197824));
        CHECK_IA_APPROX(I(-3.5).exp(), I(0.030197383422318497, 0.0301973834223185));
        CHECK_IA_APPROX(I(3.5).exp(), I(33.11545195869231, 33.11545195869232));
    }

    SECTION("log")
    {
        CHECK(E.log() == E);
        CHECK_IA_APPROX(I(1/2.).log(), I(-6.931471805599454e-01, -6.9314718055994529e-01));
        CHECK_IA_APPROX(I(0.1).log(), I(-2.3025850929940459e+00, -2.3025850929940455e+00));
        CHECK(I(-4, 0).log() == E); 
    }
}

TEST_CASE("power functions", "[IA]")
{
    SECTION("integer power")
    {
        CHECK(I(0, 3).pow(2) == I(0, 9));
        CHECK(I(2, 3).pow(2) == I(4, 9));
        CHECK(I(-3, 0).pow(2) == I(0, 9));
        CHECK(I(-3, -2).pow(2) == I(4, 9));
        CHECK(I(-3, 2).pow(2) == I(0, 9));

        CHECK(I(0, 3).pow(2) == I(0, 3).square());
        CHECK(I(2, 3).pow(2) == I(2, 3).square());
        CHECK(I(-3, 0).pow(2) == I(-3, 0).square());
        CHECK(I(-3, -2).pow(2) == I(-3, -2).square());
        CHECK(I(-3, 2).pow(2) == I(-3, 2).square());

        CHECK(I(0, 3).pow(3) == I(0, 27));
        CHECK(I(2, 3).pow(3) == I(8, 27));
        CHECK(I(-3, 0).pow(3) == I(-27, 0));
        CHECK(I(-3, -2).pow(3) == I(-27, -8));
        CHECK(I(-3, 2).pow(3) == I(-27, 8));

        CHECK_IA_APPROX(I(0, 3).pow(-2), I(1/9., fp::inf));
        CHECK_IA_APPROX(I(-3, 0).pow(-2), I(1/9., fp::inf));
        CHECK_IA_APPROX(I(-3, 2).pow(-2), I(1/9., fp::inf));
        CHECK_IA_APPROX(I(2, 3).pow(-2), I(1/9., 1./4));
        CHECK(I(1, 2).pow(-3) == I(1/8., 1.));
        CHECK_IA_APPROX(I(0, 3).pow(-3), I(1/27., fp::inf));
        CHECK(I(-1, 2).pow(-3) == F);

        CHECK(E.pow(0) == E);
    }

    SECTION("double power")
    {
        CHECK(Z.pow(1.1) == Z);
        CHECK(Z.pow(0.0) == E);
        CHECK(Z.pow(1/10.) == Z);
        CHECK(Z.pow(-1/10.) == E);
        CHECK(E.pow(0.0) == E);
        CHECK(I(2.5).pow(3) == I(15.625));
        CHECK(I(-3,4).pow(0.5) == I(0, 2));
        CHECK(I(-3,4).pow(0.5) == I(-3,4).pow(1/2.));
    }
}
