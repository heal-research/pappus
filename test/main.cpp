#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "aa.h"
#include "pappus.hpp"

#include <Eigen/Eigen>

using af = pappus::affine_form;
using ai = pappus::affine_interval;

bool operator==(ai const& lhs, AAInterval const& rhs)
{
    return lhs.lower() == rhs.getlo()
        && lhs.mid() == rhs.mid()
        && lhs.upper() == rhs.gethi()
        && lhs.radius() == rhs.radius();
}
bool operator==(AAInterval const& lhs, ai const& rhs)
{
    return rhs == lhs;
}

TEST_CASE("affine addition")
{

    ai u1(-2, 3);
    ai v1(-1, 1);

    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1(ctx, v1);

    af z1 = x1 + y1;

    std::cout << "=== pappus ===\n";
    std::cout << "u = " << u1 << "\n";
    std::cout << "v = " << v1 << "\n";

    std::cout << "x = " << x1.to_interval() << "\n";
    std::cout << "y = " << y1.to_interval() << "\n";
    std::cout << "z = " << z1.to_interval() << "\n";

    AAInterval u2(-2, 3);
    AAInterval v2(-1, 1);

    AAF x2(u2);
    AAF y2(v2);
    AAF z2 = x2 + y2;

    std::cout << "=== aaflib ===\n";
    std::cout << "u = " << u2;
    std::cout << "v = " << v2;

    std::cout << "x = " << x2.convert();
    std::cout << "y = " << y2.convert();
    std::cout << "z = " << z2.convert();

    CHECK(z1.to_interval() == z2.convert());
}

TEST_CASE("affine multiplication")
{

    ai u1(-2, 3);
    ai v1(-1, 1);

    pappus::affine_context ctx;
    af x1(ctx, u1);
    af y1(ctx, v1);

    af z1 = x1 * y1;

    std::cout << "=== pappus ===\n";
    std::cout << "u = " << u1 << "\n";
    std::cout << "v = " << v1 << "\n";

    std::cout << "x = " << x1.to_interval() << "\n";
    std::cout << "y = " << y1.to_interval() << "\n";
    std::cout << "z = " << z1.to_interval() << "\n";

    AAInterval u2(-2, 3);
    AAInterval v2(-1, 1);

    AAF x2(u2);
    AAF y2(v2);
    AAF z2 = x2 * y2;

    std::cout << "=== aaflib ===\n";
    std::cout << "u = " << u2;
    std::cout << "v = " << v2;

    std::cout << "x = " << x2.convert();
    std::cout << "y = " << y2.convert();
    std::cout << "z = " << z2.convert();

    CHECK(z1.to_interval() == z2.convert());
}
