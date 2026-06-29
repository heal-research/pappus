#include <catch2/catch_all.hpp>

#include "pappus.hpp"

TEST_CASE("ops layer")
{
    using namespace pappus;

    auto ix = ops::variable<float>(-1.0f, 1.0f);
    auto iy = ops::constant<float>(2.0f);
    CHECK(ops::add(ix, iy) == interval<float>(1.0f, 3.0f));

    ops::affine_context<float> ctx;
    ctx.max_terms = 4;
    auto ax = ops::variable(ctx, -1.0f, 1.0f);
    auto ay = ops::constant(ctx, 2.0f);
    auto az = ops::mul(ctx, ax, ay);

    CHECK(az.to_interval() == interval<float>(-2.0f, 2.0f));
    CHECK(ctx.approximation_mode() == approximation_mode::CHEBYSHEV);
}
