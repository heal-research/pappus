#ifndef PAPPUS_FPUTIL_HPP
#define PAPPUS_FPUTIL_HPP

#include <cfenv>
#include <cmath>
#include <sstream>

#define EXPECT(cond)                                                                                    \
    if (!(cond)) {                                                                                      \
        std::cout << "precondition " << #cond << " failed at " << __FILE__ << ": " << __LINE__ << "\n"; \
        std::terminate();                                                                               \
    }

#define ENSURE(cond)                                                                                     \
    if (!(cond)) {                                                                                       \
        std::cout << "postcondition " << #cond << " failed at " << __FILE__ << ": " << __LINE__ << "\n"; \
        std::terminate();                                                                                \
    }

#ifndef M_PI
#define M_PI 3.14159265358979323846 // pi
#endif

namespace pappus {
namespace fp {
    const auto inf = std::numeric_limits<double>::infinity();
    const auto nan = std::numeric_limits<double>::quiet_NaN();

    using op_add = std::plus<double>;
    using op_sub = std::minus<double>;
    using op_mul = std::multiplies<double>;
    using op_div = std::divides<double>;

    template <int ROUND_MODE>
    double from_string(std::string const& s)
    {
        static_assert(ROUND_MODE == FE_UPWARD || ROUND_MODE == FE_DOWNWARD);
        auto rounding_mode = std::fegetround();
        std::fesetround(ROUND_MODE);
        std::istringstream is(s);
        double v;
        is >> v;
        std::fesetround(rounding_mode);
        return v;
    }

    // rounded op
    template <typename OP, int ROUND_MODE>
    double rop(double a, double b)
    {
        static_assert(ROUND_MODE == FE_UPWARD || ROUND_MODE == FE_DOWNWARD);
        static_assert(std::is_invocable_r<double, OP, double, double>::value);
        auto rnd = std::fegetround();
        std::fesetround(ROUND_MODE);
        auto c = OP()(a, b);
        std::fesetround(rnd);
        return c;
    }

    // rounded op, downwards
    template <typename OP>
    const auto ropd = rop<OP, FE_DOWNWARD>;

    // rounded op, upwards
    template <typename OP>
    const auto ropu = rop<OP, FE_UPWARD>;
}
}

#endif
