#ifndef PAPPUS_FPOPS_HPP
#define PAPPUS_FPOPS_HPP

#include "util.hpp"
#include <type_traits>
#if defined(USE_CRLIBM)
#include <crlibm.h>
#endif

namespace pappus {
namespace fp {
    // binary operations
    using op_add = std::plus<double>;
    using op_sub = std::minus<double>;
    using op_mul = std::multiplies<double>;
    using op_div = std::divides<double>;
    struct op_pow { double operator()(double a, double b) { return std::pow(a, b); } };
    struct op_abs { double operator()(double a) { return std::fabs(a); } };

    // unary operations
    struct op_exp  { double operator()(double a) { return std::exp(a); } }; 
    struct op_log  { double operator()(double a) { return std::log(a); } }; 
    struct op_sin  { double operator()(double a) { return std::sin(a); } }; 
    struct op_cos  { double operator()(double a) { return std::cos(a); } }; 
    struct op_tan  { double operator()(double a) { return std::tan(a); } }; 
    struct op_asin { double operator()(double a) { return std::asin(a); } };
    struct op_acos { double operator()(double a) { return std::acos(a); } };
    struct op_atan { double operator()(double a) { return std::atan(a); } };
    struct op_fmod { double operator()(double a, double b) { return std::fmod(a, b); } };
    struct op_sinh { double operator()(double a) { return std::sinh(a); } };
    struct op_cosh { double operator()(double a) { return std::cosh(a); } };
    struct op_tanh { double operator()(double a) { return std::tanh(a); } };
    struct op_sqrt { double operator()(double a) { return std::sqrt(a); } };

    template<int ROUND_MODE, typename OP, typename... Args>
    double rop(Args&& ...args)
    {
        static_assert(std::is_invocable_r_v<double, OP, Args...>);
#if !defined(DIRECTED_ROUNDING)
        return OP()(args...);
#else
        static_assert(ROUND_MODE == FE_UPWARD || ROUND_MODE == FE_DOWNWARD);
        auto rnd = std::fegetround();
        std::fesetround(ROUND_MODE);
        auto c = OP()(args...);
        std::fesetround(rnd);
        return c;
#endif
    }

#if defined(USE_CRLIBM)
    template<typename OP = op_add>
    double ropd(double a, double b) { return rop<FE_DOWNWARD, OP>(a, b); }
    template<> double ropd<op_sub>(double a, double b) { return rop<FE_DOWNWARD, op_sub>(a, b); }
    template<> double ropd<op_mul>(double a, double b) { return rop<FE_DOWNWARD, op_mul>(a, b); }
    template<> double ropd<op_div>(double a, double b) { return rop<FE_DOWNWARD, op_div>(a, b); }
    template<> double ropd<op_fmod>(double a, double b) { return rop<FE_DOWNWARD, op_div>(a, b); }
    template<> double ropd<op_pow>(double a, double b) { return rop<FE_DOWNWARD, op_pow>(a, b); }

    template<typename OP = op_exp> double ropd(double a) { return exp_rd(a); }
    template<> double ropd<op_log>(double a) { return log_rd(a); }
    template<> double ropd<op_sin>(double a) { return sin_rd(a); }
    template<> double ropd<op_cos>(double a) { return cos_rd(a); }
    template<> double ropd<op_tan>(double a) { return tan_rd(a); }
    template<> double ropd<op_asin>(double a) { return asin_rd(a); }
    template<> double ropd<op_acos>(double a) { return acos_rd(a); }
    template<> double ropd<op_atan>(double a) { return atan_rd(a); }
    template<> double ropd<op_sinh>(double a) { return sinh_rd(a); }
    template<> double ropd<op_cosh>(double a) { return cosh_rd(a); }
    template<> double ropd<op_tanh>(double a) { return sinh_rd(a) / cosh_rd(a); }
    template<> double ropd<op_sqrt>(double a) { return rop<FE_DOWNWARD, op_sqrt>(a); }
    template<> double ropd<op_abs>(double a) { return rop<FE_DOWNWARD, op_abs>(a); }

    template<typename OP = op_add>
    double ropu(double a, double b) { return rop<FE_UPWARD, OP>(a, b); }
    template<> double ropu<op_sub>(double a, double b) { return rop<FE_UPWARD, op_sub>(a, b); }
    template<> double ropu<op_mul>(double a, double b) { return rop<FE_UPWARD, op_mul>(a, b); }
    template<> double ropu<op_div>(double a, double b) { return rop<FE_UPWARD, op_div>(a, b); }
    template<> double ropu<op_fmod>(double a, double b) { return rop<FE_UPWARD, op_div>(a, b); }
    template<> double ropu<op_pow>(double a, double b) { return rop<FE_UPWARD, op_pow>(a, b); }

    template<typename OP = op_exp> double ropu(double a) { return exp_ru(a); }
    template<> double ropu<op_log>(double a) { return log_ru(a); }
    template<> double ropu<op_sin>(double a) { return sin_ru(a); }
    template<> double ropu<op_cos>(double a) { return cos_ru(a); }
    template<> double ropu<op_tan>(double a) { return tan_ru(a); }
    template<> double ropu<op_asin>(double a) { return asin_ru(a); }
    template<> double ropu<op_acos>(double a) { return acos_ru(a); }
    template<> double ropu<op_atan>(double a) { return atan_ru(a); }
    template<> double ropu<op_sinh>(double a) { return sinh_ru(a); }
    template<> double ropu<op_cosh>(double a) { return cosh_ru(a); }
    template<> double ropu<op_tanh>(double a) { return sinh_ru(a) / cosh_rd(a); }
    template<> double ropu<op_sqrt>(double a) { return rop<FE_UPWARD, op_sqrt>(a); }
    template<> double ropu<op_abs>(double a) { return rop<FE_UPWARD, op_abs>(a); }
#else    
    template<typename OP, typename... Args>
    constexpr auto ropd = [](auto&& ...args) { return rop<FE_DOWNWARD, OP, Args...>(args...); };

    template<typename OP, typename... Args>
    constexpr auto ropu = [](auto&& ...args) { return rop<FE_UPWARD, OP, Args...>(args...); };
#endif

    namespace trig {
        const auto get_quadrant = [](auto a) { 
            static_assert(std::is_floating_point_v<decltype(a)>);
            auto x = std::fmod(a, two_pi);
            if (x == 0) return 0;
            if (x < 0) x += two_pi; 
            return (int)(x / half_pi);
        };
    }
}
}

#endif
