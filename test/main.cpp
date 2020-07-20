#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 
#include <doctest/doctest.h>

#include "pappus.hpp"

#include <Eigen/Eigen>

TEST_CASE("Eigen") 
{
    Eigen::ArrayXd a(3);
    a << 1, 2, 3;
    
    std::cout << a << "\n";

    double v = a(a.size()-1);
    std::cout << v << "\n";
}
