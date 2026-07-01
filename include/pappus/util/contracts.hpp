#ifndef PAPPUS_CONTRACTS_HPP
#define PAPPUS_CONTRACTS_HPP

#include <iostream>

#define EXPECT(cond)                                                                                    \
    if (!(cond)) {                                                                                      \
        std::cerr << "precondition " << #cond << " failed at " << __FILE__ << ": " << __LINE__ << "\n"; \
        std::terminate();                                                                               \
    }

#define ENSURE(cond)                                                                                     \
    if (!(cond)) {                                                                                       \
        std::cerr << "postcondition " << #cond << " failed at " << __FILE__ << ": " << __LINE__ << "\n"; \
        std::terminate();                                                                                \
    }

#endif
