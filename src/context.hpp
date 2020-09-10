#ifndef PAPPUS_AFFINE_CONTEXT_HPP
#define PAPPUS_AFFINE_CONTEXT_HPP

#include <cstdint>
#include <functional>
#include <vector>

namespace pappus {

class affine_form; // forward declaration

enum class approximation_mode : int {
    CHEBYSHEV, // (default)
    MINRANGE,
    SECANT
};

// a context object is necessary to keep track of shared state between
// a collection of affine_forms that interact together
class affine_context {
public:
    affine_context()
        : last_index_(0)
        , approximation_mode_(approximation_mode::CHEBYSHEV)
    {
    }

    std::size_t last_index() const
    {
        return last_index_;
    }

    void set_last_index(std::size_t last_index)
    {
        last_index_ = last_index;
    }
    // increase the highest symbol
    std::size_t increment_last()
    {
        return ++last_index_;
    }

    enum approximation_mode approximation_mode() const
    {
        return approximation_mode_;
    }

private:
    std::size_t last_index_;
    enum approximation_mode approximation_mode_;
};
}

#endif
