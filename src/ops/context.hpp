#ifndef PAPPUS_OPS_CONTEXT_HPP
#define PAPPUS_OPS_CONTEXT_HPP

#include <concepts>
#include <cstddef>

#include "affine/context.hpp"

namespace pappus::ops {

template<std::floating_point T>
struct affine_context {
    pappus::affine_context state {};
    std::size_t            max_terms {0};

    affine_context() = default;

    explicit affine_context(pappus::approximation_mode mode, std::size_t max_terms = 0)
        : max_terms(max_terms)
    {
        state.set_approximation_mode(mode);
    }

    pappus::approximation_mode approximation_mode() const
    {
        return state.approximation_mode();
    }

    void set_approximation_mode(pappus::approximation_mode mode) const
    {
        state.set_approximation_mode(mode);
    }
};

} // namespace pappus::ops

#endif
