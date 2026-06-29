#ifndef PAPPUS_AFFINE_CONTEXT_HPP
#define PAPPUS_AFFINE_CONTEXT_HPP

#include <atomic>
#include <cstddef>
#include <memory>

namespace pappus {

enum class approximation_mode : int {
    CHEBYSHEV, // (default)
    MINRANGE,
    SECANT
};

// keeps shared state (symbol counter, approximation mode) for a collection
// of affine_form objects that interact together
class affine_context {
public:
    affine_context()
        : state_(std::make_shared<state>())
    {
    }

    std::size_t last_index() const { return state_->last_index.load(); }

    void set_last_index(std::size_t last_index) const { state_->last_index.store(last_index); }

    std::size_t increment_last() const
    {
        return state_->last_index.fetch_add(1);
    }

    enum approximation_mode approximation_mode() const { return state_->approximation_mode.load(); }

    void set_approximation_mode(enum approximation_mode mode) const { state_->approximation_mode.store(mode); }

    bool shares_state_with(affine_context const& other) const { return state_ == other.state_; }

private:
    struct state {
        std::atomic<std::size_t> last_index {0};
        std::atomic<enum approximation_mode> approximation_mode {approximation_mode::CHEBYSHEV};
    };

    std::shared_ptr<state> state_;
};

} // namespace pappus

#endif
