#ifndef PAPPUS_AFFINE_CONTEXT_HPP
#define PAPPUS_AFFINE_CONTEXT_HPP

#include <cstdint>
#include <vector>
#include <functional>

namespace pappus {
// a context object is necessary to keep track of shared state between
// a collection of affine_forms that interact together

class affine_form; // forward declaration
    
class affine_context {
public:
    affine_context() {}

    void add(affine_form& af) { affines_.push_back(std::ref(af)); }

    void clear() { 
        affines_.clear();
        last_index_ = 0;
    }

    std::size_t last_index() const { return last_index_; }
    void set_last_index(std::size_t last_index) { last_index_ = last_index; }

    // increase the highest symbol
    std::size_t increment_last() { return ++last_index_; }

    affine_form& operator[](std::size_t i) { return affines_[i]; }

private:
    std::size_t last_index_;
    std::vector<std::reference_wrapper<affine_form>> affines_;
};
}

#endif
