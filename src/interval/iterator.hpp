#ifndef PAPPUS_INTERVAL_ITERATOR_HPP
#define PAPPUS_INTERVAL_ITERATOR_HPP

#include <iterator>
#include <tuple>

#include "fp/util.hpp"

namespace pappus {
class interval;

class iterator {
public:
    using value_type        = interval;
    using reference         = value_type const&;
    using pointer           = value_type const*;
    using difference_type   = std::size_t;
    using iterator_category = std::forward_iterator_tag;

    iterator& operator++() // preincrement
    {
        ++i_;
        a_ = b_;
        b_ = (i_ + 1) * h_;
        return *this;
    }
    iterator operator++(int) // postincrement
    {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    interval operator*() const;

    bool operator==(iterator other) const;
    bool operator!=(iterator other) const { return !(*this == other); }
    bool operator<(iterator other) const;

    iterator(interval const&, size_t);

    iterator begin();
    iterator end();

private:
    void advance_to_end() 
    {
        i_ = n_;
        a_ = fp::nan;
        b_ = fp::nan;
    }

    interval const& iv_; // interval
    std::size_t i_; // current index
    std::size_t const n_; // final index 
    double const h_; // increment
    double a_; // running lower bound
    double b_; // running upper bound
};

// useful as a kind of view over the interval
class subdivision {
public:
    iterator begin() const;
    iterator end() const;

    subdivision(interval const& iv, std::size_t n) : iv_(iv), n_(n) {}

private:
    interval const& iv_;
    std::size_t const n_;
};

} // namespace
#endif
