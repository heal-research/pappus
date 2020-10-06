#include "interval.hpp"
#include "iterator.hpp"

#include <algorithm>

namespace pappus {

iterator::iterator(interval const& iv, size_t n) :
    iv_(iv), i_(0ul), n_(n), h_(iv_.radius() / n_), a_(iv_.inf()), b_(a_ + h_)
{
}

interval iterator::operator*() const
{
    return interval(a_, b_);
}

iterator subdivision::begin() const
{ 
    return iterator(iv_, n_).begin();
}

iterator subdivision::end() const
{
    return iterator(iv_, n_).end();
}

iterator iterator::begin()
{
    return iterator(iv_, n_);
}

iterator iterator::end()
{
    auto tmp = *this;
    tmp.advance_to_end();
    return tmp;
}

bool iterator::operator==(iterator other) const
{
    if (iv_ != other.iv_) {
        return false;
    }

    return i_ == other.i_; // cheating for speed
    //return std::tie(i_, n_, h_, a_, b_) == std::tie(other.i_, other.n_, other.h_, other.a_, other.b_);
}

bool iterator::operator<(iterator other) const
{
    if (iv_ != other.iv_) {
        return false;
    }

    return i_ < other.i_; // cheating for speed
}
} // namespace

