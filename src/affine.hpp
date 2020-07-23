#ifndef PAPPUS_AF_HPP
#define PAPPUS_AF_HPP

#include "context.hpp"
#include "interval.hpp"

#include <vector>

#include <Eigen/Eigen>

namespace pappus {

template <typename T>
using array = Eigen::Array<T, Eigen::Dynamic, 1>;

namespace {
template<typename T>
static auto eigen_array(std::vector<T> const& vec)
{
    return Eigen::Map<const array<T>>(vec.data(), vec.size());
}

template<typename T>
static auto eigen_array(std::vector<T>& vec)
{
    return Eigen::Map<array<T>>(vec.data(), vec.size());
}
}

class affine_form {
public:
    explicit affine_form(affine_context& context, affine_interval const& interval)
        : context_(context)
        , center_((interval.upper() + interval.lower()) / 2)
        , radius_(0)
        , deviations_(1)
        , indices_(1)
        , length_(1)
    {
        deviations_[0] = (interval.upper() - interval.lower()) / 2;
        indices_[0] = this->context().increment_last(); // assign and increment global index
        update_radius();
    }

    explicit affine_form(affine_context& context, double v)
        : context_(context)
        , center_(v)
        , radius_(0)
    {
    }

    explicit affine_form(affine_context& context, double v, std::vector<double> const& deviations, std::vector<size_t> const& indices, size_t len)
        : context_(context)
        , center_(v)
        , radius_(0)
        , deviations_(deviations)
        , indices_(indices)
        , length_(len)
    {
        update_radius();
    }

    // not explicit because we want to allow implicit copy
    affine_form(const affine_form& other)
        : context_(other.context_)
        , center_(other.center_)
        , radius_(other.radius_)
        , deviations_(other.deviations_)
        , indices_(other.indices_)
        , length_(other.length_)
    {
    }

    void swap(affine_form& other) 
    {
        std::swap(center_, other.center_);
        std::swap(radius_, other.radius_);
        std::swap(length_, other.length_);
        deviations_.swap(other.deviations_);
        indices_.swap(other.indices_);
    }

    affine_context& context() const { return context_.get(); }

    double center() const { return center_; }
    double radius() const { return radius_; }

    double max() const { return center() + radius(); }
    double min() const { return center() - radius(); }

    double abs_max() const
    {
        return std::fmax(std::fabs(min()), std::fabs(max()));
    }

    double abs_min() const
    {
        auto min_ = min();
        auto max_ = max();

        auto abs_min_ = std::fabs(min_);
        auto abs_max_ = std::fabs(max_);

        return std::signbit(min_) == std::signbit(max_)
            ? std::fmin(abs_min_, abs_max_)
            : 0;
    }

    size_t size() const { return indices_.size(); }

    size_t length() const { return length_; }
    size_t last_index() const { return indices_[length_-1]; }

    // comparison operators
    bool operator<(const affine_form& other)
    {
        return max() < other.max();
    }

    bool operator<=(const affine_form& other)
    {
        return max() <= other.max();
    }

    bool operator>(const affine_form& other)
    {
        return !(*this <= other);
    }

    bool operator>=(const affine_form& other)
    {
        return !(*this < other);
    }

    bool operator==(const affine_form& other) const;
    // interval representation
    affine_interval to_interval() const
    {
        return affine_interval(min(), max());
    }

    // assignment operators
    affine_form& operator=(double);
    affine_form& operator=(const affine_form&);
    affine_form& operator=(affine_form);

    // arithmetic operations
    affine_form& operator+=(double);
    affine_form& operator-=(double);
    affine_form& operator*=(double);
    affine_form& operator/=(double);
    affine_form operator+(double) const;
    affine_form operator*(double) const;
    affine_form operator-(double) const;
    affine_form operator/(double) const;
    affine_form operator-() const;
    affine_form operator^(int) const;

    affine_form operator+(affine_form const&) const;
    affine_form operator-(affine_form const&) const;
    affine_form operator*(affine_form const&) const;
    affine_form operator/(affine_form const&) const;
    affine_form operator^(affine_form const&) const;

    affine_form& operator+=(affine_form const&);
    affine_form& operator-=(affine_form const&);
    affine_form& operator*=(affine_form const&);
    affine_form& operator/=(affine_form const&);

    // other operations
    affine_form inv() const;

    // friends
    friend affine_form operator+(double v, affine_form const& af) { return af + v; }
    friend affine_form operator-(double v, affine_form const& af) { return af - v; }
    friend affine_form operator*(double v, affine_form const& af) { return af * v; }
    friend affine_form operator/(double v, affine_form const& af) { return af / v; }

    friend std::ostream& operator<<(std::ostream& s, affine_form& af)
    {
        s << "-------------------\n";
        s << "center: " << af.center_ << "\n";
        s << "radius: " << af.radius_ << "\n";
        s << "deviations: " << eigen_array(af.deviations_).transpose() << "\n";
        s << "indices: " << eigen_array(af.indices_).transpose() << "\n";
        s << "length: " << af.length_ << "\n";
        s << "-------------------\n";
        return s;
    }

private : 
std::reference_wrapper<affine_context> context_;
double center_; // central value
double radius_; // affine radius

std::vector<double> deviations_; // vector of partial deviations
std::vector<size_t> indices_; // vector of indices

size_t length_; // the actual number of indices still in use by this affine form

void update_radius()
{
    radius_ = eigen_array(deviations_).segment(0, length_).abs().sum();
}

};

} // namespace pappus

#endif
