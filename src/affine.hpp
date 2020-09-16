#ifndef PAPPUS_AF_HPP
#define PAPPUS_AF_HPP

#include "context.hpp"
#include "interval.hpp"

#include <stdexcept>
#include <vector>

#include <Eigen/Eigen>

#ifndef PROTECTED_DIVISION
#define PROTECTED_DIVISION 1
#endif

namespace pappus {

template <typename T>
using array = Eigen::Array<T, Eigen::Dynamic, 1>;

template <typename T>
using mat = Eigen::Matrix<T, Eigen::Dynamic, 1>;

struct view {
    template <typename T>
    static auto as_array(std::vector<T> const& vec)
    {
        return Eigen::Map<const array<T>>(vec.data(), vec.size());
    }

    template <typename T>
    static auto as_array(std::vector<T>& vec)
    {
        return Eigen::Map<array<T>>(vec.data(), vec.size());
    }

    template <typename T>
    static auto as_vector(std::vector<T> const& vec) {
        return Eigen::Map<const mat<T>>(vec.data(), vec.size());
    }

    template <typename T>
    static auto as_vector(std::vector<T>& vec) {
        return Eigen::Map<mat<T>>(vec.data(), vec.size());
    }
};

// TODO: 
// in my opinion using an epsilon does not really protect from anything or bring more clarity
// maybe it would be best to avoid "protected division" and just deal with NaN's in a sensible way
struct limits {
    static constexpr double eps = std::numeric_limits<double>::epsilon();
    static constexpr double minrad = 1e-10;
};

namespace {
    enum opcode { add,
        multiply };

    template <opcode x = opcode::add>
    struct binary_op {
        template <typename T, typename U>
        auto operator()(T t, U u) { return t + u; }
    };

    template <>
    struct binary_op<opcode::multiply> {
        template <typename T, typename U>
        auto operator()(T t, U u) { return t * u; }
    };
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
        , length_(0)
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

    template <typename OP>
    static std::optional<affine_form> handle_special_cases(affine_form const& lhs, affine_form const& rhs)
    {
        if (lhs.length() == 0 && rhs.length() == 0) {
            return affine_form(lhs.context(), OP()(lhs.center(), rhs.center()));
        }

        if (lhs.length() == 0) {
            return OP()(rhs, lhs.center());
        }

        if (rhs.length() == 0) {
            return OP()(lhs, rhs.center());
        }
        return std::nullopt;
    }

    // replicates method from aaflib
    double operator[](size_t i) const {
        return i == 0 ? center_ : deviations_[i-1];
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

        if (std::signbit(min_) == std::signbit(max_))
            return std::fmin(std::fabs(min_),
                             std::fabs(max_));

        return 0;
    }

    size_t size() const { return indices_.size(); }

    size_t length() const { return length_; }
    size_t last_index() const { return indices_[length_ - 1]; }

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
    affine_form pow(int) const;    // TODO:
    affine_form pow(double) const; // decide if these two need to be separate

    affine_form operator+(affine_form const&) const;
    affine_form operator-(affine_form const&) const;
    affine_form operator*(affine_form const&) const;
    affine_form operator/(affine_form const&) const;
    affine_form pow(affine_form const&) const;

    affine_form& operator+=(affine_form const&);
    affine_form& operator-=(affine_form const&);
    affine_form& operator*=(affine_form const&);
    affine_form& operator/=(affine_form const&);

    // other operations
    affine_form inv() const;

    // friends
    friend affine_form operator+(double v, affine_form const& af) { return af + v; }
    friend affine_form operator-(double v, affine_form const& af) { return -af + v; }
    friend affine_form operator*(double v, affine_form const& af) { return af * v; }
    friend affine_form operator/(double v, affine_form const& af) { return af.inv() * v; }

    friend std::ostream& operator<<(std::ostream& s, affine_form const& af)
    {
        s << "-------------------\n";
        s << "center:     " << af.center_ << "\n";
        s << "radius:     " << af.radius_ << "\n";
        s << "deviations: " << view::as_array(af.deviations_).transpose() << "\n";
        s << "indices:    " << view::as_array(af.indices_).transpose() << "\n";
        s << "length:     " << af.length_ << "\n";
        s << "-------------------\n";
        return s;
    }

    // static methods
    static affine_form pow(double v, affine_form const& af);
    static affine_form pow(affine_form const& af, double v) { return af.pow(v); }

private:
    std::reference_wrapper<affine_context> context_;
    double center_; // central value
    double radius_; // affine radius

    std::vector<double> deviations_; // vector of partial deviations
    std::vector<size_t> indices_; // vector of indices

    size_t length_; // the actual number of indices still in use by this affine form

    void update_radius()
    {
        radius_ = view::as_array(deviations_).segment(0, length_).abs().sum();
    }
};

} // namespace pappus

#endif
