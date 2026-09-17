#ifndef PAPPUS_INTERVAL_PACKING_HPP
#define PAPPUS_INTERVAL_PACKING_HPP

#include <array>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

#include <eve/module/core.hpp>

#include "pappus/interval/box.hpp"

namespace pappus {

// Deterministic Cartesian subdivision. Schedule bit i selects the half of
// schedule[i]'s dimension for leaf bit i.
template<std::floating_point T>
class subdivision_plan {
public:
    subdivision_plan(box<T> domain, std::span<std::size_t const> schedule)
        : domain_(std::move(domain)), schedule_(schedule.begin(), schedule.end()), splits_(domain_.size())
    {
        if (schedule_.size() >= sizeof(std::size_t) * 8) {
            throw std::invalid_argument("pappus::subdivision_plan: too many subdivision levels");
        }
        for (auto dimension : schedule_) {
            if (dimension >= domain_.size()) {
                throw std::invalid_argument("pappus::subdivision_plan: dimension out of range");
            }
            ++splits_[dimension];
        }
    }

    [[nodiscard]] auto dimensions() const noexcept -> std::size_t { return domain_.size(); }
    [[nodiscard]] auto leaf_count() const noexcept -> std::size_t { return std::size_t { 1 } << schedule_.size(); }

    [[nodiscard]] auto leaf(std::size_t ordinal) const -> box<T>
    {
        if (ordinal >= leaf_count()) {
            throw std::out_of_range("pappus::subdivision_plan: leaf ordinal out of range");
        }
        auto result = domain_;
        std::vector<std::size_t> cells(domain_.size());
        std::vector<std::size_t> bits(domain_.size());
        for (std::size_t bit = 0; bit < schedule_.size(); ++bit) {
            auto const dimension = schedule_[bit];
            cells[dimension] |= ((ordinal >> bit) & std::size_t { 1 }) << bits[dimension]++;
        }
        for (std::size_t dimension = 0; dimension < domain_.size(); ++dimension) {
            if (splits_[dimension] != 0) {
                result[dimension] = domain_[dimension].segment(cells[dimension], std::size_t { 1 } << splits_[dimension]);
            }
        }
        return result;
    }

private:
    box<T> domain_;
    std::vector<std::size_t> schedule_;
    std::vector<std::size_t> splits_;
};

// One aligned SIMD pack. Lane i represents exactly leaf first_leaf()+i; lanes
// after valid_lanes() are padding and must be ignored by consumers.
template<std::floating_point T, typename Wide = eve::wide<T>>
class packed_subdomains {
public:
    static constexpr std::size_t width = eve::cardinal_v<Wide>;

private:
    struct alignas(Wide) lanes {
        std::array<T, width> values;
    };

public:
    packed_subdomains(subdivision_plan<T> const& plan, std::size_t first_leaf)
        : first_leaf_(first_leaf), dimensions_(plan.dimensions()), lower_(dimensions_), upper_(dimensions_)
    {
        for (std::size_t lane = 0; lane < width && first_leaf + lane < plan.leaf_count(); ++lane) {
            auto const leaf = plan.leaf(first_leaf + lane);
            for (std::size_t dimension = 0; dimension < dimensions_; ++dimension) {
                lower_[dimension].values[lane] = leaf[dimension].inf();
                upper_[dimension].values[lane] = leaf[dimension].sup();
            }
            ++valid_lanes_;
        }
        for (std::size_t lane = valid_lanes_; lane < width; ++lane) {
            for (std::size_t dimension = 0; dimension < dimensions_; ++dimension) {
                lower_[dimension].values[lane] = T { 0 };
                upper_[dimension].values[lane] = T { 0 };
            }
        }
    }

    [[nodiscard]] auto first_leaf() const noexcept -> std::size_t { return first_leaf_; }
    [[nodiscard]] auto valid_lanes() const noexcept -> std::size_t { return valid_lanes_; }
    [[nodiscard]] auto dimensions() const noexcept -> std::size_t { return dimensions_; }
    [[nodiscard]] auto lower(std::size_t dimension) const -> Wide { return eve::load(lower_.at(dimension).values.data(), eve::as<Wide> {}); }
    [[nodiscard]] auto upper(std::size_t dimension) const -> Wide { return eve::load(upper_.at(dimension).values.data(), eve::as<Wide> {}); }
    [[nodiscard]] auto lower_data(std::size_t dimension) const -> T const* { return lower_.at(dimension).values.data(); }
    [[nodiscard]] auto upper_data(std::size_t dimension) const -> T const* { return upper_.at(dimension).values.data(); }

private:
    std::size_t first_leaf_{};
    std::size_t dimensions_{};
    std::size_t valid_lanes_{};
    std::vector<lanes> lower_;
    std::vector<lanes> upper_;
};

} // namespace pappus
#endif
