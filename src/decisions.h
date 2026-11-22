#pragma once
#include <vector>
namespace pap_internal::search {
struct SearchInputs;
struct ConfigurationScore;
struct SideReduction;

void compute_diagnostics(
    const SearchInputs& input, int current_index,
    std::vector<ConfigurationScore>& configurations,
    const std::vector<SideReduction>* precomputed=nullptr,
    std::vector<SideReduction>* numerical_out=nullptr);


int choose_configuration(
    const SearchInputs& input, int current_index,
    const std::vector<ConfigurationScore>& configurations,
    bool& near_tie, int& rejection_reason);


}
