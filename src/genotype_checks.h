#pragma once
#include "engine_state.h"
namespace pap_internal::search {

void validate_locus_matrix(const NumericMatrix& matrix, int rows,
                                  int columns, const char* name);

int count_trio_me(const EngineState& compiled,
                         int child, int dam, int sire);

double trio_llr_ln(const EngineState& compiled, int child, int dam, int sire);
}
