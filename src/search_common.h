#pragma once
#include "pap_removed_path.h"
#include "input_checks.h"
#include "typed_age.h"
#include <atomic>
// Shared dependencies of the typed-parent search implementation.

#include <Rcpp.h>
#include <array>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <deque>
#include <functional>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <limits>
#include <list>
#include <mutex>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <stdexcept>
#include <condition_variable>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#if defined(__APPLE__) || defined(__linux__)
#include <sys/resource.h>
#endif
#ifdef __APPLE__
#include <mach/mach.h>
#endif
#include "latent_mendelian_support.h"
#include "pedigree_graph_snapshot.h"
#include "probability_log.h"
#include "partition_totals.h"
#include "score_types.h"
#include "predictive_models.h"
#include "panel_rejection.h"
#include "panel_staging.h"
#include "support_rejection.h"
#include "joint_support.h"
#include "initial_pair_check.h"
#include "predictive_error_bounds.h"
#include "result_storage.h"
#include "scorer.h"

using namespace Rcpp;
using pap_internal::DirtyRegion;
using pap_internal::PedigreeGraphSnapshot;
using pap_internal::kNoEvidenceRow;
using pap_internal::union_nodes;
