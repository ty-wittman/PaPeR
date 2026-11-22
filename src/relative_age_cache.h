#ifndef PAP_RELATIVE_AGE_CACHE_H
#define PAP_RELATIVE_AGE_CACHE_H
#include <cstddef>
#include <array>
#include <string>
#include <unordered_map>
#include <vector>
namespace pap_internal::relative_life {
// Numerical feasibility is invariant to translation of every dated constraint.
// Negative diagnostics additionally retain the absolute anchor: Bellman-Ford's
// zero initialization can select a different witness after translation.
struct Constraint {int kind=0,node=-1,parent=-1,side=-1;};
struct Outcome {bool complete=false,feasible=true;std::vector<Constraint> cycle;};
struct Template {std::vector<int> roles;std::unordered_map<std::string,Outcome> outcomes;};
struct Cache {
 std::unordered_map<std::string,Template> shapes;
 std::unordered_map<std::string,Outcome> negative;
 // Current accepted typed pedigree, one-based IDs. Updated on every move.
 std::vector<std::array<int,2>> source_parents;
 std::size_t bytes=0;
 std::size_t limit=64U*1024U*1024U;
};
inline thread_local Cache* current=nullptr;
struct Scope {Cache* previous;explicit Scope(Cache* cache):previous(current){current=cache;}~Scope(){current=previous;}};
}
#endif
