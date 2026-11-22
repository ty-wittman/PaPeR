#ifndef PAP_TOPOLOGY_CACHE_H
#define PAP_TOPOLOGY_CACHE_H
#include <set>
#include <string>
#include <cstddef>
#include <unordered_map>
namespace pap_internal::topology_life_reuse {
// Only successful exact component certificates are retained. Failed proofs
// run the original whole-model verifier, preserving its admission and reason.
struct Cache {bool reuse_builder=false,flat_topology=false,singleton_life=false;size_t singleton_proofs=0;size_t source_builds=0,source_copies=0;std::unordered_map<std::string,long long> witness;size_t witness_hits=0,witness_checks=0;std::set<std::string> valid;size_t requests=0,hits=0,builds=0,refusals=0,bytes=0;double lookup_seconds=0,verify_seconds=0;size_t limit=8U*1024U*1024U;};
inline thread_local Cache* active=nullptr;
struct Scope {Cache* old;explicit Scope(Cache* c):old(active){active=c;}~Scope(){active=old;}};
}
#endif
