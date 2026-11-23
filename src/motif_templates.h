#ifndef PAP_MOTIF_TEMPLATES_H
#define PAP_MOTIF_TEMPLATES_H
#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
namespace pap_internal::motif_templates {
// Immutable literal-role instructions only. Actual pedigree identities, slots,
// sharing, life constraints, normalization and admissions are never cached here.
struct Operation {int child,parent,side;};
struct Distinct {int first,second;std::string reason;};
struct Program {int a=0,b=0;std::size_t roles=0;std::vector<Operation> operations;std::vector<Distinct> distinct;std::string pre_failure,post_failure;};
struct Cache {std::unordered_map<std::string,std::shared_ptr<const Program>> entries;std::size_t requests=0,hits=0,builds=0,refusals=0,bytes=0,limit=1024U*1024U;};
inline thread_local Cache* active=nullptr;
struct Scope {Cache* previous;explicit Scope(Cache* next):previous(active){active=next;}~Scope(){active=previous;}};
}
#endif
