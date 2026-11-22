#ifndef PAP_TOPOLOGY_TEMPLATES_H
#define PAP_TOPOLOGY_TEMPLATES_H
#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
namespace pap_internal::topology_templates {
template<class Model>struct Cache {
 struct Geometry {std::vector<int> roots,original,order;std::vector<std::array<int,2>> parents;int focal=-1,candidate=-1;};
 struct Bound {Model model;std::string fingerprint;std::vector<int> dependencies;size_t charge=0;};
 std::unordered_map<std::string,std::shared_ptr<const Geometry>> shapes;
 std::unordered_map<std::string,Bound> bound;
 size_t shape_requests=0,shape_hits=0,shape_builds=0,bound_requests=0,bound_hits=0,invalidations=0,closure_misses=0,refusals=0,bytes=0;
 double lookup_seconds=0;
 bool shape_enabled=true,bound_enabled=true;
 size_t limit=128U*1024U*1024U;
 static inline thread_local Cache* active=nullptr;
 struct Scope {Cache* old;explicit Scope(Cache* c):old(active){active=c;}~Scope(){active=old;}};
 bool room(size_t charge){if(charge>limit||bytes>limit-charge){++refusals;return false;}return true;}
 void invalidate(int individual){for(auto at=bound.begin();at!=bound.end();)if(std::binary_search(at->second.dependencies.begin(),at->second.dependencies.end(),individual)){bytes-=at->second.charge;at=bound.erase(at);++invalidations;}else ++at;}
};
template<class T>inline void append(std::string& s,const T& value){s.append(reinterpret_cast<const char*>(&value),sizeof(value));}
inline void append_string(std::string& s,const std::string& value){append(s,value.size());s+=value;}
}
#endif
