#ifndef PAP_OPTIONAL_WORK_H
#define PAP_OPTIONAL_WORK_H
#include <array>
#include <limits>
#include <set>
namespace pap_internal::optional_demand {
// Availability owns scientific keys, not shared numerical-entry pointers.
// Keep existing whole-computed-panel counting once a key is required.
struct Ledger {
  using Key=std::array<int,5>;
  std::set<Key> extra;
  bool prescore=false,observation=false;
  void request(const Key& key,bool already_present=true) {
    if(observation)return;
    if(prescore){if(!already_present)extra.insert(key);}
    else if(!extra.empty())extra.erase(key);
  }
  void identity(int dam,int sire) {
    if(prescore||observation||extra.empty())return;
    auto at=extra.lower_bound({{5,dam,sire,std::numeric_limits<int>::min(),std::numeric_limits<int>::min()}});
    while(at!=extra.end()&&(*at)[0]==5&&(*at)[1]==dam&&(*at)[2]==sire)
      at=extra.erase(at);
  }
  bool optional(const Key& key)const{return !extra.empty()&&extra.count(key);}
  void clear(){extra.clear();prescore=observation=false;}
};
struct Scope {
  bool& value;bool prior;
  explicit Scope(bool& flag):value(flag),prior(flag){value=true;}
  ~Scope(){value=prior;}
  Scope(const Scope&)=delete;Scope& operator=(const Scope&)=delete;
};
}
#endif
