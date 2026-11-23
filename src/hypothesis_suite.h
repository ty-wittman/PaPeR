#ifndef PAP_HYPOTHESIS_SUITE_H
#define PAP_HYPOTHESIS_SUITE_H
#include <stdexcept>
#include <string>
#include <vector>
namespace pap_internal::hypothesis_suite {
enum class Suite { Legacy, S7, S37, S45 };
inline Suite parse(const std::string& value) {
  if(value=="S7")return Suite::S7;
  if(value=="S37")return Suite::S37;
  if(value=="S45")return Suite::S45;
  throw std::invalid_argument("mped_hypothesis_suite must be S7, S37 or S45");
}
inline const char* name(Suite suite) {
  switch(suite){case Suite::S7:return "S7";case Suite::S37:return "S37";
    case Suite::S45:return "S45";default:return "legacy";}
}
// Registered integer identities already encode orientation and parent slots.
// Selection precedes topology construction; no omitted models are evaluated.
inline bool contains(Suite suite,int code) {
  if(code>=0&&code<=24)return true;
  if(suite==Suite::Legacy)return code>=28&&code<=211;
  // Reserved ordinary-code gap: focal is dam/sire of the candidate.
  if(suite==Suite::S7)return code==25||code==26;
  if((code>=28&&code<=31)||(code>=36&&code<=43))return true;
  if(suite==Suite::S37&&code>=212&&code<=219)return true;
  return suite==Suite::S45&&code>=172&&code<=179;
}
inline std::vector<int> codes(Suite suite) {
  if(suite==Suite::Legacy)return {};
  std::vector<int> out;for(int code=0;code<220;++code)if(contains(suite,code))out.push_back(code);
  return out;
}
}
#endif
