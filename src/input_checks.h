#ifndef PAP_INPUT_CHECKS_H
#define PAP_INPUT_CHECKS_H
namespace pap_internal::input_checks {
struct State {bool reference=true,share=false;};
inline thread_local State* current=nullptr;
inline bool enabled(){return !current||current->reference;}
inline bool share(){return current&&current->share;}
inline bool trusted(){return current&&!current->reference;}
struct Scope {State* old;explicit Scope(State& next):old(current){current=&next;}
 ~Scope(){current=old;}};
}
#endif
