#ifndef PAP_JOINT_SUPPORT_H
#define PAP_JOINT_SUPPORT_H
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace pap_internal::joint_support_gate {
struct Result {
  std::array<double,3> support{{-INFINITY,-INFINITY,-INFINITY}};
  int required=0;
  bool available=true, pass=true, rejected=false;
};
inline double contrast(double a,double b) {
  if(std::isnan(a)||std::isnan(b)||(a==b&&!std::isfinite(a)))
    return std::numeric_limits<double>::quiet_NaN();
  return a-b;
}
template<class Score>
Result evaluate(int dam,int sire,bool incumbent_dam,bool incumbent_sire,
                double threshold,Score score,bool available_baselines=true) {
  Result out;
  if(!dam&&!sire)return out;
  if(dam&&!incumbent_dam)out.required|=1;
  if(sire&&!incumbent_sire)out.required|=2;
  if(dam&&sire&&!(incumbent_dam||incumbent_sire))out.required|=4;
  const double joint=score(dam,sire);
  if(dam&&sire) {
    const double neither=score(0,0),dam_only=score(dam,0),sire_only=score(0,sire);
    out.support[0]=contrast(joint,sire_only);
    out.support[1]=contrast(joint,dam_only);
    const double best=std::isnan(neither)||std::isnan(dam_only)||std::isnan(sire_only)?
      std::numeric_limits<double>::quiet_NaN():std::max(neither,std::max(dam_only,sire_only));
    out.support[2]=contrast(joint,best);
    if(available_baselines&&std::isnan(out.support[2])) {
      // A missing baseline cannot hide an available negative certificate.
      // A positive result still requires every baseline.
      for(double baseline:{neither,dam_only,sire_only}) {
        double x=contrast(joint,baseline);
        if(!std::isnan(x)&&x<=threshold)
          out.support[2]=std::isnan(out.support[2])?x:std::min(out.support[2],x);
      }
    }
  } else out.support[dam?0:1]=contrast(joint,score(0,0));
  for(int j=0;j<3;++j)if(out.required&(1<<j)) {
    const double x=out.support[j];
    if(std::isnan(x)){out.available=false;out.pass=false;}
    else if(x<=threshold){out.rejected=true;out.pass=false;}
  }
  // A known negative is a rejection certificate even if another required
  // comparison is unavailable. An unavailable comparison alone never rejects.
  return out;
}
}
#endif
