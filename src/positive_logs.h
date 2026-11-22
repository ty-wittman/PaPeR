#ifndef PAP_POSITIVE_LOGS_H
#define PAP_POSITIVE_LOGS_H
#include <cmath>
#include <stdexcept>
#include <limits>
#include <algorithm>
namespace pap_internal::pedigree_likelihood {
struct StorageRange {};
inline double log_add(double a,double b){
  if(a==-INFINITY)return b;if(b==-INFINITY)return a;
  const double hi=std::max(a,b);return hi+std::log1p(std::exp(std::min(a,b)-hi));
}
// Used only after the ordinary representation cannot preserve its support.
struct LogPositive {
  double log;
  LogPositive():log(-INFINITY){}
  LogPositive(double x):log(x==0.?-INFINITY:x==1.?0.:x==.5?-0.6931471805599453094:std::log(x)){}
  static LogPositive from_log(double x){LogPositive v;v.log=x;return v;}
  LogPositive& operator+=(LogPositive b){log=log_add(log,b.log);return *this;}
  LogPositive& operator*=(LogPositive b){log+=b.log;return *this;}
  friend LogPositive operator+(LogPositive a,LogPositive b){return a+=b;}
  friend LogPositive operator*(LogPositive a,LogPositive b){return a*=b;}
  friend LogPositive operator/(LogPositive a,LogPositive b){return from_log(a.log-b.log);}
  friend bool operator<(LogPositive a,LogPositive b){return a.log<b.log;}
  friend bool operator>(LogPositive a,LogPositive b){return a.log>b.log;}
  friend bool operator==(LogPositive a,LogPositive b){return a.log==b.log;}
};
static_assert(sizeof(LogPositive)==sizeof(double));
inline double positive_log(long double x){return std::log(x);}
inline double positive_log(LogPositive x){return x.log;}
inline double positive_log(double x){return std::log(x);}
inline double storage_product(double a,double b){double v=a*b;if(a>0&&b>0&&v<std::numeric_limits<double>::min())throw StorageRange();return v;}
inline LogPositive storage_product(LogPositive a,LogPositive b){return a*b;}
}
#endif
