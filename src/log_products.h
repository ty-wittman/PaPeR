#include "pap_removed_path.h"
#ifndef PAP_LOG_PRODUCTS_H
#define PAP_LOG_PRODUCTS_H
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
namespace pap_internal::log_products {
struct Value {
 std::vector<double> cells;double scale=0;
 static Value unit(size_t width){return {std::vector<double>(width,0),0};}
 void normalize(){double maximum=-INFINITY;for(double x:cells){if(std::isnan(x))throw std::invalid_argument("NaN in log product");maximum=std::max(maximum,x);}if(maximum==-INFINITY)return;for(double& x:cells)x-=maximum;scale+=maximum;}
};
inline Value multiply(const Value& a,const Value& b){
 if(a.cells.size()!=b.cells.size())throw std::invalid_argument("log product dimensions differ");
 Value out=a;out.scale+=b.scale;for(size_t i=0;i<out.cells.size();++i)out.cells[i]+=b.cells[i];out.normalize();return out;
}
// Products of disjoint ranges, including structural zeros. Neither a posterior
// nor a total product is divided by an excluded operand. Each stored object
// has exactly the original message width (G or G²), independent of degree.
class Index {
 size_t width_,count_,leaves_=1;std::vector<Value> tree_;
 public:mutable size_t work_cells=0;
 Index(size_t width,const std::vector<Value>& values):width_(width),count_(values.size()){
  while(leaves_<count_)leaves_*=2;tree_.assign(2*leaves_,Value::unit(width_));
  for(size_t i=0;i<count_;++i){if(values[i].cells.size()!=width_)throw std::invalid_argument("invalid product leaf");tree_[leaves_+i]=values[i];}
  for(size_t i=leaves_;--i;){tree_[i]=multiply(tree_[2*i],tree_[2*i+1]);work_cells+=width_;}
 }
 Value range(size_t begin,size_t end)const{
  if(begin>end||end>count_)throw std::out_of_range("invalid product range");
  Value left=Value::unit(width_),right=left;
  for(begin+=leaves_,end+=leaves_;begin<end;begin/=2,end/=2){if(begin&1){left=multiply(left,tree_[begin++]);work_cells+=width_;}if(end&1){right=multiply(tree_[--end],right);work_cells+=width_;}}
  work_cells+=width_;return multiply(left,right);
 }
 void update(size_t index,Value value){if(index>=count_||value.cells.size()!=width_)throw std::out_of_range("invalid product update");size_t at=leaves_+index;tree_[at]=std::move(value);while((at/=2)){tree_[at]=multiply(tree_[2*at],tree_[2*at+1]);work_cells+=width_;}}
};
}
#endif
