#ifndef PAP_PROGRAM_REGISTRY_H
#define PAP_PROGRAM_REGISTRY_H
// Included after Program/resource-plan definitions, inside pedigree_likelihood namespace.
// Charges vector capacities and conservative allocator/node overhead. Shared
// descendants are charged once per entry; sharing between entries is charged
// again, so the retention certificate never relies on external ownership.
template<class T> inline std::size_t retained_vector(const std::vector<T>& x){return x.capacity()*sizeof(T)+(x.capacity()?64:0);}
inline std::size_t program_storage(const Program& p,std::set<const Program*>& seen){
 if(!seen.insert(&p).second)return 0;
 std::size_t n=sizeof(Program)+64+retained_vector(p.audit_domains)+retained_vector(p.output_dims)+retained_vector(p.slots)+retained_vector(p.steps)+retained_vector(p.scalars)+retained_vector(p.unary_offset);
 n+=retained_vector(p.audit_scopes);for(const auto& row:p.audit_scopes)n+=retained_vector(row);
 for(const auto& s:p.steps)n+=retained_vector(s.inputs)+retained_vector(s.dims)+retained_vector(s.axis0)+retained_vector(s.strides);
 if(p.exact){const auto& e=*p.exact;n+=sizeof(e)+64+retained_vector(e.domains)+retained_vector(e.original_domains)+retained_vector(e.restricted)+retained_vector(e.factors)+retained_vector(e.steps)+retained_vector(e.scalars);
  for(const auto& f:e.factors)n+=retained_vector(f.scope);
  for(const auto& s:e.steps)n+=retained_vector(s.inputs)+retained_vector(s.axes)+retained_vector(s.scatter)+retained_vector(s.axis0)+retained_vector(s.axis1)+retained_vector(s.strides)+retained_vector(s.out_strides);
 }
 return n;
}
inline std::size_t program_storage(const Program& p){std::set<const Program*> seen;return program_storage(p,seen);}

class ProgramRegistry {
 using Key=std::vector<int>;
 using Order=std::list<const Key*>;
 struct Entry {std::shared_ptr<const Program> value;std::size_t bytes;Order::iterator recent;};
 using Table=std::unordered_map<Key,Entry,ShapeHash>;
 Table table;Order order;std::size_t entries_bytes=0;
 void release(){order.clear();Table empty;table.swap(empty);entries_bytes=0;}
 void evict(){const auto* key=order.back();auto at=table.find(*key);entries_bytes-=at->second.bytes;order.pop_back();table.erase(at);++evictions;}
 public:
 bool bounded=true;std::size_t limit=16U*1024U*1024U;
 std::size_t hits=0,misses=0,evictions=0,declines=0,insertions=0,clearings=0,peak_bytes=0;
 std::size_t bytes()const{return entries_bytes+(table.bucket_count()>1?table.bucket_count()*sizeof(void*):0);}
 std::size_t size()const{return table.size();}
 bool empty()const{return table.empty();}
 void reset_counters(){hits=misses=evictions=declines=insertions=clearings=0;peak_bytes=bytes();}
 void clear(){release();reset_counters();}
 void configure(bool on,std::size_t cap){clear();bounded=on;limit=cap;}
 std::shared_ptr<const Program> lookup(const Key& key){
  auto at=table.find(key);if(at==table.end()){++misses;return {};}
  ++hits;if(bounded)order.splice(order.begin(),order,at->second.recent);return at->second.value;
 }
 void insert(Key key,std::shared_ptr<const Program> value){
  const std::size_t charge=bounded?sizeof(Entry)+sizeof(Key)+128+key.capacity()*sizeof(int)+program_storage(*value):0;
  if(bounded){
   // A cache miss is never a numerical refusal: the caller already owns and
   // executes this immutable program even when retaining it is uneconomical.
   if(charge>limit){++declines;return;}
   while(!order.empty()&&(bytes()>limit-charge||table.size()>=50000))evict();
  }else if(table.size()>=50000){release();++clearings;}
  auto [at,fresh]=table.emplace(std::move(key),Entry{std::move(value),charge,{}});
  if(!fresh)return;
  if(bounded){order.push_front(&at->first);at->second.recent=order.begin();}entries_bytes+=charge;++insertions;
  if(bounded)while(!order.empty()&&bytes()>limit)evict();
  peak_bytes=std::max(peak_bytes,bytes());
 }
};
inline ProgramRegistry& program_registry(){thread_local ProgramRegistry registry;return registry;}
inline std::size_t& program_registry_reused_entries(){thread_local std::size_t n=0;return n;}
// One bounded registry per native caller thread. A fresh public-call generation
// clears it; calibration restarts carrying the same generation retain programs.
// Keys own every symbolic domain/resource input, and never contain rates.
inline void configure_program_registry(bool bounded,std::size_t cap,bool reuse,double generation){
 thread_local double previous=0;
 auto& registry=program_registry();
 const bool valid=reuse&&std::isfinite(generation)&&generation>0&&generation==std::floor(generation);
 const bool keep=valid&&generation==previous&&registry.bounded==bounded&&registry.limit==cap;
 program_registry_reused_entries()=keep?registry.size():0;
 if(keep)registry.reset_counters();else registry.configure(bounded,cap);
 previous=valid?generation:0;
}
#endif
