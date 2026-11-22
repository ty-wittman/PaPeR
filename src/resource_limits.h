#pragma once

#include <cstddef>
#include <stdexcept>

namespace pap_internal::pedigree_likelihood {
struct ResourceLimits {
  int mib = 64, work_millions = 400;
  bool live_arena = false; // opt-in static last-use packing; arithmetic unchanged
  bool shared_piece_scratch = false;
  bool barren_pruning=true; // normalized generative leaf elimination
  bool ancestral_unit_pairs=true; // elide certified unit exterior-brood factors
  bool domain_order=true; // order is a pure function of this program domain
  bool total_work_limit=false; // reviewed aggregate program budget
  bool contraction_log_retry=true; // positive contractions cannot become exact zeros
  bool s37_refusal_repair=false; // unvalidated S37-only repair; OFF for other suites
  bool resource_planning=false;
  int resource_min_alleles=20;
  std::size_t cells() const { return std::size_t(mib)*1024U*1024U/sizeof(double); }
  double work() const { return double(work_millions)*1e6; }
  bool operator==(const ResourceLimits& x)const{return total_work_limit==x.total_work_limit&&contraction_log_retry==x.contraction_log_retry&&domain_order==x.domain_order&&s37_refusal_repair==x.s37_refusal_repair&&mib==x.mib&&work_millions==x.work_millions&&live_arena==x.live_arena&&shared_piece_scratch==x.shared_piece_scratch&&resource_planning==x.resource_planning&&resource_min_alleles==x.resource_min_alleles;}
  void validate()const {
    if(resource_min_alleles<2||resource_min_alleles>20)throw std::invalid_argument("mped2_resource_min_alleles must be between 2 and 20");
    if(mib!=64&&mib!=128&&mib!=256&&mib!=512)throw std::invalid_argument("mped2_memory_mib must be 64, 128, 256 or 512");
    if(work_millions!=400&&work_millions!=800&&work_millions!=1600&&work_millions!=2000&&work_millions!=5000)throw std::invalid_argument("mped2_work_millions must be 400, 800, 1600, 2000 or 5000");
  }
  const char* table_reason()const{return mib==64?"individual scalar table exceeds 64 MiB":mib==128?"individual scalar table exceeds 128 MiB":mib==256?"individual scalar table exceeds 256 MiB":"individual scalar table exceeds 512 MiB";}
  const char* arena_reason()const{return mib==64?"aggregate scalar arena exceeds 64 MiB":mib==128?"aggregate scalar arena exceeds 128 MiB":mib==256?"aggregate scalar arena exceeds 256 MiB":"aggregate scalar arena exceeds 512 MiB";}
  const char* work_reason()const{return work_millions==400?"elimination work exceeds 400M operations":work_millions==800?"elimination work exceeds 800M operations":work_millions==1600?"elimination work exceeds 1600M operations":work_millions==2000?"elimination work exceeds 2000M operations":"elimination work exceeds 5000M operations";}
};
}  // namespace pap_internal::pedigree_likelihood
