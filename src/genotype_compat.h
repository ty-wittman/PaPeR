#ifndef PAP_GENOTYPE_COMPAT_H
#define PAP_GENOTYPE_COMPAT_H

inline bool typed_diploid(int a, int b) {
  return a > 0 && b > 0;
}

inline bool diploid_disjoint(int a, int b, int c, int d) {
  return typed_diploid(a, b) && typed_diploid(c, d) &&
    a != c && a != d && b != c && b != d;
}

#endif
