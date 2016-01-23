// Lightweight Newick tree format parser
// Adapted from https://github.com/attractivechaos/klib

#ifndef KNHX_H_
#define KNHX_H_

#include <stddef.h>

#define KNERR_MISSING_LEFT   0x01
#define KNERR_MISSING_RGHT   0x02
#define KNERR_BRACKET        0x04
#define KNERR_COLON          0x08

typedef struct {
  int parent, n;
  int *child;
  char *name;
  double d;
} knode_t;

typedef struct {
  int error, n, max;
  knode_t *node;
} ktree_t;

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif
  
  ktree_t *kn_parse(const char *nhx);
  void kn_format(const knode_t *node, int root, kstring_t *s);

  void kn_free(ktree_t *tree);

#ifdef __cplusplus
}
#endif

#endif
