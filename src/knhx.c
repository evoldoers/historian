// Lightweight Newick tree format parser
// Adapted from https://github.com/attractivechaos/klib

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "knhx.h"

static inline char *add_node(const char *s, ktree_t *aux, int x)
{
  char *p, *nbeg, *nend = 0;
  knode_t *r;
  if (aux->n == aux->max) {
    aux->max = aux->max? aux->max<<1 : 8;
    aux->node = (knode_t*)realloc(aux->node, sizeof(knode_t) * aux->max);
  }
  r = aux->node + (aux->n++);
  r->n = x; r->parent = -1;
  for (p = (char*)s, nbeg = p, r->d = -1.0; *p && *p != ',' && *p != ')' && *p != ';'; ++p) {
    if (*p == '[') {
      if (nend == 0) nend = p;
      do ++p; while (*p && *p != ']');
      if (*p == 0) {
	aux->error |= KNERR_BRACKET;
	break;
      }
    } else if (*p == ':') {
      if (nend == 0) nend = p;
      r->d = strtod(p + 1, &p);
      --p;
    } else if (!isgraph(*p)) if (nend == 0) nend = p;
  }
  if (nend == 0) nend = p;
  if (nend != nbeg) {
    r->name = (char*)calloc(nend - nbeg + 1, 1);
    strncpy(r->name, nbeg, nend - nbeg);
  } else r->name = strdup("");
  return p;
}

ktree_t *kn_parse(const char *nhx)
{
  char *p;
  int *stack, top, max;
  ktree_t *aux;
  knode_t *ret;

#define __push_back(y) do {				\
    if (top == max) {					\
      max = max? max<<1 : 16;				\
      stack = (int*)realloc(stack, sizeof(int) * max);	\
    }							\
    stack[top++] = (y);					\
  } while (0);

  stack = 0; top = max = 0;
  p = (char*)nhx;
  aux = (ktree_t*)calloc(1, sizeof(ktree_t));
  while (*p) {
    while (*p && !isgraph(*p)) ++p;
    if (*p == 0 || *p == ';') break;
    if (*p == ',') ++p;
    else if (*p == '(') {
      __push_back(-1);
      ++p;
    } else if (*p == ')') {
      int x = aux->n, m, i;
      for (i = top - 1; i >= 0; --i)
	if (stack[i] < 0) break;
      m = top - 1 - i;
      p = add_node(p + 1, aux, m);
      aux->node[x].child = (int*)calloc(m, sizeof(int));
      for (i = top - 1, m = m - 1; m >= 0; --m, --i) {
	aux->node[x].child[m] = stack[i];
	aux->node[stack[i]].parent = x;
      }
      top = i;
      __push_back(x);
    } else {
      __push_back(aux->n);
      p = add_node(p, aux, 0);
    }
  }
  ret = aux->node;
  free(stack);
  return aux;
}

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static inline int kputsn(const char *p, int l, kstring_t *s)
{
  if (s->l + l + 1 >= s->m) {
    s->m = s->l + l + 2;
    kroundup32(s->m);
    s->s = (char*)realloc(s->s, s->m);
  }
  memcpy(s->s + s->l, p, l);
  s->l += l; s->s[s->l] = 0;
  return l;
}

static inline int kputc(int c, kstring_t *s)
{
  if (s->l + 1 >= s->m) {
    s->m = s->l + 2;
    kroundup32(s->m);
    s->s = (char*)realloc(s->s, s->m);
  }
  s->s[s->l++] = c; s->s[s->l] = 0;
  return c;
}

static void format_node_recur(const knode_t *node, const knode_t *p, kstring_t *s, char *numbuf)
{
  if (p->n) {
    int i;
    kputc('(', s);
    for (i = 0; i < p->n; ++i) {
      if (i) kputc(',', s);
      format_node_recur(node, &node[p->child[i]], s, numbuf);
    }
    kputc(')', s);
    if (p->name) kputsn(p->name, (int)strlen(p->name), s);
    if (p->d >= 0) {
      sprintf(numbuf, ":%g", p->d);
      kputsn(numbuf, (int)strlen(numbuf), s);
    }
  } else {
    kputsn(p->name, (int)strlen(p->name), s);
    // code to show branch length for leaf branches added by IH, 6/26/2015
    if (p->d >= 0) {
      sprintf(numbuf, ":%g", p->d);
      kputsn(numbuf, (int)strlen(numbuf), s);
    }
  }
}

void kn_format(const knode_t *node, int root, kstring_t *s) // TODO: get rid of recursion
{
  char numbuf[128];
  format_node_recur(node, &node[root], s, numbuf);
  kputc(';',s);
}

void kn_free(ktree_t *tree) {
  int n;
  for (n = 0; n < tree->n; ++n)
    free (tree->node[n].name);
  free (tree->node);
  free ((void*) tree);
}

#ifdef KNHX_MAIN


int main(int argc, char *argv[])
{
  const char *s = "((a[abc],d1)x:0.5,((b[&&NHX:S=MOUSE],h2)[&&NHX:S=HUMAN:B=99][blabla][&&NHX:K=foo],c));";
  ktree_t *aux;
  int i, j;
  kstring_t str;
  aux = kn_parse(s);
  for (i = 0; i < aux->n; ++i) {
    knode_t *p = aux->node + i;
    printf("[%d] %s\t%d\t%d\t%g", i, p->name, p->parent, p->n, p->d);
    for (j = 0; j < p->n; ++j)
      printf("\t%d", p->child[j]);
    putchar('\n');
  }
  str.l = str.m = 0; str.s = 0;
  kn_format(aux->node, aux->n - 1, &str);
  puts(str.s);
  free(str.s);
  kn_free(aux);
  return 0;
}
#endif
