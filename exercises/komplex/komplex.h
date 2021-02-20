#ifndef HAVE_KOMPLEX_H
#define HAVE_KOMPLEX_H

struct komplex {double re; double im;};
typedef struct komplex komplex;

void komplex_print (char* s, komplex z);
void komplex_set (komplex* z, double x, double y);
komplex komplex_new (double x, double y);
komplex komplex_add (komplex a, komplex b);
komplex komplex_sub (komplex a, komplex b);

#endif

