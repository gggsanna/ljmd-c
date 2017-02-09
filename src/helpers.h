#ifndef __HELPERS_H__
#define __HELPERS_H__

/* helper function: zero out an array */
void azzero(double *d, const int n);

/* helper function: apply minimum image convention */
double pbc(double x, const double boxby2);

#endif
