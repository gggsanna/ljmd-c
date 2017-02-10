#ifndef __MDSYS_T_H__
#define __MDSYS_T_H__

#include "param.h"

/* structure to hold the complete information
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
    void (*force)(struct _mdsys *);
};
typedef struct _mdsys mdsys_t;

/* memory management functions */
void allocate_mdsys_mem(mdsys_t *sys);
void free_mdsys_mem(mdsys_t *sys);

/* physical quantity computation */
void ekin(mdsys_t *sys);

#endif
