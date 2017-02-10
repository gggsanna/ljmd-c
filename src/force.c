#include <math.h>

#include "helpers.h" /* contains azzero, pbc */
#include "force.h"

void force(mdsys_t *sys)
{
    double rsq,rinv,r6,ffac;
    double rx,ry,rz;
    int i,j;

    /* pre-compute expensive operations whose result is used inside the loop */
    double c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
    double c6=4.0*sys->epsilon*pow(sys->sigma, 6.0);
    double rcsq = sys->rcut * sys->rcut;
    double half_box = 0.5*sys->box;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    for(i=0; i < (sys->natoms); ++i) {
/* if not using openmp, only loop over pairs and update
 * the force of two particles at each iteration */
#ifndef _OPENMP
        for(j=i+1; j < (sys->natoms); ++j) {
#else
        for(j=0; j < (sys->natoms); ++j) {
            /* particles have no interactions with themselves */
            if (i==j) continue;
#endif

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], half_box);
            ry=pbc(sys->ry[i] - sys->ry[j], half_box);
            rz=pbc(sys->rz[i] - sys->rz[j], half_box);

            /* compute force and energy if within cutoff */
            rsq = rx*rx + ry*ry + rz*rz;
            if (rsq < rcsq) {
                rinv=1.0/rsq;
                r6=rinv*rinv*rinv;

                ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                sys->epot += r6*(c12*r6 - c6);

                sys->fx[i] += rx*ffac;
                sys->fy[i] += ry*ffac;
                sys->fz[i] += rz*ffac;

#ifndef _OPENMP
                sys->fx[j] += -rx*ffac;
                sys->fy[j] += -ry*ffac;
                sys->fz[j] += -rz*ffac;
#endif
            }
        }
    }
}
