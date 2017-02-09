#include "mdsys_t.h"
#include <stdlib.h>

extern const double mvsq2e;
extern const double kboltz;

/* compute kinetic energy */
void ekin(mdsys_t *sys)
{
    int i;

    sys->ekin=0.0;
    for (i=0; i<sys->natoms; ++i) {
        sys->ekin += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}


void allocate_mdsys_mem(mdsys_t *sys){

    /* allocate memory */
    sys->rx=(double *)malloc(sys->natoms*sizeof(double));
    sys->ry=(double *)malloc(sys->natoms*sizeof(double));
    sys->rz=(double *)malloc(sys->natoms*sizeof(double));
    sys->vx=(double *)malloc(sys->natoms*sizeof(double));
    sys->vy=(double *)malloc(sys->natoms*sizeof(double));
    sys->vz=(double *)malloc(sys->natoms*sizeof(double));
    sys->fx=(double *)malloc(sys->natoms*sizeof(double));
    sys->fy=(double *)malloc(sys->natoms*sizeof(double));
    sys->fz=(double *)malloc(sys->natoms*sizeof(double));

}


void free_mdsys_mem(mdsys_t *sys){

    free(sys->rx);
    free(sys->ry);
    free(sys->rz);
    free(sys->vx);
    free(sys->vy);
    free(sys->vz);
    free(sys->fx);
    free(sys->fy);
    free(sys->fz);

}

