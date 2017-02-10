#include "mdsys_t.h"
#include "force.h"
#include "verlet_time_integration.h"

//extern const double mvsq2e;

/* velocity verlet */
void velverlet(mdsys_t *sys )
{
    propagate_velocity_half_step(sys);
    propagate_position(sys);
    sys->force(sys); /* update force field */
    propagate_velocity_half_step(sys);
}

/* propagate velocities by half step */
void propagate_velocity_half_step(mdsys_t *sys)
{
    /* precompute a common factor which appears in the for loop */
    double fcoeff = 0.5*sys->dt / mvsq2e / sys->mass;
    int i;
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += fcoeff * sys->fx[i];
        sys->vy[i] += fcoeff * sys->fy[i];
        sys->vz[i] += fcoeff * sys->fz[i];
    }
}

/* propagate positions by full step */
void propagate_position(mdsys_t *sys)
{
    int i;
    for (i=0; i<sys->natoms; ++i) {
        sys->rx[i] += sys->dt*sys->vx[i];
        sys->ry[i] += sys->dt*sys->vy[i];
        sys->rz[i] += sys->dt*sys->vz[i];
    }
}
