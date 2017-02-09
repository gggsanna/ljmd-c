#include "mdsys_t.h"
#include "force.h"
#include "verlet_time_integration.h"

//extern const double mvsq2e;

/* velocity verlet */
void velverlet(mdsys_t *sys)
{
    propagate_velocity_half_step(sys);
    propagate_position(sys);
    force(sys); /* update force field */
    propagate_velocity_half_step(sys);
}

/* propagate velocities by half step */
void propagate_velocity_half_step(mdsys_t *sys)
{
    int i;
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass;
        sys->vy[i] += 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass;
        sys->vz[i] += 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass;
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
