#include <stdio.h>
#include "mdsys_t.h"
#include "verlet_time_integration.h"

int main(){
  mdsys_t sys;

  sys.natoms=1;
  sys.nfi=0;
  sys.nsteps=2;
  sys.dt=2.0;
  sys.mass=10.0;
  sys.epsilon=0.2379;
  sys.sigma=3.405;
  sys.box=20.0;
  sys.rcut=10.0;

  allocate_mdsys_mem( &sys );

  sys.rx[0]=0.0;
  sys.ry[0]=0.0;
  sys.rz[0]=0.0;
  sys.vx[0]=1.0;
  sys.vy[0]=2.0;
  sys.vz[0]=3.0;

  sys.fx[0]=3.0*mvsq2e;
  sys.fy[0]=4.0*mvsq2e;
  sys.fz[0]=0.0*mvsq2e;

  propagate_position( &sys);
  /* rx = rx + dt*vx */
  double expected_rx=2.0;
  double expected_ry=4.0;
  double expected_rz=6.0;

  printf("propagate_position:\t %f %f %f \nexpected position:\t %f %f %f\n",
    sys.rx[0], sys.ry[0], sys.rz[0], expected_rx, expected_ry, expected_rz);

  propagate_velocity_half_step( &sys);
  /* vx = vx + dt/2*fx/mass/mvsq2e */
  double expected_vx=1.3;
  double expected_vy=2.4;
  double expected_vz=3.0;

  printf("propagate_velocity_half_step:\t %f %f %f \nexpected velocity:\t\t %f %f %f\n",
    sys.vx[0], sys.vy[0], sys.vz[0], expected_vx, expected_vy, expected_vz);

}
