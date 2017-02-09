#include <stdio.h>
#include "mdsys_t.h"

int main(){
  mdsys_t sys;

  sys.natoms=2;
  sys.nfi=0;
  sys.nsteps=2;
  sys.dt=10.0;
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
  sys.vz[0]=2.0;

  sys.rx[1]=1.0;
  sys.ry[1]=1.0;
  sys.rz[1]=1.0;
  sys.vx[1]=3.0;
  sys.vy[1]=4.0;
  sys.vz[1]=0.0;

  double expected_ekin = 170.0*mvsq2e; /* 1/2*mass*mvsq2e*(v[0]^2 + v[1]^2) */

  ekin( &sys );
  printf("Ekin = %f . Expected ekin = %f\n",sys.ekin, expected_ekin);

  free_mdsys_mem( &sys );

  return 0;
}
