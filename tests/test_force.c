#include <stdio.h>
#include "mdsys_t.h"
#include "force.h"

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

  sys.rx[0]=1.0;
  sys.ry[0]=0.0;
  sys.rz[0]=0.0;
  sys.vx[0]=1.0;
  sys.vy[0]=2.0;
  sys.vz[0]=2.0;

  sys.rx[1]=0.0;
  sys.ry[1]=0.0;
  sys.rz[1]=0.0;
  sys.vx[1]=1.0;
  sys.vy[1]=1.0;
  sys.vz[1]=0.0;

  /* r=1, copy-paste into wolfram alpha: -4.0*0.2379*(-12.0*pow(3.405/1,12.0)/1 +6*pow(3.405/1,6.0)/1) */

  double expected_fx=2.7726925774361321e7;
  double expected_fy=0;
  double expected_fz=0;

  sys.force=&force_Newton_3rd;
  sys.force(&sys);
  printf("force_Newton \t= %f %f %f \nExpected force \t= %f 0 0 \n",
    sys.fx[0] , sys.fy[0] , sys.fz[0] , expected_fx);

  sys.force=&force_OpenMP;
  sys.force(&sys);
  printf("force_OpenMP \t= %f %f %f \nExpected force \t= %f 0 0 \n",
    sys.fx[0] , sys.fy[0] , sys.fz[0] , expected_fx);



  sys.rx[0]=0.0;
  sys.ry[0]=-1.0;

  expected_fx=0;
  expected_fy=-2.7726925774361321e7;

  sys.force=&force_Newton_3rd;
  sys.force(&sys);
  printf("force_Newton \t= %f %f %f \nExpected force \t= 0 %f 0 \n",
    sys.fx[0] , sys.fy[0] , sys.fz[0] , expected_fy);

  sys.force=&force_OpenMP;
  sys.force(&sys);
  printf("force_OpenMP \t= %f %f %f \nExpected force \t= 0 %f 0 \n",
    sys.fx[0] , sys.fy[0] , sys.fz[0] , expected_fy);



  sys.ry[0]=0.0;
  sys.rz[0]=1.0;

  expected_fy=0;
  expected_fz=2.7726925774361321e7;

  sys.force=&force_Newton_3rd;
  sys.force(&sys);
  printf("force_Newton \t= %f %f %f \nExpected force \t= 0 0 %f \n",
    sys.fx[0] , sys.fy[0] , sys.fz[0] , expected_fz);

  sys.force=&force_OpenMP;
  sys.force(&sys);
  printf("force_OpenMP \t= %f %f %f \nExpected force \t= 0 0 %f \n",
    sys.fx[0] , sys.fy[0] , sys.fz[0] , expected_fz);

  sys.rx[0]=5.0;
  sys.ry[0]=5.0;
  sys.rz[0]=5.0;

  sys.rx[1]=-5.0;
  sys.ry[1]=-5.0;
  sys.rz[1]=-5.0;

  expected_fz=0;

  sys.force=&force_Newton_3rd;
  sys.force(&sys);
  printf("force_Newton \t= %f %f %f \nExpected force \t= 0 0 0 (distance greater than cutoff) \n",
    sys.fx[0] , sys.fy[0] , sys.fz[0]);

  sys.force=&force_Newton_3rd;
  sys.force(&sys);
  printf("force_OpenMP \t= %f %f %f \nExpected force \t= 0 0 0 (distance greater than cutoff) \n",
    sys.fx[0] , sys.fy[0] , sys.fz[0]);


  free_mdsys_mem( &sys );

  return 0;
}


// ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
//                         +6*pow(sys->sigma/r,6.0)/r);
