#include <stdio.h>
#include "input.h"

int main()
{

    int nprint;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN];
    FILE *traj,*erg;
    mdsys_t sys;

    get_mdsys_stdin(&sys, restfile, trajfile, ergfile, &nprint);

    allocate_mdsys_mem(&sys);

#if 0
    printf("natoms: %i\nnfi: %i\nnsteps: %i\n", sys.natoms, sys.nfi, sys.nsteps);
    printf("dt: %f\nmass: %f\nepsilon: %f\nsigma: %f\nbox: %f\nrcut: %f\n",
           sys.dt, sys.mass, sys.epsilon, sys.sigma, sys.box, sys.rcut);
    printf("ekin: %f\nepot: %f\ntemp: %f\n", sys.ekin, sys.epot, sys.temp);
#endif

    printf("\n%i\t\t\t# natoms", sys.natoms);
    printf("\n%f\t\t# mass in AMU", sys.mass);
    printf("\n%f\t\t# epsilon in kcal/mol", sys.epsilon);
    printf("\n%f\t\t# sigma in angstrom", sys.sigma);
    printf("\n%f\t\t# rcut in angstrom", sys.rcut);
    printf("\n%f\t\t# box length (in angstrom)", sys.box);

    printf("\n%s\t\t# restart", restfile);
    printf("\n%s\t\t# trajectory", trajfile);
    printf("\n%s\t\t# energies", ergfile);

    printf("\n%i\t\t\t# nr MD steps", sys.nsteps);
    printf("\n%f\t\t# MD time step", sys.dt);
    printf("\n%i\t\t\t# output print frequency", nprint);


    printf("\n");
    printf("\n");

#if 0
    get_rest_file(restfile, &sys);
    int i;
    for (i=0; i<sys.natoms; ++i)
      printf("%e\t%e\t%e\n",*(sys.rx+i), *(sys.ry+i), *(sys.rz+i));

    for (i=0; i<sys.natoms; ++i)
      printf("%e\t%e\t%e\n",*(sys.vx+i), *(sys.vy+i), *(sys.vz+i));

    printf("\n");
#endif

    return 0;

}
