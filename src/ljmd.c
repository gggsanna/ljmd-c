/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * refactored c version.
 */

#include <stdio.h>

#include "mdsys_t.h"
#include "helpers.h"
#include "force.h"
#include "verlet_time_integration.h"
#include "input.h"
#include "output.h"

/* generic file- or pathname buffer length */
#define BLEN 200

/* main */
int main(int argc, char **argv)
{
    int nprint;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN];
    FILE *traj,*erg;
    mdsys_t sys;

#ifdef _OPENMP
    sys.force = &force_OpenMP;
#else
    sys.force = &force_Newton_3rd;
#endif

    if( get_mdsys_stdin(&sys, restfile, trajfile, ergfile, &nprint) != 0 )
    {
      perror("error in stdin");
      return 1;
    }

    allocate_mdsys_mem(&sys);

    if( get_rest_file(restfile, &sys) != 0 )
      return 3;

    /* initialize forces and energies.*/
    sys.nfi=0;
    sys.force(&sys);
    ekin(&sys);

    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);

    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);

        /* propagate system and recompute energies */
        velverlet(&sys);
        ekin(&sys);
    }
    /**************************************************/

    /* clean up: close files, free memory */
    printf("Simulation Done.\n");
    fclose(erg);
    fclose(traj);

    free_mdsys_mem(&sys);

    return 0;
}
