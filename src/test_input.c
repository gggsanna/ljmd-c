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




}


