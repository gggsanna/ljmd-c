#include "select_force.h"
#include "force.h"
#include "mdsys_t.h"

/* force_id | force
 * 0 : force_Old, the original function
 * 1 : force_OpenMP, avoid costly math operations + OpenMP
 * 2 : force_Newton_3rd, avoid costly math operations + only compute once per pairs
 * 3 : force_Morse, use Morse potential instead of Lennard-Jones
 * other values all default to force_Newton_3rd
 */

void select_force( mdsys_t * sys, int force_id){
#ifdef _OPENMP
  sys->force=&force_OpenMP;
  return;
#endif
  if(force_id==0)
    sys->force=&force_Old;
  else if(force_id==1)
    sys->force=&force_OpenMP;
  else if(force_id==2)
    sys->force=&force_Newton_3rd;
  else if(force_id==3)
    sys->force=&force_Morse;
  else
    sys->force=&force_Newton_3rd;
  return;
}
