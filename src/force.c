#include <math.h>

#include "helpers.h" /* contains azzero, pbc */
#include "force.h"

/* IF NOT USING OPENMP, only loop over pairs and update
 * the force of two particles at each iteration */
void force_Newton_3rd(mdsys_t *sys)
{
    double rsq,rinv,r6,ffac;
    double rx,ry,rz;
    int i,j;

    /* pre-compute expensive operations whose result is used inside the loop */
    double c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
    double c6=4.0*sys->epsilon*pow(sys->sigma, 6.0);
    double rcsq = sys->rcut * sys->rcut;
    double half_box = 0.5*sys->box;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    for(i=0; i < (sys->natoms); ++i) {
        for(j=i+1; j < (sys->natoms); ++j) {
            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], half_box);
            ry=pbc(sys->ry[i] - sys->ry[j], half_box);
            rz=pbc(sys->rz[i] - sys->rz[j], half_box);

            /* compute force and energy if within cutoff */
            rsq = rx*rx + ry*ry + rz*rz;
            if (rsq < rcsq) {
                rinv=1.0/rsq;
                r6=rinv*rinv*rinv;

                ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                sys->epot += r6*(c12*r6 - c6);

                sys->fx[i] += rx*ffac;
                sys->fy[i] += ry*ffac;
                sys->fz[i] += rz*ffac;

                sys->fx[j] += -rx*ffac;
                sys->fy[j] += -ry*ffac;
                sys->fz[j] += -rz*ffac;
            }
        }
    }
}


void force_OpenMP(mdsys_t *sys)
{
    double rsq,rinv,r6,ffac;
    double rx,ry,rz;
    int i,j;

    /* pre-compute expensive operations whose result is used inside the loop */
    double c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
    double c6=4.0*sys->epsilon*pow(sys->sigma, 6.0);
    double rcsq = sys->rcut * sys->rcut;
    double half_box = 0.5*sys->box;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    for(i=0; i < (sys->natoms); ++i) {
        for(j=0; j < (sys->natoms); ++j) {
            /* particles have no interactions with themselves */
            if (i==j) continue;

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], half_box);
            ry=pbc(sys->ry[i] - sys->ry[j], half_box);
            rz=pbc(sys->rz[i] - sys->rz[j], half_box);

            /* compute force and energy if within cutoff */
            rsq = rx*rx + ry*ry + rz*rz;
            if (rsq < rcsq) {
                rinv=1.0/rsq;
                r6=rinv*rinv*rinv;

                ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                sys->epot += r6*(c12*r6 - c6);

                sys->fx[i] += rx*ffac;
                sys->fy[i] += ry*ffac;
                sys->fz[i] += rz*ffac;
            }
        }
    }
}

/* D_e*( exp(-2a(r-r_e)) - 2*exp(a(r - r_e)) )
 D_e is epsilon in lennard-jones
 a is 0.15
 r_e is 1.122462*sigma = 44.84
 exp_r is exp(-a(r-r_e))
 V(r)
 */
void force_Morse( mdsys_t *sys){
  double exp_r, r_e, alpha, ffac;
  double rx,ry,rz,r;
  int i,j;

  /* pre-compute operations whose result is used inside the loop */
  alpha = 0.15;
  r_e = 1.122462*sys->sigma;
  double half_box = 0.5*sys->box;
  double twice_a_D = 2*alpha*sys->epsilon;


  /* zero energy and forces */
  sys->epot=0.0;
  azzero(sys->fx,sys->natoms);
  azzero(sys->fy,sys->natoms);
  azzero(sys->fz,sys->natoms);

  for(i=0; i < (sys->natoms); ++i) {
      for(j=0; j < (sys->natoms); ++j) {
          /* particles have no interactions with themselves */
          if (i==j) continue;

          /* get distance between particle i and j */
          rx=pbc(sys->rx[i] - sys->rx[j], half_box);
          ry=pbc(sys->ry[i] - sys->ry[j], half_box);
          rz=pbc(sys->rz[i] - sys->rz[j], half_box);
          r = sqrt(rx*rx + ry*ry + rz*rz);

          /* compute force and energy if within cutoff */
          if (r < sys->rcut) {
              exp_r = exp(-alpha*(r - r_e));

              ffac = twice_a_D*exp_r*(1 - exp_r);
              sys->epot += sys->epsilon*exp_r*(1 - 2*exp_r);

              sys->fx[i] += rx*ffac;
              sys->fy[i] += ry*ffac;
              sys->fz[i] += rz*ffac;
          }
      }
  }
}
