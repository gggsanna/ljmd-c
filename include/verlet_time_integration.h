#ifndef __VERLET_TIME_INTEGRATION_H__
#define __VERLET_TIME_INTEGRATION_H__

#include "mdsys_t.h"

void velverlet(mdsys_t *sys, void (*force)(mdsys_t *));
void propagate_velocity_half_step(mdsys_t *sys);
void propagate_position(mdsys_t *sys);

#endif
