# -*- Makefile -*-
SHELL=/bin/sh

default: serial

serial:
	$(MAKE) $(MFLAGS) -C Obj-serial ../ljmd-serial.x

OMP:
	$(MAKE) $(MFLAGS) -C Obj-OMP ../ljmd-OMP.x


test_input:
	$(MAKE) $(MFLAGS) -C Obj-serial ../test_input.x

test_force:
	$(MAKE) $(MFLAGS) -C Obj-serial ../test_force.x

test_ekin:
	$(MAKE) $(MFLAGS) -C Obj-serial ../test_ekin.x




clean:
	$(MAKE) $(MFLAGS) -C Obj-serial clean

