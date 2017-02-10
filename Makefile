# -*- Makefile -*-
SHELL=/bin/sh

default: serial


serial:
	$(MAKE) $(MFLAGS) -C Obj-serial ../ljmd-serial.x

OMP:
	$(MAKE) $(MFLAGS) -C Obj-OMP ../ljmd-OMP.x


test_serial:
	$(MAKE) $(MFLAGS) -C Obj-serial tests

test_OMP:
	$(MAKE) $(MFLAGS) -C Obj-OMP tests


clean:
	$(MAKE) $(MFLAGS) -C Obj-serial clean
	$(MAKE) $(MFLAGS) -C Obj-OMP clean

