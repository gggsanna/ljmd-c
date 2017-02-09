# -*- Makefile -*-
SHELL=/bin/sh

default: serial

serial:
	$(MAKE) $(MFLAGS) -C Obj-serial ../ljmd-serial.x

test_input:
	$(MAKE) $(MFLAGS) -C Obj-serial ../test_input.x


clean:
	$(MAKE) $(MFLAGS) -C Obj-serial clean

