# -*- Makefile -*-
SHELL=/bin/sh

default: serial

serial:
	$(MAKE) $(MFLAGS) -C Obj-$@

clean:
	$(MAKE) $(MFLAGS) -C Obj-serial clean

