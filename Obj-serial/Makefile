# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer
LDLIBS=-lm

# list of source files
SRC = ljmd.c helpers.c force.c verlet_time_integration.c input.c output.c mdsys_t.c
SPATH=../src/
IPATH=../include/
OBJ=$(SRC:%.c=%.o)

default: ../ljmd-serial.x

clean:
	rm -f *.mod *.o ../ljmd-serial.x

# linker rule
../ljmd-serial.x: $(OBJ)
	$(CC) -I$(IPATH) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: $(SPATH)%.c
	$(CC) -I$(IPATH) -c $(CFLAGS) $<

