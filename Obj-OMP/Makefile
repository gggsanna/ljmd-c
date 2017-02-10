# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -fopenmp
LDLIBS=-lm

# list of source files for ljmd-serial
SRC = ljmd.c helpers.c force.c verlet_time_integration.c input.c output.c mdsys_t.c
SPATH=../src/
IPATH=../include/
OBJ=$(SRC:%.c=%.o)

default: ../ljmd-OMP.x

clean:
	rm -f *.o ../*.x 

# linker rule
../ljmd-OMP.x: $(OBJ)
	$(CC) -I$(IPATH) -o $@ $(CFLAGS) $^ $(LDLIBS)


# tests
../test_ekin-OMP.x:  ../tests/test_ekin.c mdsys_t.o
	$(CC) -I$(IPATH) -o $@ $(CFLAGS) $^ $(LDLIBS) 

../test_force-OMP.x: ../tests/test_force.c mdsys_t.o force.o helpers.o
	$(CC) -I$(IPATH) -o $@ $(CFLAGS) $^ $(LDLIBS) 

../test_verlet-OMP.x: ../tests/test_verlet.c mdsys_t.o verlet_time_integration.o force.o helpers.o
	$(CC) -I$(IPATH) -o $@ $(CFLAGS) $^ $(LDLIBS) 


# compilation pattern rule for objects
%.o: $(SPATH)%.c
	$(CC) -I$(IPATH) -c $(CFLAGS) $<

