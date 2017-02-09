# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer
LDLIBS=-lm

# list of source files for ljmd-serial
SRC = ljmd.c helpers.c force.c verlet_time_integration.c input.c output.c mdsys_t.c
SPATH=../src/
IPATH=../include/
OBJ=$(SRC:%.c=%.o)

default: ../ljmd-serial.x

clean:
	rm -f *.o ../*.x 

# linker rule
../ljmd-serial.x: $(OBJ)
	$(CC) -I$(IPATH) -o $@ $(CFLAGS) $^ $(LDLIBS)


# tests
../test_input.x:  mdsys_t.o input.o helpers.o
	$(CC) -I$(IPATH) -o $@ $(CFLAGS) $^ $(LDLIBS) ../tests/test_input.c

../test_ekin.x:  mdsys_t.o
	$(CC) -I$(IPATH) -o $@ $(CFLAGS) $^ $(LDLIBS) ../tests/test_ekin.c

../test_force.x:  mdsys_t.o force.o helpers.o
	$(CC) -I$(IPATH) -o $@ $(CFLAGS) $^ $(LDLIBS) ../tests/test_force.c


# compilation pattern rule for objects
%.o: $(SPATH)%.c
	$(CC) -I$(IPATH) -c $(CFLAGS) $<

