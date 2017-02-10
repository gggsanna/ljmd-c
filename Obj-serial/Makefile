# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer
LDLIBS=-lm

# list of source files for ljmd-serial
SRC = ljmd.c helpers.c force.c verlet_time_integration.c input.c output.c mdsys_t.c select_force.c
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
tests: ../test_ekin.x ../test_force.x ../test_input.x ../test_verlet.x
	../test_ekin.x 
	../test_force.x 
	../test_input.x < ../examples/argon_108.inp | diff - ../examples/argon_108.inp -y -b 
	../test_verlet.x 

../test_ekin.x:  ../tests/test_ekin.c mdsys_t.o
	$(CC) -I$(IPATH) -o $@ $(CFLAGS) $^ $(LDLIBS) 

../test_force.x: ../tests/test_force.c mdsys_t.o force.o helpers.o
	$(CC) -I$(IPATH) -o $@ $(CFLAGS) $^ $(LDLIBS) 

../test_input.x: ../tests/test_input.c mdsys_t.o input.o helpers.o
	$(CC) -I$(IPATH) -o $@ $(CFLAGS) $^ $(LDLIBS) 

../test_verlet.x: ../tests/test_verlet.c mdsys_t.o verlet_time_integration.o force.o helpers.o
	$(CC) -I$(IPATH) -o $@ $(CFLAGS) $^ $(LDLIBS) 


# compilation pattern rule for objects
%.o: $(SPATH)%.c
	$(CC) -I$(IPATH) -c $(CFLAGS) $<

