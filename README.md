This package contains simplified MD code with multi-threading
parallelization for simulating atoms with a Lennard-Jones potential.

The bundled makefiles are set up to compile the executable once
with OpenMP disabled and once with OpenMP enabled with each build
placing the various object files in separate directories.

The examples directory contains 3 sets of example input decks
and the reference directory the corresponding outputs.

Type: make
to compile everything and: make clean
to remove all compiled objects


*************************

## How to use this repo

To compile the serial code: "make serial"
To compile the OpenMP-parallelized code: "make OMP"
To run the tests: "make test_serial" or "make test_OMP"
See the makefile for more specific test targets.

To run the executable (ljmd-serial.x or ljmd-OMP.x):

"../ljmd-ZZZ.x <n> < <input_file>" from the input files folder (ZZZ is serial or OMP, n is an integer)

## Serial optimization

Expensive operations as powers or square roots have been moved outside of loops, or replaced altogether
