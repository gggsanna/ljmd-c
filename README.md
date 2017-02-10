
This package contains simplified MD code with multi-threading
parallelization for simulating atoms with a Lennard-Jones or Morse potential.

The examples directory contains 2 sets of example input decks
and the reference directory the corresponding outputs.

## How to use this repository

To compile the serial code: "make serial".
To compile the OpenMP-parallelized code: "make OMP".

To run the executable (ljmd-serial.x or ljmd-OMP.x): "../ljmd-ZZZ < <input_file>".
You can provide an integer "n" as optional argument "../ljmd-ZZZ.x <n> < <input_file>": this will select a different implementation of the force function according to the following list.


 n | force
---|---
 0 | force_Old, the original function
 1 | force_OpenMP, avoid costly math operations + OpenMP
 2 | force_Newton_3rd, avoid costly math operations + only compute once per pairs
 3 | force_Morse, use Morse potential instead of Lennard-Jones

all other values default to force_Newton_3rd


To compile and run all the unit tests: "make test_serial" or "make test_OMP".
See the makefile for more specific test targets.

## Serial optimization

Expensive operations such as powers, exponentials, or square roots have been moved outside of loops, or avoided altogether when possible.
In the force_Newton_3rd implementation of the force function, the force is computed for each pair of atoms instead (half as many times as in the original code).

### Serial profiling

The profiling was done with gprof on the cluster cosilt.

Output of the original (refactored) code:
```
Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls  ms/call  ms/call  name    
 69.99    126.79   126.79     1001   126.66   181.15  force_Old
 28.22    177.90    51.11 12641018532     0.00     0.00  pbc
  1.90    181.33     3.43     3006     1.14     1.14  azzero
  0.02    181.37     0.04     1000     0.04   181.19  velverlet
  0.00    181.37     0.00     1001     0.00     0.00  ekin
  0.00    181.37     0.00      101     0.00     0.00  output
```

Output of the serial optimized code:
```
Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls  ms/call  ms/call  name    
 55.55     34.57    34.57     1001    34.54    62.23  force_Newton_3rd
 43.22     61.48    26.90 12762960210     0.00     0.00  pbc
  1.32     62.30     0.82     3006     0.27     0.27  azzero
  0.03     62.32     0.02     1000     0.02    62.25  velverlet
  0.00     62.32     0.00     1001     0.00     0.00  ekin
  0.00     62.32     0.00      101     0.00     0.00  output
```

The profiling shows that the next step should be hardcoding pbc inside the force function.

## OpenMP optimization

![time](OMPtime.png)

![scal](OMPspeedup.png)
