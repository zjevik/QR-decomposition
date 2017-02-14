#This file gives basic manual for running program to compute QR decomposition in MPI.

module add prog/pgi
## Compilation
g++ ./gen.c -o gen
pgCC -mp ./qr_omp.cc -g -o omp -O2

## Generating file with size of matrix 20 and save it in file named "input"
./gen 5 > input

### Few example how to call program
## Print matrix & solution & time
./omp -n 4 -file input

## Print time only - for benchmark purposes
#./omp -n 4 -file input -silent

## Delete files
rm gen input omp