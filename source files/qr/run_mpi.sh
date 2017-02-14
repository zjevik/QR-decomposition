#This file gives basic manual for running program to compute QR decomposition in MPI.

## Compilation
g++ ./gen.c -o gen
mpic++ ./qr_mpi.cc -o mpi

## Generating file with size of matrix 20 and save it in file named "input"
./gen 5 > input

### Few example how to call program
## Print matrix & solution & time
mpiexec -n 4 ./mpi -file input

## Print time only - for benchmark purposes
#mpiexec -n 4 ./mpi -file input -silent

## Delete files
rm gen input mpi