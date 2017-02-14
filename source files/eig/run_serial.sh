#This file gives basic manual for running program to compute eigenvalues with QR decomposition.

## Compilation
g++ ./gen.c -o gen
g++ ./ser.cc -o serial -O1

## Generating file with size of matrix 20 and save it in file named "input"
./gen 20 > input

### Few example how to call program
## Print matrix A in each step & solution & time
#./serial -file input -solution

## Print solution & time
#./serial -file input -solution -silent

## Print time only - for benchmark purposes
./serial -file input -silent

## Delete files
rm gen input serial