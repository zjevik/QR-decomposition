#This file gives basic manual for running program to compute QR decomposition.

## Compilation
g++ ./gen.c -o gen
g++ ./qr_ser.cc -o serial -O1

## Generating file with size of matrix 20 and save it in file named "input"
./gen 5 > input

### Few example how to call program
## Print matrix & solution & time
./serial -file input

## Print time only - for benchmark purposes
#./serial -file input -silent

## Delete files
rm gen input serial