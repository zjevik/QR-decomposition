/**
 * 
 * \brief Generating a symmetric matrix for qr_mpi/omp/ser
 *
 * Generate an input for qr_mpi, qr_omp or qr_ser scripts
 *
 * \author Ondrej Zjevik, (c) 2012
 *
 *
 * call example: ./gen 4
 * 
 * output:
 * 
 *   4 3 2 1
 *   3 4 3 2
 *   2 3 4 3
 *   1 2 3 4
 *   
*/

#include <iostream>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>     
#include <math.h>

                    
int main(int argc, char* argv[])
{
  int lines = atoi(argv[1]);
  printf("%d,\n",lines); 
  for(int i = 0; i < lines; i++){
    for(int j = 0; j < lines; j++){
      printf("%d,",lines-abs(j-i));
    }
    printf("\n");
  }
  printf("\n");
}
