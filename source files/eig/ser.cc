/**
 * 
 * \brief Computing Eigen values with QR decomposition
 * Uses Housholder's QR decomposition      
 *
 * \author Ondrej Zjevik, (c) 2012
 *
 *
 *   
 * call example: mpiexec -n 1 ./program -file file_name [-silent]
*/

#include <iostream>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <sys/time.h>      
#include <math.h>
#include "cTimer.h"  

void computeQR(int argc, char* argv[], bool verbose, bool definition, float **matrixR, float **matrixQ, int rank, int lines);

//creating continuous 2d array
float** create_matrix( int numrows, int numcols, bool verbose){
	float *buffer=new float[numrows*numcols];
	float **data=new float*[numrows];
	for(int i=0;i<numrows;++i) data[i]=buffer+i*numcols;
	
	return *&data;
}

int main(int argc, char* argv[])
{ 
  bool verbose = true, definition = true, solution = false; 
  float **matrixA,**matrixR,**matrixQ;
  int rank=0,i,j,k,l,m,lines;
  FILE *file;
  arg::cTimer timer;
  
  if(rank != 0) verbose = false;
  
  if(rank == 0){
    //using extra arguments from command line
    for(i = 0; i < argc; i++){
      if(strcmp(argv[i],"-silent") == 0) verbose = false;
      if(strcmp(argv[i],"-solution") == 0) solution = true;
      if(strcmp(argv[i],"-file") == 0){
        if((file=fopen(argv[++i], "r")) == NULL) {
          printf("Cannot open file.\n");
          //MPI_Finalize();
          return 0;
        }
        else{
          definition = false;
        }
      }
    }
    
    //the input isn't complete
    if(definition){
      printf("Not enought input parameters.\nUse call like: mpiexec -n 1 ./program -file file_name [-silent]\n");
      //MPI_Finalize();
      return 0;
    }
    
    //scanning file for number of lines
    fscanf(file,"%d,\n",&lines);
    
    
    matrixQ = create_matrix(lines,lines, verbose);
    matrixR = create_matrix(lines,lines, verbose);
    matrixA = create_matrix(lines,lines, verbose);
    //loading of matrixR from file to array
    for(i = 0; i < lines; i++){
      for(j = 0; j < lines; j++){
        fscanf(file,"%f,",&matrixR[i][j]);
      }
    }
    for(i = 0; i < lines; i++){
      for(j = 0; j < lines; j++){
        matrixA[i][j] = matrixR[i][j];
      }
    }
    if(verbose){
      printf("matrix was sucesfuly loaded to memory.\nMartix:\n");
      for(i = 0; i < lines; i++){
        for(j = 0; j < lines; j++){
          if(verbose) printf("%f ",matrixR[i][j]);
        }
        if(verbose) printf("\n");
      }
    }
  }
  int a = 20;
  int b = (int)ceil(sqrt( lines ));
  //start time  
  timer.CpuStart();   
  for(i = 0; i < std::min(a,b); i++)
  {
    for(j = 0; j < lines; j++){
      for(k = 0; k < lines; k++){
        matrixR[j][k] = matrixA[j][k];
        matrixQ[j][k] = 0;
      }
    }
    computeQR(argc, argv, verbose, definition, matrixR, matrixQ, rank, lines);    
      
    //martix multiplication Q'*A*Q
    // A = Q'*R
    for(k = 0; k < lines; k++){
      for(l = 0; l < lines; l++){
        float tm = 0;
        for(m = 0; m < lines; m++){
          tm += matrixQ[m][k]*matrixA[m][l];
        }
        matrixR[k][l] = tm;
      }
    }      
        
    // R = A*Q
    for(k = 0; k < lines; k++){
      for(l = 0; l < lines; l++){
        float tm = 0;
        for(m = 0; m < lines; m++){
          tm += matrixR[k][m]*matrixQ[m][l];
        }
        matrixA[k][l] = tm;
      }
    }
    
    if(verbose){
      printf("Matrix A in %d round.\n",i);
      for(k = 0; k < lines; k++){
        for(j = 0; j < lines; j++){
          if(verbose) printf("%f ",matrixA[k][j]);
        }
        if(verbose) printf("\n");
      
      }
    }
  }
  double time;
  time = timer.CpuStop().CpuSeconds();//end time
  if(solution){
    printf("Eig values:\n--------------------------------------\n",i);
    for(k = 0; k < lines; k++){
      printf("%f\n",matrixA[k][k]);
    }
    printf("--------------------------------------\n");
  }
  //Print time
  printf("%f\n",time);
  return 0;
}

void computeQR(int argc, char* argv[], bool verbose, bool definition, float **matrixR, float **matrixQ, int rank, int lines)
{ 
  bool root = false;   
  int i,j,size;
  float **mat,**p;
  float *vec, coef;
  if(rank == 0) root = true;


  

  for(int i = 0; i < lines; i++){
    matrixQ[i][i] = 1;
  }
  if(rank == 0){    
 
  }
  
  int k,l,m;
  
  //tmp matrixR for parallel computing
  mat = create_matrix(lines,lines, false);
   
  for(i = 0; i < lines; i++){
    
    float x = 0;
    vec = new float[lines-i]; 
    for(j = i; j < lines; j++){
      vec[j-i] = -matrixR[j][i];
      x += vec[j-i]*vec[j-i];
    }
    x = sqrt(x);
     
    if(vec[0] > 0) x = -x; 
    vec[0] = vec[0] + x;
    x = 0;
    for(j = 0; j < lines-i; j++){
      x += vec[j]*vec[j];
    }
    x = sqrt(x);     
    
    if(x > 0){
      //normalizovat vec  
      for(j = 0; j < lines-i; j++){
        vec[j] /= x;
      }   
      
      //sestavit matici P
      p = create_matrix(lines-i,lines-i, false);
      for(k = 0; k < lines-i; k++){
        for(l = 0; l < lines-i; l++){
          if(k == l) p[k][k] = 1 - 2*vec[k]*vec[l];
          else p[k][l] = -2*vec[k]*vec[l];
        }
      }     
      
      //nasobeni matic (paralelizace)
      //R      
      for(k = i; k < lines; k++){
        for(l = i; l < lines; l++){
          float tm = 0;
          for(m = i; m < lines; m++){
            tm += p[k-i][m-i]*matrixR[m][l];
          }
          mat[k][l] = tm;
        }
      }          
      for(k = i; k < lines; k++){
        for(l = i; l < lines; l++){
          matrixR[k][l] = mat[k][l];
        }
      }        
      
      //Q
      for(k = 0; k < lines; k++){
        for(l = i; l < lines; l++){
          float tm = 0;
          for(m = i; m < lines; m++){
            tm += matrixQ[k][m]*p[m-i][l-i];
          }
          mat[k][l] = tm;
        }
      }
      for(k = 0; k < lines; k++){
        for(l = i; l < lines; l++){
          matrixQ[k][l] = mat[k][l];
        }
      }
    }
    
  }    
  return;                                                                                                
}
