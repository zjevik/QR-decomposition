/**
 * 
 * \brief Housholder's QR decomposition for MPI
 *
 *
 * \author Ondrej Zjevik, (c) 2012
 *
 *
 *   
 * call example: mpiexec -n 1 ./program -file file_name [-silent]
*/

#include "mpi.h"
#include <iostream>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include "cTimer.h"      

                                                                    
#define MAX_PROCESSES 200

//creating continuous 2d array
float** create_matrix( int numrows, int numcols, bool verbose){
	float *buffer=new float[numrows*numcols];
	float **data=new float*[numrows];
	for(int i=0;i<numrows;++i) data[i]=buffer+i*numcols;
	
	return *&data;
}

int main(int argc, char* argv[])
{ 
  bool verbose = true, definition = true, root = false;   
  int i,j,rank=0,size,lines;
  FILE *file;
  float **matrixA,**matrixQ,**mat,**p,**matTmp,**matTmp2;
  float *vec, coef;
  struct timeval atime,ztime;
  int displs[MAX_PROCESSES], send_counts[MAX_PROCESSES], displs2[MAX_PROCESSES], send_counts2[MAX_PROCESSES];
  arg::cTimer timer;
  
  MPI_Init( &argc, &argv );
  
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if(rank == 0) root = true;

  MPI_Comm_size( MPI_COMM_WORLD, &size );
  if(rank != 0) verbose = false;
                                                        
  if(rank == 0){
    //using extra arguments from command line
    for(i = 0; i < argc; i++){
      if(strcmp(argv[i],"-silent") == 0) verbose = false;
      if(strcmp(argv[i],"-file") == 0){
        if((file=fopen(argv[++i], "r")) == NULL) {
          printf("Cannot open file.\n");
          MPI_Finalize();
          return(0);
        }
        else{
          definition = false;
        }
      }
    }
    
    //the input isn't complete
    if(definition){
      printf("Not enought input parameters.\nUse call like: mpiexec -n 1 ./program -file file_name [-silent]\n");
      MPI_Finalize();
      return(0);
    }
    
    //scanning file for number of lines
    fscanf(file,"%d,\n",&lines);
  }
    
  MPI_Bcast((void *)&lines, 1, MPI_INT, 0, MPI_COMM_WORLD);
  matrixA = create_matrix(lines,lines, verbose);
  matrixQ = create_matrix(lines,lines, verbose);
  if(rank == 0){
    for(i = 0; i < lines; i++){
      matrixQ[i][i] = 1;
    }
    //loading of matrix from file to array
    for(i = 0; i < lines; i++){
      for(j = 0; j < lines; j++){
        fscanf(file,"%f,",&matrixA[i][j]);
      }
    }
    if(verbose) printf("Matrix was sucesfuly loaded to memory.\nMartix:\n");
    for(i = 0; i < lines; i++){
      for(j = 0; j < lines; j++){
        if(verbose) printf("%f ",matrixA[i][j]);
      }
      if(verbose) printf("| %f \n",matrixA[i][lines]);
    }
    
    //start time  
    timer.CpuStart();  
  }
  
  int k,l,m,tmpLines,tmpLines2;
  
   
  for(i = 0; i < lines; i++){ 
    tmpLines = (lines-i)/size;
    if(rank == (size-1) && size > 1) tmpLines += (lines-i)%size; 
    tmpLines2 = (lines)/size;
    if(rank == (size-1) && size > 1) tmpLines2 += (lines)%size;  
    //tmp matrix for parallel computing  
    if(i > 0 && i < lines-size){
      delete [] mat[0];
    }
    mat = create_matrix(tmpLines,lines-i, false);
    
    MPI_Bcast((void *)&matrixA[0][0], lines*lines, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast((void *)&matrixQ[0][0], lines*lines, MPI_FLOAT, 0, MPI_COMM_WORLD);      
    // nastaveni pro posilani matice o rozmeru [lines-i][lines-i]
    int piece = (lines-i)/size;
    int radek = lines-i;
    send_counts[0] = (piece)*(radek);
    displs[0] = 0;
    if(size > 1){
      for(k = 1; k < size-1; k++){
        send_counts[k] = piece*(radek);
        displs[k] = displs[k-1] + piece*(radek);
      }
      displs[size-1] = displs[size-2] + piece*(radek);
      send_counts[size-1] = (((lines-i)*(radek)) - displs[size-1]);
    }
    // nastaveni pro posilani matice o rozmeru [lines][lines-i]
    piece = (lines)/size;
    radek = lines-i;
    send_counts2[0] = (piece)*(radek);
    displs2[0] = 0;
    if(size > 1){
      for(k = 1; k < size-1; k++){
        send_counts2[k] = piece*(radek);
        displs2[k] = displs2[k-1] + piece*(radek);
      }
      displs2[size-1] = displs2[size-2] + piece*(radek);
      send_counts2[size-1] = (((lines)*(radek)) - displs2[size-1]);
    }
            
    float x = 0;         
    if(i > 0 && i < lines-size){
      delete [] vec;
      vec = NULL;
    } 
    vec = new float[lines-i]; 
    if(rank == 0){
      for(j = 0; j < lines-i; j++){
        vec[j] = -matrixA[j+i][i];
        x += vec[j]*vec[j];
      }
    
      x = sqrt(x);
      
       
      if(vec[0] > 0) x = -x; 
      vec[0] = vec[0] + x;
      x = 0;
      for(j = 0; j < lines-i; j++){
        x += vec[j]*vec[j];
      }
      x = sqrt(x); 
    }       
    MPI_Bcast((void *)&x, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast((void *)&vec[0], lines-i, MPI_FLOAT, 0, MPI_COMM_WORLD);      
    if(x > 0){
      if(rank == 0){
        //normalizovat vec  
        for(j = 0; j < lines-i; j++){
          vec[j] /= x;
        }
      } 
        
      MPI_Bcast((void *)&vec[0], lines-i, MPI_FLOAT, 0, MPI_COMM_WORLD);  
      //sestavit matici P
      if(i > 0 && i < lines-size){
        delete [] p[0];
        p = NULL;
      }
      p = create_matrix(lines-i,lines-i, false);
      
      MPI_Scatterv(&p[0][0],send_counts,displs,MPI_FLOAT,&mat[0][0],tmpLines*(lines-i),MPI_FLOAT,0,MPI_COMM_WORLD);
      for(k = 0; k < send_counts[rank]/(lines-i); k++){
        for(l = 0; l < lines-i; l++){                              
          if((k+(displs[rank]/(lines-i))) == l) mat[k][l] = 1 - 2*vec[k+displs[rank]/(lines-i)]*vec[l];
          else mat[k][l] = -2*vec[k+displs[rank]/(lines-i)]*vec[l];
        }
      }                            
      MPI_Gatherv(&mat[0][0],send_counts[rank],MPI_FLOAT,&p[0][0],send_counts,displs,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast((void *)&p[0][0], (lines-i)*(lines-i), MPI_FLOAT, 0, MPI_COMM_WORLD); 
      
      //nasobeni matic (paralelizace)
      //R      
      for(k = 0; k < send_counts[rank]/(lines-i); k++){
        for(l = 0; l < lines-i; l++){    
          float tm = 0;
          for(m = i; m < lines; m++){                  
            tm += p[k+displs[rank]/(lines-i)][m-i]*matrixA[m][l+i];   
          }
          mat[k][l] = tm;
        }
      }
      if(i > 0 && i < lines-size){
        delete [] matTmp[0];
        matTmp = NULL;
      }
      matTmp = create_matrix(lines-i,lines-i, false);
      MPI_Gatherv(&mat[0][0],send_counts[rank],MPI_FLOAT,&matTmp[0][0],send_counts,displs,MPI_FLOAT,0,MPI_COMM_WORLD);
      if(rank == 0){
        for(k = i; k < lines; k++){
          for(l = i; l < lines; l++){   
            matrixA[k][l] = matTmp[k-i][l-i];
          }                                
        } 
      }
      
      //Q
      if(i > 0 && i < lines-size){
        delete [] mat[0];
        mat = NULL;
      }
      mat = create_matrix(tmpLines2,lines-i, false);
      for(k = 0; k < send_counts2[rank]/(lines-i); k++){
        for(l = 0; l < lines-i; l++){
          float tm = 0;
          for(m = i; m < lines; m++){
            tm += matrixQ[k+displs2[rank]/(lines-i)][m]*p[m-i][l];
          }
          mat[k][l] = tm;
        }
      }

      if(i > 0 && i < lines-size){
        delete [] matTmp[0];
        matTmp = NULL;
      }
      matTmp = create_matrix(lines,lines-i, false);
      MPI_Gatherv(&mat[0][0],send_counts2[rank],MPI_FLOAT,&matTmp[0][0],send_counts2,displs2,MPI_FLOAT,0,MPI_COMM_WORLD);
      if(rank == 0){
        for(k = 0; k < lines; k++){
          for(l = i; l < lines; l++){
            matrixQ[k][l] = matTmp[k][l-i];
          }
        } 
      }
    }
  }       
    
  //Print solution 
  if(verbose){
    printf("\nSolution is:\n");
    if(verbose) printf("Matrix Q.\n");
    for(k = 0; k < lines; k++){
      for(j = 0; j < lines; j++){
        if(verbose) printf("%f ",matrixQ[k][j]);
      }
      if(verbose) printf("\n");
    } 
    printf("\n");    
    if(verbose) printf("Matrix R.\n");
    for(k = 0; k < lines; k++){
      for(j = 0; j < lines; j++){
        if(verbose) printf("%f ",matrixA[k][j]);
      }
      if(verbose) printf("\n");
    }   
  } 
    
    //Print time     
  if(rank == 0){
    double time;
    time = timer.CpuStop().CpuSeconds();//end time
    if(verbose) printf("Time: ");
    printf("%f\n",time);
  }                   
  MPI_Finalize();
  return(0);                                                                                                
}
