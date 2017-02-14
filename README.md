# QR-decomposition
Parallelization of QR decomposition with Householder transformation

/*
    GQ Decomposition & Eigenvalues by Ondrej Zjevik, 2012
*/

Project structure

[eig]
  - contains everything needed to run the code to find eigenvalues of a matrix
  
  run_mpi.sh 
  
  run_serial.sh
    - these scripts run the programs to find the eigenvalues
    
[qr]
  - contains everything needed to run the code to find a QR decomposition of a matrix
  
  run_mpi.sh 
  
  run_omp.sh 
  
  run_serial.sh
  
    - these scripts contains examples how to run and they run the programs to find the QR decomposition of a matrix  
