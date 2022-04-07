//Pthread Matrix Multiplication 
//COMS 480
//Michael Davis

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

void matMult(double **A, double **B, double**C, int n)
{
  
  for( int r = 0; r < n; r++)
      for(int c = 0; c < n; c++)
	for(int k = 0; k < n; k++)
	  C[r][c] += A[r][k] * B[k][c];       
}

void blockMatMulti(double **A, double **B, double**C, int n, int blockSize)
{
  
  double sum;
  int br, bc, r, c, k;
  int myMin(int a, int b)
  {
    return (a > b) ? b : a;
  }
  
  for (br = 0; br < n; br += blockSize){
      for(bc = 0; bc < n; bc += blockSize){
	  for(r = 0; r < n; ++r){
	      for(c = br; c < myMin(br + blockSize, n); ++c){
		  sum = 0.0;
		  for(k = bc; k < myMin(bc + blockSize, n); ++k)
		    sum += A[r][k] * B[k][c];
		    C[r][c] += sum;
		}
	    }
	}
    }
}

typedef struct {
  double** A;
  double** B;
  double** C;
  int n;
  int threadID;
  int numThreads;
  
}myData;

void* matMultPthreadsFunct(void* arg)
{
  int r, c;
  myData* theMats = (myData*)arg;
  double** A = theMats->A;
  double** B = theMats->B;
  double** C = theMats->C;
  long threadid = theMats->threadID;
  long numThreads = theMats->numThreads;
  long size = theMats->n;

  long sRow = threadid * size / numThreads;
  long eRow = (threadid + 1) * size / numThreads - 1;

  for (r = sRow; r <= eRow; r++)
    for(c = 0; c < size; c++)
      C[r][c] = 0.0;
      for( int i = 0; i < size; i++)
	C[r][c] += A[r][i] * B[i][c];
      
  return NULL;
}

double** pthreadMatMulti(double** mat1, double** mat2, double** mat5, int n, int numThreads)
{
  pthread_t tid[numThreads];
  myData Matrices[numThreads];
  
  for(int i = 0; i < numThreads; i++)
    {
      Matrices[i].A = mat1;
      Matrices[i].B = mat2;
      Matrices[i].C = mat5;
      Matrices[i].n = n;
      Matrices[i].threadID = i;
      Matrices[i].numThreads = numThreads;

      
      pthread_create(&tid[i],NULL , matMultPthreadsFunct,(void*)(&Matrices[i-1]));
    }

  for(int i = 0; i < numThreads; i++)
    pthread_join(tid[i], NULL);
 }

int main()
{
  struct timespec start, end;
  
  int size, r, c;
  clock_t t;
  int nthread;
  
  printf("Enter size: ");
  scanf("%d",&size);

  int nRows = size;
  int nCols = size;

  printf("Enter Thread Count: ");
  scanf("%d", &nthread);
  
    //Allocate Matrix 1
      double** mat1 = (double**)malloc(nRows * sizeof(double*));
      mat1[0] = (double*)malloc(nRows * nCols * sizeof(double));
      for(r = 1; r < nRows; ++r){
	mat1[r] = &(mat1[0][r * nCols]);
      }

  //Allocate Matrix 2
      double** mat2 = (double**)malloc(nRows * sizeof(double*));
      mat2[0] = (double*)malloc(nRows * nCols * sizeof(double));
      for(r = 1; r < nRows; ++r){
	mat2[r] = &(mat2[0][r * nCols]);
      }

  //Allocate Matrix 3
      double** mat3 = (double**)malloc(nRows * sizeof(double*));
      mat3[0] = (double*)malloc(nRows * nCols * sizeof(double));
      for(r = 1; r < nRows; ++r){
	mat3[r] = &(mat3[0][r * nCols]);
      }
 
  //Allocate Matrix 4
      double** mat4 = (double**)malloc(nRows * sizeof(double*));
      mat4[0] = (double*)malloc(nRows * nCols * sizeof(double));
      for(r = 1; r < nRows; ++r){
	mat4[r] = &(mat4[0][r * nCols]);
      }

  //Allocate Matrix 5
      double** mat5 = (double**)malloc(nRows * sizeof(double*));
      mat5[0] = (double*)malloc(nRows * nCols * sizeof(double));
      for(r = 1; r < nRows; ++r){
	mat5[r] = &(mat5[0][r * nCols]);
      }
  //Fill Matrix
      for(int r = 0; r < nRows; ++r)
	for(int c = 0; c < nCols; ++c){
	  mat1[r][c] = r * nCols + c + 1;
	}
  // Fill Matrix 2
      for(int r = nRows - 1 ; r >= 0; r--)
	for(int c = nCols - 1; c >= 0; c--)
	  mat2[nRows - 1 - r][nCols - 1 - c] = r * nCols + c + 1;


	  // Testing all functions 

      clock_gettime(CLOCK_MONOTONIC, &start);
      
      matMult(mat1, mat2, mat3, size);

      clock_gettime(CLOCK_MONOTONIC, &end);

      double matmulti_timetaken;

      matmulti_timetaken = (end.tv_sec - start.tv_sec) * 1e9;
      matmulti_timetaken = (matmulti_timetaken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
      
      printf("\nSerial Matrix Mutli: %f seconds\n", matmulti_timetaken);
      
      clock_gettime(CLOCK_MONOTONIC, &start);
      
      blockMatMulti(mat1, mat2, mat4, size, nthread);

      clock_gettime(CLOCK_MONOTONIC, &end);

      double blockMat_time;

      blockMat_time = (end.tv_sec - start.tv_sec) * 1e9;
      blockMat_time = (blockMat_time + (end.tv_nsec - start.tv_nsec)) * 1e-9;

      printf("\nBlock Matrix Mutli: %f seconds\n", blockMat_time);

      clock_gettime(CLOCK_MONOTONIC, &start);

      pthreadMatMulti(mat1, mat2, mat5 , size, nthread);

      clock_gettime(CLOCK_MONOTONIC, &end);

      double pthread_time;

      pthread_time = (end.tv_sec - start.tv_sec) * 1e9;
      pthread_time = (pthread_time + (end.tv_nsec - start.tv_nsec)) * 1e-9;

      printf("\nP Thread Multi: %f seconds\n", pthread_time);

      free(mat1);
      free(mat2);
      free(mat3);
      free(mat4);
      free(mat5);
  
  return 0;

}
