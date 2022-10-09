#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
// M × N × O matrix
#define M 3
#define N 4
#define O 6

int menor(int *** A, int row, int col)
{
  int menor_num = 9999;

    for (int k = 0; k < O; k++) 
    {
      if(A[row][col][k] < menor_num) 
      {
        menor_num = A[row][col][k];
      }
    }

  return menor_num;
}

int maior(int *** A, int row, int col)
{
  int maior_num = -1;

    for (int k = 0; k < O; k++) 
    {
      if(A[row][col][k] > maior_num) 
      {
        maior_num = A[row][col][k];
      }
    }

  return maior_num;
}

double mediana(int *** A, int row, int col)
{
  double sum = 0;

    for (int k = 0; k < O; k++) 
    {
      sum += A[row][col][k];
    }

  return sum/O;
}

double media(int *** A, int row, int col)
{
  double sum = 0;

    for (int k = 0; k < O; k++) 
    {
      sum += A[row][col][k];
    }

  return sum/O;
}

double calculateSD(int *** A, int row, int col) 
{
    double sum = 0.0, mean, SD = 0.0;
    for (int k = 0; k < O; k++) 
    {
      sum += A[row][col][k];
    }
    mean = sum / O;
    for (int k = 0; k < O; k++) 
    {
      SD += pow(A[row][col][k] - mean, 2);
    }
    return sqrt(SD / (O-1));
}
 
// Dynamically allocate memory for 3D Array
int main(void)
{
    int*** A = (int***)malloc(M * sizeof(int**));
 
    if (A == NULL)
    {
        fprintf(stderr, "Out of memory");
        exit(0);
    }
 
    for (int i = 0; i < M; i++)
    {
        A[i] = (int**)malloc(N * sizeof(int*));
 
        if (A[i] == NULL)
        {
            fprintf(stderr, "Out of memory");
            exit(0);
        }
 
        for (int j = 0; j < N; j++)
        {
            A[i][j] = (int*)malloc(O * sizeof(int));
            if (A[i][j] == NULL)
            {
                fprintf(stderr, "Out of memory");
                exit(0);
            }
        }
    }
 
    srand(7);
    // assign values to the allocated memory
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < O; k++) 
            {
                A[i][j][k] = rand() % 101;
            }
        }
    }

    // print the 3D array
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < O; k++) {
                printf("%d ", A[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
 
    // print custom the 3D array
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < O; k++) 
            {
                if(k==0) printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.2lf, media: %.2lf e DP: %.2lf", i, j, menor(A, i, j), maior(A, i, j), mediana(A, i, j), media(A, i, j), calculateSD(A, i, j));
            }
            printf("\n");
        }
        printf("\n");
    }
 
    // deallocate memory
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++) {
            free(A[i][j]);
        }
        free(A[i]);
    }
    free(A);
 
    return 0;
}
