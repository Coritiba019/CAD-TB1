#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// R × C × A matrix
int R, C, A, SEED;

void readConstraints(int *pR, int *pC, int *pA, int *pSEED)
{
    scanf("%d %d %d %d", pR, pC, pA, pSEED);
}

int minimum(int ***mat, int row, int col)
{
    int min = 9999;

    for (int k = 0; k < A; k++)
    {
        if (mat[row][col][k] < min)
            min = mat[row][col][k];
    }

    return min;
}

int maximum(int ***mat, int row, int col)
{
    int max = -1;

    for (int k = 0; k < A; k++)
    {
        if (mat[row][col][k] > max)
            max = mat[row][col][k];
    }

    return max;
}

double median(int ***mat, int row, int col)
{
    double sum = 0;

    for (int k = 0; k < A; k++)
    {
        sum += mat[row][col][k];
    }

    return sum / A;
}

double mean(int ***mat, int row, int col)
{
    double sum = 0;

    for (int k = 0; k < A; k++)
    {
        sum += mat[row][col][k];
    }

    return sum / A;
}

double standardDeviation(int ***mat, int row, int col)
{
    double sum = 0.0, mean, SD = 0.0;

    for (int k = 0; k < A; k++)
    {
        sum += mat[row][col][k];
    }

    mean = sum / A;

    for (int k = 0; k < A; k++)
    {
        SD += pow(mat[row][col][k] - mean, 2);
    }

    return sqrt(SD / (A - 1));
}

int ***alloc_3d_matrix()
{
    int ***mat = (int ***)malloc(R * sizeof(int **));

    if (mat == NULL)
    {
        fprintf(stderr, "Out of memory");
        exit(0);
    }

    for (int i = 0; i < R; i++)
    {
        mat[i] = (int **)malloc(C * sizeof(int *));

        if (mat[i] == NULL)
        {
            fprintf(stderr, "Out of memory");
            exit(0);
        }

        for (int j = 0; j < C; j++)
        {
            mat[i][j] = (int *)malloc(A * sizeof(int));
            if (mat[i][j] == NULL)
            {
                fprintf(stderr, "Out of memory");
                exit(0);
            }
        }
    }
    return mat;
}

void populateMatrix(int ***mat)
{
    srand(SEED);
    // assign values to the allocated memory
    for (int i = 0; i < R; i++)
    {
        for (int j = 0; j < C; j++)
        {
            for (int k = 0; k < A; k++)
            {
                mat[i][j][k] = rand() % 101;
            }
        }
    }
}

void print3dMatrix(int ***mat)
{
    // print the 3D array
    for (int i = 0; i < R; i++)
    {
        for (int j = 0; j < C; j++)
        {
            for (int k = 0; k < A; k++)
            {
                printf(" %3d ", mat[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

void printStats(int ***mat)
{
    // print custom the 3D array
    for (int i = 0; i < R; i++)
    {
        for (int j = 0; j < C; j++)
        {
            for (int k = 0; k < A; k++)
            {
                if (k == 0)
                    printf("Reg %d - Cid %d: menor: %2d, maior: %3d, mediana: %.2lf, media: %.2lf e DP: %.2lf", i, j,
                           minimum(mat, i, j), maximum(mat, i, j), median(mat, i, j), mean(mat, i, j),
                           standardDeviation(mat, i, j));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void free3dMatrix(int ***mat)
{
    // deallocate memory
    for (int i = 0; i < R; i++)
    {
        for (int j = 0; j < C; j++)
        {
            free(mat[i][j]);
        }
        free(mat[i]);
    }
    free(mat);
}

int main(void)
{
    readConstraints(&R, &C, &A, &SEED);
    int ***mat = alloc_3d_matrix();
    populateMatrix(mat);
    print3dMatrix(mat);
    printStats(mat);
    free(mat);
    return 0;
}
