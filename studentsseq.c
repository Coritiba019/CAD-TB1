#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// omp is included JUST to use omp_get_wtime
// and have consistency between the two programs.
#include <omp.h>

// R × C × A matrix
int R, C, A, SEED;
const int MAX_GRADE = 101;

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

int ***alloc3dMatrix(int x, int y, int z)
{
    int ***mat = (int ***)malloc(x * sizeof(int **));

    if (mat == NULL)
    {
        fprintf(stderr, "Out of memory");
        exit(0);
    }

    for (int i = 0; i < x; i++)
    {
        mat[i] = (int **)malloc(y * sizeof(int *));

        if (mat[i] == NULL)
        {
            fprintf(stderr, "Out of memory");
            exit(0);
        }

        for (int j = 0; j < y; j++)
        {
            mat[i][j] = (int *)malloc(z * sizeof(int));
            if (mat[i][j] == NULL)

            {
                fprintf(stderr, "Out of memory");
                exit(0);
            }

            for (int k = 0; k < z; k++)
                mat[i][j][k] = 0;
        }
    }
    return mat;
}

int **alloc2dMatrix(int x, int y)
{
    int **mat = (int **)malloc(x * sizeof(int *));

    if (mat == NULL)
    {
        fprintf(stderr, "Out of memory");
        exit(0);
    }

    for (int i = 0; i < x; i++)
    {
        mat[i] = (int *)malloc(y * sizeof(int));

        if (mat[i] == NULL)
        {
            fprintf(stderr, "Out of memory");
            exit(0);
        }

        for (int j = 0; j < y; j++)
            mat[i][j] = 0;
    }
    return mat;
}

int *alloc1dMatrix(int x)
{
    int *mat = (int *)malloc(x * sizeof(int));

    if (mat == NULL)
    {
        fprintf(stderr, "Out of memory");
        exit(0);
    }

    for (int i = 0; i < x; i++)
        mat[i] = 0;

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

// TODO: Accept parameters.
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

int *randomFreqArray()
{
    srand(omp_get_wtime());
    int *freqArr = malloc(sizeof(int) * 101);
    for (int i = 0; i <= 100; i++)
        freqArr[i] = rand() % 100;
    return freqArr;
}

double meanFreqArray(int *freqArr, long long numberOfStudents)
{
    double sum = 0;
    for (int i = 0; i <= 100; i++)
    {
        sum += freqArr[i] * i;
    }
    double mean = sum / numberOfStudents;
    return mean;
}

double sdFreqArray(int *freq, double mean)
{
    // todo
}

int maxFreqArray(int *freqArr)
{
    for (int i = 100; i >= 0; i--)
    {
        if (freqArr[i] != 0)
        {
            return i;
        }
    }
    fprintf(stderr, "Frequency array is empty");
    exit(0);
}

int minFreqArray(int *freqArr)
{
    for (int i = 0; i <= 100; i++)
    {
        if (freqArr[i] != 0)
        {
            return i;
        }
    }
    fprintf(stderr, "Frequency array is empty");
    exit(0);
}

int medianFreqArray(int *freqArr, long long numberOfStudents)
{
    // Calculate the number of students in the frequency array
    long long currentNumberOfStudents = 0;
    for (int i = 0; i <= 100; i++)
    {
        currentNumberOfStudents += freqArr[i];
        if (currentNumberOfStudents * 2 >= numberOfStudents)
        {
            return i;
        }
    }
}

int ***calculateCityFreqArray(int ***mat)
{
    int ***freqArray = alloc3dMatrix(R, C, MAX_GRADE);

    for (int reg = 0; reg < R; reg++)
    {
        for (int city = 0; city < C; city++)
        {
            for (int student = 0; student < A; student++)
            {
                freqArray[reg][city][mat[reg][city][student]]++;
            }
        }
    }

    return freqArray;
}

int **calculateRegionFreqArray(int ***cityFreqArray)
{
    int **freqArray = alloc2dMatrix(R, MAX_GRADE);

    for (int reg = 0; reg < R; reg++)
    {
        for (int city = 0; city < C; city++)
        {
            for (int ind = 0; ind <= 100; ind++)
            {
                freqArray[reg][ind] += cityFreqArray[reg][city][ind];
            }
        }
    }

    return freqArray;
}

int *calculateNationalFreqArray(int **regionFreqArray)
{
    int *freqArray = alloc1dMatrix(MAX_GRADE);

    for (int reg = 0; reg < R; reg++)
    {
        for (int ind = 0; ind <= 100; ind++)
        {
            freqArray[ind] += regionFreqArray[reg][ind];
        }
    }

    return freqArray;
}

int main(void)
{
    readConstraints(&R, &C, &A, &SEED);
    int ***mat = alloc3dMatrix(R, C, A);
    populateMatrix(mat);
    print3dMatrix(mat);
    double tic = omp_get_wtime();
    printStats(mat);
    double toc = omp_get_wtime();
    // TODO: Shouldn't count time printing
    printf("Time elapsed: %fs\n", toc - tic);

    int ***cityFreqArray = calculateCityFreqArray(mat);
    int **regionFreqArray = calculateRegionFreqArray(cityFreqArray);
    int *nationalFreqArray = calculateNationalFreqArray(regionFreqArray);

    free3dMatrix(mat);
    free3dMatrix(cityFreqArray);
    // TODO: free2dMatrix(regionFreqArray);
    free(nationalFreqArray);
    return 0;
}
