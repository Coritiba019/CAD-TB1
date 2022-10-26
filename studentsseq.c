#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// omp is included JUST to use omp_get_wtime
// and have consistency between the two programs.
#include <omp.h>

typedef struct stats {
    double mean;
    double stddev;
    double median;
    int32_t min;
    int32_t max;
    int32_t region;
    int32_t city;
} stats;

static const char noMemory[] = "Not enough memory available.";

static const int32_t MAX_GRADE = 100;

void readConstraints(int32_t *pR, int32_t *pC, int32_t *pA, int32_t *pSEED) {
    scanf(" %d %d %d %d", pR, pC, pA, pSEED);
}

void *alloc1dMatrix(int32_t x, size_t elemSize) {
    void *mat = (void *)calloc(x, elemSize);

    if (mat == NULL) {
        fprintf(stderr, noMemory);
        exit(ENOMEM);
    }

    return mat;
}

void **alloc2dMatrix(int32_t x, int32_t y, size_t elemSize) {
    void **mat = (void **)malloc(x * sizeof(void *));

    if (mat == NULL) {
        fprintf(stderr, noMemory);
        exit(ENOMEM);
    }

    for (int32_t i = 0; i < x; i++) {
        mat[i] = alloc1dMatrix(y, elemSize);
    }
    return mat;
}

void ***alloc3dMatrix(int32_t x, int32_t y, int32_t z, size_t elemSize) {
    void ***mat = (void ***)malloc(x * sizeof(void **));

    if (mat == NULL) {
        fprintf(stderr, noMemory);
        exit(ENOMEM);
    }

    for (int32_t i = 0; i < x; i++) {
        mat[i] = alloc2dMatrix(y, z, elemSize);
    }
    return mat;
}

void populateMatrix(int32_t ***matrix, int32_t R, int32_t C, int32_t A) {
    // assign values to the allocated memory
    for (int32_t i = 0; i < R; i++) {
        for (int32_t j = 0; j < C; j++) {
            for (int32_t k = 0; k < A; k++) {
                matrix[i][j][k] = rand() % (MAX_GRADE + 1);
            }
        }
    }
}

void print3dMatrix(int32_t ***mat, int32_t R, int32_t C, int32_t A) {
    // print the 3D array
    for (int32_t i = 0; i < R; i++) {
        for (int32_t j = 0; j < C; j++) {
            for (int32_t k = 0; k < A; k++) {
                printf(" %3d ", mat[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

void free2dMatrix(void **mat, int32_t x) {
    for (int32_t i = 0; i < x; i++) {
        free(mat[i]);
    }
    free(mat);
}

void free3dMatrix(void ***mat, int32_t x, int32_t y) {
    for (int32_t i = 0; i < x; i++) {
        free2dMatrix(mat[i], y);
    }
    free(mat);
}

double meanFreqArray(int64_t *freqArr, int64_t n) {
    double sum = 0;

    for (int32_t i = 0; i <= MAX_GRADE; i++) {
        sum += freqArr[i] * i;
    }

    return sum / n;
}

double stdDevFreqArray(int64_t *freq, double mean, int64_t n) {
    double sum = 0;

    for (int32_t i = 0; i <= MAX_GRADE; i++) {
        sum += freq[i] * pow(i - mean, 2);
    }

    return sqrt(sum / (n - 1));
}

int32_t maxFreqArray(int64_t *freqArr) {
    for (int32_t i = MAX_GRADE; i >= 0; i--) {
        if (freqArr[i] != 0) {
            return i;
        }
    }
    return -1;
}

int32_t minFreqArray(int64_t *freqArr) {
    for (int32_t i = 0; i <= MAX_GRADE; i++) {
        if (freqArr[i] != 0) {
            return i;
        }
    }
    return -1;
}

double medianFreqArray(int64_t *freqArr, int64_t n) {

    int64_t currN = 0;
    int64_t botMid = (n - 1) / 2;
    int64_t topMid = n / 2;

    int32_t botMidGrade = -1;
    int32_t topMidGrade = -1;

    for (int32_t i = 0; i <= MAX_GRADE; i++) {
        currN += freqArr[i];

        // Minus 1 to align with the array index (starts at 0)

        if ((currN - 1) >= botMid && botMidGrade == -1) {
            botMidGrade = i;
        }
        if ((currN - 1) >= topMid) {
            topMidGrade = i;
            break;
        }
    }

    return n % 2 == 0 ? (botMidGrade + topMidGrade) / 2.0 : (double)topMidGrade;
}

int64_t ***calculateCityFreqArray(int32_t ***infoMatrix, int32_t R, int32_t C,
                                  int32_t A) {

    int64_t ***cityFreqArray =
        (int64_t ***)alloc3dMatrix(R, C, MAX_GRADE + 1, sizeof(int64_t));

    int32_t grade;

    for (int32_t reg = 0; reg < R; reg++) {
        for (int32_t city = 0; city < C; city++) {
            for (int32_t student = 0; student < A; student++) {
                grade = infoMatrix[reg][city][student];
                cityFreqArray[reg][city][grade]++;
            }
        }
    }

    return cityFreqArray;
}

int64_t **calculateRegionFreqArray(int64_t ***cityFreqArray, int32_t R,
                                   int32_t C) {

    int64_t **regionFreqArray =
        (int64_t **)alloc2dMatrix(R, MAX_GRADE + 1, sizeof(int64_t));

    for (int32_t reg = 0; reg < R; reg++) {
        for (int32_t city = 0; city < C; city++) {
            for (int32_t grade = 0; grade <= MAX_GRADE; grade++) {
                regionFreqArray[reg][grade] += cityFreqArray[reg][city][grade];
            }
        }
    }

    return regionFreqArray;
}

int64_t *calculateCountryFreqArray(int64_t **regionFreqArray, int32_t R) {

    int64_t *countryFreqArray = alloc1dMatrix(MAX_GRADE + 1, sizeof(int64_t));

    for (int32_t reg = 0; reg < R; reg++) {
        for (int32_t grade = 0; grade <= MAX_GRADE; grade++) {
            countryFreqArray[grade] += regionFreqArray[reg][grade];
        }
    }

    return countryFreqArray;
}

stats *getStatFromFreqArray(int64_t *freqArr, int32_t n) {
    stats *stat = (stats *)malloc(sizeof(stats));
    if (stat == NULL) {
        fprintf(stderr, noMemory);
        exit(ENOMEM);
    }

    stat->mean = meanFreqArray(freqArr, n);
    stat->stddev = stdDevFreqArray(freqArr, stat->mean, n);
    stat->max = maxFreqArray(freqArr);
    stat->min = minFreqArray(freqArr);
    stat->median = medianFreqArray(freqArr, n);

    return stat;
}

void getStats(int32_t R, int32_t C, int32_t A, int64_t ***cityFreqArray,
              int64_t **regFreqArray, int64_t *countryFreqArray,
              stats **statsVector, stats **bestReg, stats **bestCity) {
    int32_t idx = 0;
    for (int32_t reg = 0; reg < R; reg++) {
        for (int32_t city = 0; city < C; city++) {
            statsVector[idx] =
                getStatFromFreqArray(cityFreqArray[reg][city], A);
            statsVector[idx]->region = reg;
            statsVector[idx]->city = city;

            if ((*bestCity) == NULL ||
                statsVector[idx]->mean > (*bestCity)->mean) {
                *bestCity = statsVector[idx];
            }

            idx++;
        }
    }
    for (int32_t reg = 0; reg < R; reg++) {
        statsVector[idx] = getStatFromFreqArray(regFreqArray[reg], A * C);
        statsVector[idx]->region = reg;
        statsVector[idx]->city = -1;

        if ((*bestReg) == NULL || statsVector[idx]->mean > (*bestReg)->mean) {
            *bestReg = statsVector[idx];
        }

        idx++;
    }
    statsVector[idx] = getStatFromFreqArray(countryFreqArray, A * C * R);
    statsVector[idx]->region = -1;
    statsVector[idx]->city = -1;
}

void printStats(stats **stats, int32_t R, int32_t C) {

    int32_t idx = 0;

    for (int32_t i = 0; i < R; i++) {
        for (int32_t j = 0; j < C; j++) {

            printf("Reg %d - Cid %d: menor: %d, maior: %d, mediana: "
                   "%.2lf, média: %.2lf e DP: %.2lf\n",
                   stats[idx]->region, stats[idx]->city, stats[idx]->min,
                   stats[idx]->max, stats[idx]->median, stats[idx]->mean,
                   stats[idx]->stddev);

            idx++;
        }
        printf("\n\n");
    }

    for (int32_t i = 0; i < R; i++) {
        printf("Reg %d: menor: %d, maior: %d, mediana: %.2lf, média: "
               "%.2lf e DP: %.2lf\n",
               stats[idx]->region, stats[idx]->min, stats[idx]->max,
               stats[idx]->median, stats[idx]->mean, stats[idx]->stddev);

        idx++;
    }
    printf("\n\n");

    printf("Brasil: menor: %d, maior: %d, mediana: %.2lf, média: %.2lf e "
           "DP: %.2lf\n",
           stats[idx]->min, stats[idx]->max, stats[idx]->median,
           stats[idx]->mean, stats[idx]->stddev);

    printf("\n\n");
}

int main(void) {
    // R × C × A matrix
    int32_t R, C, A, SEED;
    double timeStart, timeEnd;
    readConstraints(&R, &C, &A, &SEED);

    srand(SEED);

    int32_t ***infoMatrix =
        (int32_t ***)alloc3dMatrix(R, C, A, sizeof(int32_t));
    populateMatrix(infoMatrix, R, C, A);
    // print3dMatrix(infoMatrix, R, C, A);

    timeStart = omp_get_wtime();
    stats *bestCity = NULL;
    stats *bestReg = NULL;

    int64_t ***cityFreqArray = calculateCityFreqArray(infoMatrix, R, C, A);
    int64_t **regFreqArray = calculateRegionFreqArray(cityFreqArray, R, C);
    int64_t *countryFreqArray = calculateCountryFreqArray(regFreqArray, R);

    stats **statsVec = alloc1dMatrix((R * C) + R + 1, sizeof(stats *));
    getStats(R, C, A, cityFreqArray, regFreqArray, countryFreqArray, statsVec,
             &bestReg, &bestCity);

    timeEnd = omp_get_wtime();

    printStats(statsVec, R, C);

    printf("Melhor região: Região %d\n", bestReg->region);
    printf("Melhor cidade: Região %d, Cidade %d\n", bestCity->region,
           bestCity->city);

    printf("Tempo de resposta sem considerar E/S, em segundos: %lfs\n",
           timeEnd - timeStart);

    free3dMatrix((void ***)infoMatrix, R, C);
    free3dMatrix((void ***)cityFreqArray, R, C);
    free2dMatrix((void **)regFreqArray, R);
    free(countryFreqArray);
    free2dMatrix((void **)statsVec, (R * C) + R + 1);
    return 0;
}
