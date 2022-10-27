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

typedef struct maxMeanIdx {
    double mean;
    int64_t idx;
} maxMeanIdx;

#define flat2d(i, j, Y) ((i) * (Y) + (j))
#define flat3d(i, j, k, Y, Z) ((i) * (Y) * (Z) + (j) * (Z) + (k))
#define flat4d(i, j, k, l, Y, Z, W)                                            \
    ((i) * (Y) * (Z) * (W) + (j) * (Z) * (W) + (k) * (W) + (l))

static const char noMemory[] = "Not enough memory available.";

static const int32_t MAX_GRADE = 100;

void readConstraints(int32_t *pR, int32_t *pC, int32_t *pA, int32_t *pSEED) {
    scanf(" %d %d %d %d", pR, pC, pA, pSEED);
}

void *allocMatrix(int32_t x, size_t elemSize) {
    void *mat = (void *)calloc(x, elemSize);

    if (mat == NULL) {
        fprintf(stderr, noMemory);
        exit(ENOMEM);
    }

    return mat;
}

void populateMatrix(int32_t *matrix, int32_t R, int32_t C, int32_t A) {
    // assign values to the allocated memory
    for (int32_t i = 0; i < R; i++) {
        for (int32_t j = 0; j < C; j++) {
            for (int32_t k = 0; k < A; k++) {
                matrix[flat3d(i, j, k, C, A)] = rand() % (MAX_GRADE + 1);
            }
        }
    }
}

void print3dMatrix(int32_t *mat, int32_t R, int32_t C, int32_t A) {
    // print the 3D array
    for (int32_t i = 0; i < R; i++) {
        for (int32_t j = 0; j < C; j++) {
            for (int32_t k = 0; k < A; k++) {
                printf(" %3d ", mat[flat3d(i, j, k, C, A)]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

double meanFreqArray(int64_t *freqArr, int64_t n) {
    double sum = 0;

#pragma omp simd reduction(+ : sum)
    for (int32_t i = 0; i <= MAX_GRADE; i++) {
        sum += freqArr[i] * i;
    }

    return sum / n;
}

double stdDevFreqArray(int64_t *freq, double mean, int64_t n) {
    double sum = 0;

#pragma omp simd reduction(+ : sum)
    for (int32_t i = 0; i <= MAX_GRADE; i++) {
        sum += freq[i] * (i - mean) * (i - mean);
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

int64_t *calculateCityFreqArray(int32_t *infoMatrix, int32_t R, int32_t C,
                                int32_t A) {

    int32_t nT = omp_get_max_threads();

    int64_t *cityFreqArray =
        (int64_t *)allocMatrix(R * C * (MAX_GRADE + 1), sizeof(int64_t));

    int64_t *threadFreqArray =
        (int64_t *)allocMatrix(nT * R * C * (MAX_GRADE + 1), sizeof(int64_t));

#pragma omp parallel
    {
        int32_t tid = omp_get_thread_num();

#pragma omp for collapse(3)
        for (int32_t reg = 0; reg < R; reg++) {
            for (int32_t city = 0; city < C; city++) {
                for (int32_t student = 0; student < A; student++) {

                    int32_t grade =
                        infoMatrix[flat3d(reg, city, student, C, A)];

                    threadFreqArray[flat4d(tid, reg, city, grade, R, C,
                                           MAX_GRADE + 1)]++;
                }
            }
        }

#pragma omp for collapse(3) nowait
        for (int32_t reg = 0; reg < R; reg++) {
            for (int32_t city = 0; city < C; city++) {
                for (int32_t grade = 0; grade <= MAX_GRADE; grade++) {
                    for (int32_t thread = 0; thread < nT; thread++) {
                        cityFreqArray[flat3d(reg, city, grade, C,
                                             MAX_GRADE + 1)] +=
                            threadFreqArray[flat4d(thread, reg, city, grade, R,
                                                   C, MAX_GRADE + 1)];
                    }
                }
            }
        }
    }

    free(threadFreqArray);

    return cityFreqArray;
}

int64_t *calculateRegionFreqArray(int64_t *cityFreqArray, int32_t R,
                                  int32_t C) {

    int64_t *regionFreqArray =
        (int64_t *)allocMatrix(R * (MAX_GRADE + 1), sizeof(int64_t));

#pragma omp parallel for collapse(2) reduction(+ : regionFreqArray[:R*(MAX_GRADE+1)])
    for (int32_t reg = 0; reg < R; reg++) {
        for (int32_t city = 0; city < C; city++) {
            for (int32_t grade = 0; grade <= MAX_GRADE; grade++) {
                regionFreqArray[flat2d(reg, grade, MAX_GRADE + 1)] +=
                    cityFreqArray[flat3d(reg, city, grade, C, MAX_GRADE + 1)];
            }
        }
    }

    return regionFreqArray;
}

int64_t *calculateCountryFreqArray(int64_t *regionFreqArray, int32_t R) {

    int64_t *countryFreqArray = allocMatrix(MAX_GRADE + 1, sizeof(int64_t));

#pragma omp parallel for reduction(+ : countryFreqArray[:MAX_GRADE + 1])
    for (int32_t reg = 0; reg < R; reg++) {
        for (int32_t grade = 0; grade <= MAX_GRADE; grade++) {
            countryFreqArray[grade] +=
                regionFreqArray[flat2d(reg, grade, MAX_GRADE + 1)];
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

void getStats(int32_t R, int32_t C, int32_t A, int64_t *cityFreqArray,
              int64_t *regFreqArray, int64_t *countryFreqArray,
              stats **statsVector, stats **bestReg, stats **bestCity) {

#pragma omp parallel for
    for (int64_t i = 0; i < R * C + R + 1; i++) {
        if (i < R * C) {
            int32_t reg = i / C;
            int32_t city = i % C;

            statsVector[i] = getStatFromFreqArray(
                &cityFreqArray[flat3d(reg, city, 0, C, MAX_GRADE + 1)], A);
            statsVector[i]->region = reg;
            statsVector[i]->city = city;

        } else if (i < R * C + R) {
            int32_t reg = i - R * C;

            statsVector[i] = getStatFromFreqArray(
                &regFreqArray[flat2d(reg, 0, MAX_GRADE + 1)], C * A);
            statsVector[i]->region = reg;
            statsVector[i]->city = -1;

        } else {
            statsVector[i] = getStatFromFreqArray(countryFreqArray, A * C * R);
            statsVector[i]->region = -1;
            statsVector[i]->city = -1;
        }
    }

#pragma omp declare reduction(                                                 \
    maxMean                                                                    \
    : struct maxMeanIdx                                                        \
    : omp_out = omp_in.mean > omp_out.mean                                     \
                    ? omp_in                                                   \
                    : (omp_out.mean == omp_in.mean                             \
                           ? omp_in.idx < omp_out.idx ? omp_in : omp_out       \
                           : omp_out))

    maxMeanIdx bestCityIdx = {-1.0, -1};

#pragma omp parallel for reduction(maxMean : bestCityIdx)
    for (int64_t i = 0; i < R * C; i++) {
        if (statsVector[i]->mean > bestCityIdx.mean) {
            bestCityIdx.mean = statsVector[i]->mean;
            bestCityIdx.idx = i;
        }
    }

    maxMeanIdx bestRegIdx = {-1.0, -1};

#pragma omp parallel for reduction(maxMean : bestRegIdx)
    for (int64_t i = R * C; i < R * C + R; i++) {
        if (statsVector[i]->mean > bestRegIdx.mean) {
            bestRegIdx.mean = statsVector[i]->mean;
            bestRegIdx.idx = i;
        }
    }

    *bestCity = statsVector[bestCityIdx.idx];
    *bestReg = statsVector[bestRegIdx.idx];
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

    int32_t *infoMatrix = (int32_t *)allocMatrix(R * C * A, sizeof(int32_t));
    populateMatrix(infoMatrix, R, C, A);
    // print3dMatrix(infoMatrix, R, C, A);

    timeStart = omp_get_wtime();

    stats *bestCity = NULL;
    stats *bestReg = NULL;

    int64_t *cityFreqArray = calculateCityFreqArray(infoMatrix, R, C, A);
    int64_t *regFreqArray = calculateRegionFreqArray(cityFreqArray, R, C);
    int64_t *countryFreqArray = calculateCountryFreqArray(regFreqArray, R);

    stats **statsVec = allocMatrix((R * C) + R + 1, sizeof(stats *));
    getStats(R, C, A, cityFreqArray, regFreqArray, countryFreqArray, statsVec,
             &bestReg, &bestCity);

    timeEnd = omp_get_wtime();

    printStats(statsVec, R, C);

    printf("Melhor região: Região %d Média: %lf\n", bestReg->region,
           bestReg->mean);
    printf("Melhor cidade: Região %d, Cidade %d Média: %lf\n", bestCity->region,
           bestCity->city, bestCity->mean);

    // printf("Tempo de resposta sem considerar E/S, em segundos: %lfs\n",
    //        timeEnd - timeStart);

    FILE *timeFile = fopen("timePar", "w");

    fprintf(timeFile, "%lf", timeEnd - timeStart);

    fclose(timeFile);

    free(infoMatrix);
    free(cityFreqArray);
    free(regFreqArray);
    free(countryFreqArray);
    for (int32_t i = 0; i < (R * C) + R + 1; i++) {
        free(statsVec[i]);
    }
    free(statsVec);
    return 0;
}
