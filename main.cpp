#include <iostream>
#include <string>
#include <vector>
#include <string.h>
#include <random>
#include <omp.h>
#include <unistd.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include "BoostBarrier.h"
#include "ind.h"
#include "Worker.h"
#include "Worker2.h"
#include "Worker3.h"
#include "uint.h"
#include "ConcurrentQueue.h"


int minimum(int a, int b, int c) {
    int min = a;

    if (b < min) min = b;
    if (c < min) min = c;

    return min;
}

int minimum_(const int a, const int b, const int c) {
    return std::min(std::min(a, b), c);
}


int editDistanceST(const char *x, const char *y) {
    const int M = strlen(x);
    const int N = strlen(y);

    int distance;

    uint16 **D = new uint16 *[M + 1];

    int i, j;

    for (i = 0; i <= M; i++)
        D[i] = new uint16[N + 1];

    for (i = 0; i <= M; i++)
        D[i][0] = i;

    for (j = 1; j <= N; j++)
        D[0][j] = j;

    float t2, t1;
    t1 = omp_get_wtime();
    for (i = 1; i <= M; i++) {
        for (j = 1; j <= N; j++) {
            if (x[i - 1] != y[j - 1]) {
                uint16 k = minimum_(D[i][j - 1], //insertion
                                 D[i - 1][j], //insertion
                                 D[i - 1][j - 1]); //substitution
                D[i][j] = k + 1;
            } else {
                D[i][j] = D[i - 1][j - 1];
            }
        }
    }
    t2 = omp_get_wtime();

    //printf("Pure ST computing time: %f\n", t2-t1);
    distance = D[M][N];
/*
    for(i = 0; i < M+1; i++){
        for(j = 0; j < N+1; j++){
            std::cout << D[i][j] << ", ";
        }
        std::cout << std::endl;
    }
*/
    for (i = 0; i <= M; i++)
        delete D[i];

    delete[] D;
    return distance;
}

void printMatrix(uint16 const *D, const uint16 M, const uint16 N) {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            printf("%i, ", D[i * M + j]);
        }
        std::cout << std::endl;
    }
}

const int MAX_NUM_THREADS = 12;
const uint TW = 512;

void computeSubMatrix(uint I, uint J, const int M, const int N, const char *x, const char *y, uint16 *D) {
    uint M_ = M + 1;
    uint N_ = N + 1;
    I = I * TW + 1;
    J = J * TW + 1;
    for (int i = I; i < N_ && i < I + TW; i++) {
        for (int j = J; j < M_ && j < J + TW; j++) {
            if (x[i - 1] != y[j - 1]) {
                uint16 k = minimum(D[i * M_ + j - 1], //insertion
                                D[(i - 1) * M_ + j], //insertion
                                D[(i - 1) * M_ + j - 1]); //substitution
                D[i * M_ + j] = k + 1;
            } else {
                D[i * M_ + j] = D[(i - 1) * M_ + j - 1];
            }
        }
    }
}

int editDistanceOMP(const char *x, const char *y) {
    const uint M = strlen(x);
    const uint N = strlen(y);

    const uint M_ = M + 1;
    const uint N_ = N + 1;

    //std::vector<uint16> D(M_*N_);
    const uint dim = (uint)M_*(uint)N_;
    uint16 * D = new uint16[dim];
    //D.resize(M_ * N_);

    uint i, j;

    for (i = 0; i < M_; i++)
        D[i] = i;

    for (j = 1; j < N_; j++)
        D[j * M_] = j;

    ////////////////

    int M_Tiles = ceil((float) M / TW);
    int N_Tiles = ceil((float) N / TW);

    const int DStart = -M_Tiles + 1;
    const int DFinish = N_Tiles;
    const int TILE_WIDTH = Worker::TILE_WIDTH;


    const int TOTAL_TILES = M_Tiles * N_Tiles;
    float t1, t2;

#pragma omp parallel num_threads(MAX_NUM_THREADS)
    {
#pragma omp master
        {
            for (int d = DStart; d < DFinish; d++) {
                const int iMax = std::min(M_Tiles + d, N_Tiles);
                const int iMin = std::max(d, 0);
                for (i = iMin; i < iMax; i++) {
#pragma omp task
                    {
                        j = M_Tiles - i + d - 1;
                        computeSubMatrix(i, j, M, N, x, y, &D[0]);
                    }
                }
#pragma omp taskwait
            }
        }
    }

    //b.count_down_and_wait();
    int distance = D[N * M_ + M];

    delete[] D;
    return distance;
}

int editDistanceCPPT(const char *x, const char *y) {
    const int M = strlen(x);
    const int N = strlen(y);

    const uint M_ = M + 1;
    const uint N_ = N + 1;

    const uint dim = (uint)M_*(uint)N_;
    uint16 * D = new uint16[dim];

    uint i, j;

    uint16 * D_ptr = &D[0];

    for (i = 0; i < M_; i++)
        D_ptr[i] = i;

    for (j = 1; j < N_; j++)
        D_ptr[j * M_] = j;

    ////////////////

    int M_Tiles = ceil((float) M / Worker::TILE_WIDTH);
    int N_Tiles = ceil((float) N / Worker::TILE_WIDTH);

    const int DStart = -M_Tiles + 1;
    const int DFinish = N_Tiles;

    const int TOTAL_TILES = M_Tiles * N_Tiles;
    float t1, t2;


    std::atomic_bool * tileComputed = new std::atomic_bool[TOTAL_TILES];
    for(i = 0; i < TOTAL_TILES; i++){
        tileComputed[i] = false;
    }

    barrier threadBarrier(Worker::MAX_THREAD_COUNT + 1);
    ConcurrentQueue<ind> indexQueue;
    std::vector<std::thread> threadVec(Worker::MAX_THREAD_COUNT);

    for (int d = DStart; d < DFinish; d++) {
        const int iMax = std::min(M_Tiles + d, N_Tiles);
        const int iMin = std::max(d, 0);
        for (i = iMin; i < iMax; i++) {
            j = M_Tiles - i + d - 1;
            indexQueue.push_unsafe(ind(i,j));
        }
    }

    //poison pills
    for(i = 0; i < Worker::MAX_THREAD_COUNT; i++){
        indexQueue.push_unsafe(ind(-1,-1));
    }
    t1 = omp_get_wtime();
    //launch threads asynchronously
    for(i = 0; i < Worker::MAX_THREAD_COUNT; i++){
        threadVec[i] = std::thread(Worker(x,y,D, &threadBarrier, indexQueue, tileComputed));
        threadVec[i].detach();
    }
    threadBarrier.count_down_and_wait(); // wait threads for the completion
    t2 = omp_get_wtime();
    //printf("Pure computing time: %f\n", t2-t1);

    int distance = D[N * M_ + M];

    delete[] D;
    delete[] tileComputed;
    return distance;
}

int editDistanceCPPT2(const char *x, const char *y) {
    const int M = strlen(x);
    const int N = strlen(y);

    const uint M_ = M + 1;
    const uint N_ = N + 1;

    const uint dim = (uint)M_*(uint)N_;
    uint16 * D = new uint16[dim];

    uint i, j;

    for (i = 0; i < M_; i++)
        D[i] = i;

    for (j = 1; j < N_; j++)
        D[j * M_] = j;

    ////////////////

    int M_Tiles = ceil((float) M / Worker2::TILE_WIDTH);
    int N_Tiles = ceil((float) N / Worker2::TILE_WIDTH);

    const int DStart = -M_Tiles + 1;
    const int DFinish = N_Tiles;

    const int TOTAL_TILES = M_Tiles * N_Tiles;
    float t1, t2;

    barrier threadBarrier(Worker2::MAX_THREAD_COUNT + 1);
    ConcurrentQueue<ind> indexQueue;
    std::vector<std::thread> threadVec(Worker2::MAX_THREAD_COUNT);

    for(i = 0; i < Worker2::MAX_THREAD_COUNT; i++){
        threadVec[i] = std::thread(Worker2(x,y,D, &threadBarrier, indexQueue));
        threadVec[i].detach();
    }


    for (int d = DStart; d < DFinish; d++) {
        const int iMax = std::min(M_Tiles + d, N_Tiles); //handling the case when we are filling the left up corner of the matrix
        const int iMin = std::max(d, 0); //handling the case when we are filling the right low corner of the matrix

        for (i = iMin; i < iMax; i++) {
            j = M_Tiles - i + d - 1;
            indexQueue.push(ind(i,j));
        }
         for (i = 0; i < Worker2::MAX_THREAD_COUNT; i++) {
                indexQueue.push(ind(-2, -2));
         }
         threadBarrier.count_down_and_wait();
    }

    //poison pills
    for(i = 0; i < Worker2::MAX_THREAD_COUNT; i++){
        indexQueue.push(ind(-1,-1));
    }


    threadBarrier.count_down_and_wait();
    int distance = D[N * M_ + M];

    delete[] D;
    return distance;
}
int editDistanceCPPT3(const char *x, const char *y) {
    const int M = strlen(x);
    const int N = strlen(y);

    const uint M_ = M + 1;
    const uint N_ = N + 1;

    const uint dim = (uint)M_*(uint)N_;
    uint16 * D = new uint16[dim];

    uint i, j;

    for (i = 0; i < M_; i++)
        D[i] = i;

    for (j = 1; j < N_; j++)
        D[j * M_] = j;

    ////////////////

    barrier threadBarrier(Worker3::MAX_THREAD_COUNT + 1);
    std::vector<std::thread> threadVec(Worker3::MAX_THREAD_COUNT);

    for(i = 0; i < Worker3::MAX_THREAD_COUNT; i++){
        threadVec[i] = std::thread(Worker3(x,y,D, &threadBarrier));
        threadVec[i].detach();
    }

    threadBarrier.count_down_and_wait();
    Worker3::nThreads = 0;
    int distance = D[N * M_ + M];
    delete[] D;
    return distance;
}

int main() {

    std::cout << "uint16 is: " << sizeof(uint16)*8 << " bit" << std::endl;

    const int N = 50000+1;
    std::cout << "N is: " << N-1 << std::endl;
    char *A = new char[N];
    char *B = new char[N];

    char LUT[4] = {'A', 'T', 'C', 'G'};

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, 3);
    std::uniform_real_distribution<double> uniformProb(0.0, 1.0);

    const float pChange = 0.5;
    int priorDistance = 0;

    for (int i = 0; i < N-1; i++) {
        int dice = distribution(generator);
        B[i] = A[i] = LUT[dice];

        if (uniformProb(generator) < pChange) {
            int dice2 = distribution(generator);
            if (dice2 != dice) {
                B[i] = LUT[dice2];
                ++priorDistance;
            }
        }

    }

    A[N - 1] = B[N - 1] = '\0';


    if (strlen(A) < strlen(B)) {
        std::swap(A, B);
    }


    printf("computing edit distance\n");
    int d;
    double t1, t2;
    double ST_Time = 0.0, CPPT_Time = 0.0, CPPT2_Time = 0.0, CPPT3_Time = 0.0, OMP_Time = 0.0;
    int iterations = 25;
    //ST_Time = 0.0014*iterations;
for(int i = 0; i < iterations; i++) {
    printf("Iteration n: %i\n", i+1);

    t1 = omp_get_wtime();
    d = editDistanceST(A, B);
    t2 = omp_get_wtime();
    ST_Time += t2 - t1;
    std::cout << "The edit distance is: " << d << "; Computed in (ST): " << t2 - t1 << std::endl;

    t1 = omp_get_wtime();
    d = editDistanceCPPT(A, B);
    t2 = omp_get_wtime();
    CPPT_Time += t2 - t1;
    std::cout << "The edit distance is: " << d << "; Computed in (Dynamic Scheduling): " << t2 - t1 << std::endl;

    t1 = omp_get_wtime();
    d = editDistanceCPPT2(A, B);
    t2 = omp_get_wtime();
    CPPT2_Time += t2 - t1;

    std::cout << "The edit distance is: " << d << "; Computed in (Dynamic Barrier): " << t2 - t1 << std::endl;

    t1 = omp_get_wtime();
    d = editDistanceCPPT3(A, B);
    t2 = omp_get_wtime();
    CPPT3_Time += t2 - t1;

    std::cout << "The edit distance is: " << d << "; Computed in (BlockCyclic): " << t2 - t1 << std::endl;
    t1 = omp_get_wtime();
    d = editDistanceOMP(A, B);
    t2 = omp_get_wtime();
    OMP_Time += t2 - t1;
    std::cout << "The edit distance is: " << d << "; Computed in (OMP): " << t2 - t1 << std::endl;
 }
    std::cout << "ST time " << ST_Time/iterations << std::endl;
    std::cout << "CPPT time " << CPPT_Time/iterations << " SpeedUp: " << ST_Time/CPPT_Time << std::endl;
    std::cout << "CPPT2 time " << CPPT2_Time/iterations << " SpeedUp: " << ST_Time/CPPT2_Time << std::endl;
    std::cout << "CPPT3 time " << CPPT3_Time/iterations << " SpeedUp: " << ST_Time/CPPT3_Time << std::endl;
    std::cout << "OMP time " << OMP_Time/iterations << " SpeedUp: " << ST_Time/OMP_Time <<std::endl;


    return 0;
}