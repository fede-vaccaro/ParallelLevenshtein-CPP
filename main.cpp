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
#include "ctpl_stl.h"
#include "BoostBarrier.h"
#include "ind.h"
#include "Worker.h"
#include "Worker2.h"
#include "Worker3.h"
#include <boost/lockfree/queue.hpp>

int minimum(int a, int b, int c) {
    int min = a;

    if (b < min) min = b;
    if (c < min) min = c;

    return min;
}

int minimum_(const int a, const int b, const int c) {
    return std::min(std::min(a, b), c);
}

typedef unsigned short int uint16;

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

    delete D;
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

const int MAX_NUM_THREADS = 8;

const int TW = 512;

void computeSubMatrix(int I, int J, const int M, const int N, const char *x, const char *y, int *D) {
    int M_ = M + 1;
    int N_ = N + 1;
    I = I * TW + 1;
    J = J * TW + 1;
    for (int i = I; i < N_ && i < I + TW; i++) {
        for (int j = J; j < M_ && j < J + TW; j++) {
            if (x[i - 1] != y[j - 1]) {
                int k = minimum(D[i * M_ + j - 1], //insertion
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
    const int M = strlen(x);
    const int N = strlen(y);

    const int M_ = M + 1;
    const int N_ = N + 1;

    //int ** D = new int *[M+1];
    std::vector<int> D;
    D.resize(M_ * N_);

    int i, j;

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

#pragma omp parallel num_threads(8)
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

    return distance;
}

int editDistanceCPPT(const char *x, const char *y) {
    const int M = strlen(x);
    const int N = strlen(y);

    const int M_ = M + 1;
    const int N_ = N + 1;

    std::vector<int> D;
    D.resize(M_ * N_);

    int i, j;

    for (i = 0; i < M_; i++)
        D[i] = i;

    for (j = 1; j < N_; j++)
        D[j * M_] = j;

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
    boost::lockfree::queue<ind> indexQueue(TOTAL_TILES + Worker::MAX_THREAD_COUNT);

    std::vector<std::thread> threadVec(Worker::MAX_THREAD_COUNT);


    for (int d = DStart; d < DFinish; d++) {
        const int iMax = std::min(M_Tiles + d, N_Tiles);
        const int iMin = std::max(d, 0);
        for (i = iMin; i < iMax; i++) {
            j = M_Tiles - i + d - 1;
            indexQueue.push(ind(i,j));
        }
    }

    for(i = 0; i < Worker::MAX_THREAD_COUNT; i++){
        indexQueue.push(ind(-1,-1));
    }

    for(i = 0; i < Worker::MAX_THREAD_COUNT; i++){
        threadVec[i] = std::thread(Worker(x,y,&D[0], &threadBarrier, indexQueue, tileComputed));
        threadVec[i].detach();
    }

    threadBarrier.count_down_and_wait();
    int distance = D[N * M_ + M];

    delete tileComputed;
    return distance;
}

/*
int editDistanceMT2(const char *x, const char *y) {
    const int M = strlen(x);
    const int N = strlen(y);

    const int M_ = M + 1;
    const int N_ = N + 1;

    //int ** D = new int *[M+1];
    std::vector<int> D;
    D.resize(M_ * N_);

    int i, j;

    for (i = 0; i < M_; i++)
        D[i] = i;

    for (j = 1; j < N_; j++)
        D[j * M_] = j;

    ////////////////

    int M_Tiles = ceil((float) M / Worker2::TILE_WIDTH);
    int N_Tiles = ceil((float) N / Worker2::TILE_WIDTH);

    const int TOTAL_TILES = M_Tiles * N_Tiles;

    barrier b(Worker2::MAX_THREAD_COUNT + 1);
    std::vector<std::thread> threadVec;
    const int tilesPerThreadA = ceil((float) M_Tiles / Worker2::MAX_THREAD_COUNT);
    const int tilesPerThreadB = ceil((float) N_Tiles / Worker2::MAX_THREAD_COUNT);
    for (i = 0; i < Worker2::MAX_THREAD_COUNT; i++) {

        ind startStopA = ind(i * tilesPerThreadA, (i + 1) * tilesPerThreadA);
        ind startStopB = ind(i * tilesPerThreadB, (i + 1) * tilesPerThreadB);

        if (startStopA.i > M_Tiles || startStopB.i > N_Tiles)
            break;

        threadVec.push_back(std::thread(
                Worker2(x, y, &D[0], &b,
                        startStopA, startStopB)
        ));
        threadVec[i].detach();
    }


    b.count_down_and_wait();
    int distance = D[N * M_ + M];

    //for (i = 0; i <= M; i++)
    //    delete D[i];

    return distance;
}

*/
/*
int editDistanceMT3(const char *x, const char *y) {
    const int M = strlen(x);
    const int N = strlen(y);

    const int M_ = M + 1;
    const int N_ = N + 1;

    //int ** D = new int *[M+1];
    std::vector<int> D;
    D.resize(M_ * N_);

    int i, j;

    for (i = 0; i < M_; i++)
        D[i] = i;

    for (j = 1; j < N_; j++)
        D[j * M_] = j;

    ////////////////

    int M_Tiles = ceil((float) M / Worker2::TILE_WIDTH);
    int N_Tiles = ceil((float) N / Worker2::TILE_WIDTH);

    const int TOTAL_TILES = M_Tiles * N_Tiles;

    barrier b(Worker3::MAX_THREAD_COUNT + 1);
    std::vector<std::thread> threadVec;
    const int tilesPerThreadA = ceil((float) M_Tiles / Worker2::MAX_THREAD_COUNT);
    const int tilesPerThreadB = ceil((float) N_Tiles / Worker2::MAX_THREAD_COUNT);
    for (i = 0; i < Worker3::MAX_THREAD_COUNT; i++) {

        ind startStopA = ind(i * tilesPerThreadA, (i + 1) * tilesPerThreadA);
        ind startStopB = ind(i * tilesPerThreadB, (i + 1) * tilesPerThreadB);

        if (startStopA.i > M_Tiles || startStopB.i > N_Tiles)
            break;

        threadVec.push_back(std::thread(
                Worker2(x, y, &D[0], &b,
                        startStopA, startStopB)
        ));
        threadVec[i].detach();
    }


    b.count_down_and_wait();
    int distance = D[N * M_ + M];

    //printMatrix(&D[0], M_, N_);

    //for (i = 0; i <= M; i++)
    //    delete D[i];

    return distance;
}
*/

int main() {

    const int N = 40000;
    std::cout << "N is: " << N << std::endl;
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

    //std::cout << "Edit distance of the generated string: " << priorDistance << std::endl;


    //const char *A = "aaabbaaa";
    //const char *B = "aaabcaaa";


    if (strlen(A) < strlen(B)) {
        std::swap(A, B);
    }


    printf("computing edit distance\n");
    int d;
    double t1, t2;

    t1 = omp_get_wtime();
    d = editDistanceST(A, B);
    t2 = omp_get_wtime();
    std::cout << "The edit distance is: " << d << "; Computed in (ST): " << t2 - t1 << std::endl;

    t1 = omp_get_wtime();
    d = editDistanceCPPT(A, B);
    t2 = omp_get_wtime();
    std::cout << "The edit distance is: " << d << "; Computed in (CPPT): " << t2 - t1 << std::endl;

    t1 = omp_get_wtime();
    d = editDistanceOMP(A, B);
    t2 = omp_get_wtime();
    std::cout << "The edit distance is: " << d << "; Computed in (OMP): " << t2 - t1 << std::endl;

    return 0;
}