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
#include "MyBarrier.h"
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

int editDistanceST(const char *x, const char *y) {
    const int M = strlen(x);
    const int N = strlen(y);

    int distance;

    int **D = new int *[M + 1];

    int i, j;

    for (i = 0; i <= M; i++)
        D[i] = new int[N + 1];

    for (i = 0; i <= M; i++)
        D[i][0] = i;

    for (j = 1; j <= N; j++)
        D[0][j] = j;

    for (i = 1; i <= M; i++) {
        for (j = 1; j <= N; j++) {
            if (x[i - 1] != y[j - 1]) {
                int k = minimum_(D[i][j - 1], //insertion
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

struct Worker {

public:
    void operator()() {
        acquirePermission();
        //printf("Hello from thread %i, %i!\n", tx, ty);
        computeSubMatrix();
        releasePermission();
        thread_barrier->count_down_and_wait();

        //printf("Quitting from thread %i, %i!\n", tx, ty);

    }

    Worker(int i, int j, const char *x, const char *y, int *D, barrier *b, int tx, int ty,
            std::vector<std::mutex>& mutexVec, std::vector<std::condition_variable>& condVec, std::atomic_bool* tileComputed) :
    I(i), J(j), D(D), x(x), y(y), tx(tx), ty(ty), mutexVec(mutexVec), condVec(condVec), tileComputed(tileComputed), thread_barrier(b) {
        M = strlen(x);
        N = strlen(y);
    }

    static const int TILE_WIDTH = 512;

private:

    void computeSubMatrix() {
        int M_ = M + 1;
        int N_ = N + 1;
        for (int i = I; i < N_ && i < I + TILE_WIDTH; i++) {
            for (int j = J; j < M_ && j < J + TILE_WIDTH; j++) {
                if (x[i - 1] != y[j - 1]) {
                    int k = minimum_(D[i * M_ + j - 1], //insertion
                                     D[(i - 1) * M_ + j], //insertion
                                     D[(i - 1) * M_ + j - 1]); //substitution
                    D[i * M_ + j] = k + 1;
                } else {
                    D[i * M_ + j] = D[(i - 1) * M_ + j - 1];
                }
            }
        }
    }

    void acquirePermission(){
        int M_Tiles = ceil((float) M / Worker::TILE_WIDTH);
        int N_Tiles = ceil((float) N / Worker::TILE_WIDTH);

        std::unique_lock<std::mutex> lock(mutexVec[tx*M_Tiles + ty]);

        if(tx > 0 && ty > 0) {
            while (!(tileComputed[(tx - 1) * M_Tiles + ty] &&
                   tileComputed[tx * M_Tiles + ty - 1] &&
                   tileComputed[(tx - 1) * M_Tiles + ty - 1])) {
                condVec[tx * M_Tiles + ty].wait(lock);
            }
        }else if(tx > 0){
            while (!tileComputed[(tx - 1) * M_Tiles + ty]) {
                condVec[tx * M_Tiles + ty].wait(lock);
            }
        }else if(ty > 0){
            while (!tileComputed[tx * M_Tiles + ty - 1]) {
                condVec[tx * M_Tiles + ty].wait(lock);
            }
        }

    }

    void releasePermission(){
        int M_Tiles = ceil((float) M / Worker::TILE_WIDTH);
        int N_Tiles = ceil((float) N / Worker::TILE_WIDTH);

        tileComputed[tx*M_Tiles + ty] = true;

        if(tx < N_Tiles - 1) {
            condVec[(tx+1) * M_Tiles + ty].notify_one();
        }
        if(ty < M_Tiles - 1) {
            condVec[tx * M_Tiles + ty + 1].notify_one();
        }
        if(tx < N_Tiles - 1 && ty < M_Tiles - 1){
            condVec[(tx+1)*M_Tiles + ty + 1].notify_one();
        }
    }

    int *D;
    int tx, ty;
    int M, N;
    int I, J;
    const char *x, *y;
    std::vector<std::mutex>& mutexVec;
    std::vector<std::condition_variable>& condVec;
    std::atomic_bool* tileComputed;
    barrier *thread_barrier;
};

class Dummy {
public:
    void operator()() {
        printf("Dummy class\n");
    }
};

typedef struct ind{
    ind(int i, int j): i(i), j(j){}
    int i, j;
} ind;

int editDistanceST2(const char *x, const char *y) {
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

    int M_Tiles = ceil((float) M / Worker::TILE_WIDTH);
    int N_Tiles = ceil((float) N / Worker::TILE_WIDTH);

    const int DStart = -M_Tiles+1;
    const int DFinish = N_Tiles;
    const int TILE_WIDTH = Worker::TILE_WIDTH;

    std::vector<std::thread> threadVec;
    threadVec.resize(M_Tiles * N_Tiles);

    float t2, t1;

    t2 = omp_get_wtime();
    std::vector<std::mutex> mutexVec(M_Tiles*N_Tiles);
    std::vector<std::condition_variable> condVec(M_Tiles*N_Tiles);
    std::atomic_bool * tileComputed = new std::atomic_bool[M_Tiles*N_Tiles];
    for(i = 0; i < M_Tiles*N_Tiles; i++){
        tileComputed[i] = false;
    }
    int TOTAL_THREAD = M_Tiles * N_Tiles;

    boost::lockfree::queue<struct ind> indexQueue;
    for (int d = DStart; d < DFinish; d++) {
        const int iMax = std::min(M_Tiles + d, N_Tiles);
        const int iMin = std::max(d, 0);
        int nThread = iMax - iMin;
        for (i = iMin; i < iMax; i++) {
            j = M_Tiles - i + d - 1;
            indexQueue.push(ind(i,j));
            printf("Pushed: (%i,%i)\n",i,j);
        }
    }
    printf("TOTAL_THREAD: %i\n", TOTAL_THREAD);
    barrier b(TOTAL_THREAD+1);
    t1 = omp_get_wtime();
    printf("Time to initialize vectors: %f\n", t2 - t1);

    t1 = omp_get_wtime();
    /*for (int d = DStart; d < DFinish; d++) {
        const int iMax = std::min(M_Tiles + d, N_Tiles);
        const int iMin = std::max(d, 0);
        int nThread = iMax - iMin;
        for (i = iMin; i < iMax; i++) {
            j = M_Tiles - i + d - 1;
            //printf("i: %i, j: %i, nThread: %i\n",i,j,nThread);
            std::thread t = std::thread(Worker(i*TILE_WIDTH + 1, j*TILE_WIDTH + 1,x,y,&D[0], &b, i, j, mutexVec, condVec, tileComputed));
            t.detach();
            /*if (i < 0 || i >= N_Tiles)
                printf("ERROR i out of bound: %i\n", i);
            if (j < 0 || j >= M_Tiles)
                printf("ERROR j out of bound %i\n", j);*/

        }
        /*for(i = iMin; i < iMax; i++){
            j = M_Tiles - i + d - 1;
            if(threadVec[i*M_Tiles + j].joinable())
                threadVec[i*M_Tiles + j].join();
        }
        //b.count_down_and_wait();
        //printf("Cycle %i finished\n", d - DStart);
    }*/
    t2 = omp_get_wtime();
    printf("Time to launch threads: %f\n", t2-t1);
/*
    std::thread tx = std::thread();

    printf("Number of iterations: %i\n", DFinish - DStart);
    for (int d = DStart; d < DFinish; d++) {
        const int iMax = std::min(M_Tiles + d - 1, N_Tiles);
        const int iMin = std::max(d, 1);
        for (i = iMin; i < iMax; i++) {
            j = M_Tiles - i + d - 1;
            int nThread = iMax - iMin;
            my_barrier threadBarrier(nThread + 1);


            threadBarrier.wait();

            if (i < 0 || i >= N_Tiles)
                printf("ERROR i out of bound: %i\n", i);
            if (j < 0 || j >= M_Tiles)
                printf("ERROR j out of bound %i\n", j);

        }
    }*/
/*
    //boost::barrier b(2);
    my_barrier thread_barrier(2);
    Worker w1 = Worker(1, 1, x, y, &D[0], &thread_barrier);

    std::thread th1(w1);

    th1.detach();
    thread_barrier.wait();
    printf("T0 finished\n");

    my_barrier thread_barrier2(3);
    Worker w2 = Worker(1, 5, x, y, &D[0], &thread_barrier2);
    Worker w3 = Worker(5, 1, x, y, &D[0], &thread_barrier2);
    std::thread th2(w2), th3(w3);
    th2.detach();
    th3.detach();
    thread_barrier2.wait();
    printf("T1 and T2 finished\n");


    Worker w4 = Worker(5, 5, x, y, &D[0], &thread_barrier);
    std::thread th4(w4);
    th4.detach();
    thread_barrier.wait();
    printf("T3 finished\n");
*/

    b.count_down_and_wait();
    int distance = D[N * M_ + M];

    /*
    for (i = 0; i < N_; i++) {
        for (j = 0; j < M_; j++) {
            printf("%i, ", D[i * M_ + j]);
        }
        std::cout << std::endl;
    }*/

    //for (i = 0; i <= M; i++)
    //    delete D[i];
    delete tileComputed;
    return distance;
}


int main() {

    const int N = 40000;
    std::cout << "N is: " << N << std::endl;
    char* A = new char[N];
    char* B = new char[N - 20];

    for(int i = 0; i < N; i++) {

        A[i] = 'a';

    }

    for(int i = 0; i < N - 20; i++){
        B[i] = 'a';
    }

    B[(N - 20)/2] = 'b';
        /*
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> dis(0, 255);
        char value = (char)(dis(gen));
        A[i] = value;
        B[i] = value;
        if(B[i] % 5){
            B[i] = B[i] + 1;
        }*/



    A[N - 1] = B[N - 20 - 1] = '\0';




    //const char *A = "aaabbaaa";
    //const char *B = "aaabcaaa";

    float t1 = omp_get_wtime();

    if (strlen(A) < strlen(B)) {
        std::swap(A, B);
    }


    printf("computing edit distance\n");
    int d = editDistanceST(A, B);
    float t2 = omp_get_wtime();
    std::cout << "The edit distance is: " << d << "; Computed in: " << t2 - t1 << std::endl;
    return 0;
}