//
// Created by fede on 21/01/19.
//

#ifndef LEVENSHTEIN_WORKER2_H
#define LEVENSHTEIN_WORKER2_H
#include <iostream>
#include <string>
#include <vector>
#include <string.h>
#include <random>
#include <omp.h>
#include <unistd.h>
#include <thread>
#include <atomic>
#include "BoostBarrier.h"
#include "ind.h"
#include "uint.h"
#include "ConcurrentQueue.h"

class Worker2 {

public:
    void operator()() {

        ind ind1;
        while(1) {

                indexQueue.pop(ind1);
                if(ind1.i == -1){
                    //printf("Thread %i is quitting...\n", tid);
                    break;
                }
                if(ind1.i == -2){
                    thread_barrier->count_down_and_wait();
                    continue;
                }
                I = ind1.i * TILE_WIDTH + 1;
                J = ind1.j * TILE_WIDTH + 1;

                computeSubMatrix();
        }

        thread_barrier->count_down_and_wait();
    }

    Worker2(const char *x, const char *y, uint16 *D, barrier *b, ConcurrentQueue<ind>& indexQueue) :
            D(D), x(x), y(y), thread_barrier(b), indexQueue(indexQueue) {
        M = strlen(x);
        N = strlen(y);
        M_Tiles = ceil((float) M / Worker2::TILE_WIDTH);
        N_Tiles = ceil((float) N / Worker2::TILE_WIDTH);
        tid = nThreads++;
    }

    static const int TILE_WIDTH = 512;
    static const int MAX_THREAD_COUNT = 12;

private:

    void computeSubMatrix() {
        uint M_ = M + 1;
        uint N_ = N + 1;
        for (uint i = I; i < N_ && i < I + TILE_WIDTH; i++) {
            for (uint j = J; j < M_ && j < J + TILE_WIDTH; j++) {
                if (x[i - 1] != y[j - 1]) {
                    uint16 k = (uint16)minimum_(D[i * M_ + j - 1], //insertion
                                     D[(i - 1) * M_ + j], //insertion
                                     D[(i - 1) * M_ + j - 1]); //substitution
                    D[i * M_ + j] = k + 1;
                } else {
                    D[i * M_ + j] = D[(i - 1) * M_ + j - 1];
                }
            }
        }
    }

    int minimum_(const int a, const int b, const int c) {
        return std::min(std::min(a, b), c);
    }


    uint16 *D;
    int tx, ty;
    uint M, N;
    uint I, J;
    int M_Tiles, N_Tiles;
    const char *x, *y;
    barrier *thread_barrier;
    ConcurrentQueue<ind>& indexQueue;
    int tid;
    static int nThreads;
};




#endif //LEVENSHTEIN_WORKER_H
