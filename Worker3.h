//
// Created by fede on 21/01/19.
//

#ifndef LEVENSHTEIN_Worker3_H
#define LEVENSHTEIN_Worker3_H

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

class Worker3 {

public:
    void operator()() {

        const int DStart = -M_Tiles + 1;
        const int DFinish = N_Tiles;

        for (int d = DStart; d < DFinish; d++) {
            const uint16 iMax = std::min(M_Tiles + d,
                                      N_Tiles); //handling the case when we are filling the left up corner of the matrix
            const uint16 iMin = std::max(d, 0); //handling the case when we are filling the right low corner of the matrix

            uint16 diagLength = iMax - iMin;
            if (diagLength > MAX_THREAD_COUNT) {
                const uint16 tilePerThread = diagLength / MAX_THREAD_COUNT;
                lower = tid * tilePerThread + iMin;
                upper = lower + tilePerThread;

            } else if (tid < diagLength) {
                lower = tid + iMin;
                upper = lower + 1;
            } else {
                lower = upper = -1;
            }

            uint16 remainder = diagLength % MAX_THREAD_COUNT;

            for (int i = lower; i < upper; i++) {

                uint16 j = M_Tiles - i + d - 1;
                I = i * TILE_WIDTH + 1;
                J = j * TILE_WIDTH + 1;

                computeSubMatrix();
            }

            if(remainder && tid < remainder){
                int i = iMax - tid - 1;
                uint16 j = M_Tiles - i + d - 1;
                I = i * TILE_WIDTH + 1;
                J = j * TILE_WIDTH + 1;

                computeSubMatrix();
            }

            privateBarrier.count_down_and_wait();
        }
        thread_barrier->count_down_and_wait();
    }


    Worker3(const char *x, const char *y, uint16 *D, barrier *b) :
            D(D), x(x), y(y), thread_barrier(b) {
        M = strlen(x);
        N = strlen(y);
        M_Tiles = ceil((float) M / Worker3::TILE_WIDTH);
        N_Tiles = ceil((float) N / Worker3::TILE_WIDTH);
        tid = nThreads++;
        lower = upper = 0;
    }

    void setLUd(uint lower, uint upper, uint d) {
        this->lower = lower;
        this->upper = upper;
    }

    static const int TILE_WIDTH = 512;
    static const int MAX_THREAD_COUNT = 12;
    static barrier privateBarrier;
    static int nThreads;

private:

    void computeSubMatrix() {
        uint M_ = M + 1;
        uint N_ = N + 1;
        for (uint i = I; i < N_ && i < I + TILE_WIDTH; i++) {
            for (uint j = J; j < M_ && j < J + TILE_WIDTH; j++) {
                if (x[i - 1] != y[j - 1]) {
                    uint16 k = (uint16) minimum_(D[i * M_ + j - 1], //insertion
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

    uint lower, upper;
    uint16 *D;
    int tx, ty;
    uint M, N;
    uint I, J;
    int M_Tiles, N_Tiles;
    const char *x, *y;
    barrier *thread_barrier;
    int tid;
};


#endif //LEVENSHTEIN_WORKER_H
