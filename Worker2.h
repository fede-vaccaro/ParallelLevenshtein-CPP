//
// Created by fede on 21/01/19.
//

#ifndef LEVENSHTEIN_WORKER2_H
#define LEVENSHTEIN_WORKER2_H

#include "BoostBarrier.h"
#include <vector>
#include "ind.h"
#include <string.h>
#include <iostream>
#include <cmath>

class Worker2 {

public:
    Worker2(const char *x, const char *y, int *D, barrier *b, ind startStopA, ind startStopB) :
            D(D), x(x), y(y), b(b), startStopA(startStopA), startStopB(startStopB){
        M = strlen(x);
        N = strlen(y);
        M_Tiles = ceil((float) M / Worker2::TILE_WIDTH);
        N_Tiles = ceil((float) N / Worker2::TILE_WIDTH);
        tid = threadCount++;
    }

    void operator()(){
        for(int i = 0; i < N_Tiles; ++i ){
            for(int j = startStopA.i; j < startStopA.j; ++j){
                // (i,j) is the submatrix index
                computeSubMatrixA(i, j);
            }
            sharedBarrier.count_down_and_wait();
        }

        for(int j = 0; j < M_Tiles; ++j){
            for(int i = startStopB.i; i < startStopB.j; ++i){
                computeSubMatrixB(i, j);
            }
            sharedBarrier.count_down_and_wait();
        }
        b->count_down_and_wait();
    }

    static const int TILE_WIDTH = 128;
    static const int MAX_THREAD_COUNT = 12;

private:

    void computeSubMatrixA(int I, int J) {
        int M_ = M + 1;
        int N_ = N + 1;

        // now (i,j) is the element index, of the (I,J) sub matrix
        int iStart = I*TILE_WIDTH + 1;
        int jStart = J*TILE_WIDTH + 1;



        for (int i = iStart; i < N_ && i < iStart + TILE_WIDTH; i++) {
            for (int j = jStart; j < M_ && j < jStart + TILE_WIDTH; j++) {
                if (x[i - 1] != y[j - 1]) {
                    int k = std::min(D[i * M_ + j - 1],
                            D[(i - 1) * M_ + j - 1]);
                    D[i * M_ + j] = k + 1;
                } else {
                    D[i * M_ + j] = D[(i - 1) * M_ + j - 1];
                }
            }
        }
    }

    void computeSubMatrixB(int I, int J) {
        int M_ = M + 1;
        int N_ = N + 1;

        int iStart = I*TILE_WIDTH + 1;
        int jStart = J*TILE_WIDTH + 1;

        for (int i = iStart; i < N_ && i < iStart + TILE_WIDTH; i++) {
            for (int j = jStart; j < M_ && j < jStart + TILE_WIDTH; j++) {
                if (x[i - 1] != y[j - 1]) {
                    int k = std::min(D[i * M_ + j],
                                     D[i * M_ + j - 1]);
                    D[i * M_ + j] = k + 1;
                } //else {
                //    D[i * M_ + j] = D[(i - 1) * M_ + j - 1];
                //}
            }
        }
    }
    int *D;
    int tx, ty;
    int M, N;
    int tid;
    ind startStopA, startStopB;
    int M_Tiles, N_Tiles;
    barrier * b;
    const char *x, *y;
    static barrier sharedBarrier;
    static int threadCount;
};


#endif //LEVENSHTEIN_WORKER2_H
