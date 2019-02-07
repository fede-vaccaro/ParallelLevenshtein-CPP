//
// Created by fede on 21/01/19.
//

#ifndef LEVENSHTEIN_WORKER_H
#define LEVENSHTEIN_WORKER_H
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

class LockBased {

public:
    void operator()() {
        ind ind1;
        while(1) {
                indexQueue.pop(ind1);
                if(ind1.i == -1){
                    break;
                }
                I = ind1.i * TILE_WIDTH + 1;
                J = ind1.j * TILE_WIDTH + 1;
                tx = ind1.i;
                ty = ind1.j;

                checkPermission(); // check if the assigned (i,j) tile is ready to be computed, checking if (i-1,j), (i,j-1) and (i-1,j-1) are been computed
                computeSubMatrix(); // compute the (i,j) submatrix
                releasePermission(); // set the tile (i,j) as computed
        }
        thread_barrier->count_down_and_wait(); // to synchronize each thread after the completion
    }

    LockBased(const char *x, const char *y, uint16 *D, barrier *b, ConcurrentQueue<ind>& indexQueue, std::atomic_bool* tileComputed) :
            D(D), x(x), y(y), tileComputed(tileComputed), thread_barrier(b), indexQueue(indexQueue) {
        M = strlen(x);
        N = strlen(y);
        M_Tiles = ceil((float) M / LockBased::TILE_WIDTH);
        N_Tiles = ceil((float) N / LockBased::TILE_WIDTH);
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
                    uint16 k = (uint16)minimum_(D[i * M_ + j - 1],
                                     D[(i - 1) * M_ + j],
                                     D[(i - 1) * M_ + j - 1]);
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

    void checkPermission(){

        if(tx > 0 && ty > 0) {
            while (!(tileComputed[(tx - 1) * M_Tiles + ty] &&
                     tileComputed[tx * M_Tiles + ty - 1] &&
                     tileComputed[(tx - 1) * M_Tiles + ty - 1])) {
            }
        }else if(tx > 0){
            while (!tileComputed[(tx - 1) * M_Tiles + ty]) {
            }
        }else if(ty > 0){
            while (!tileComputed[tx * M_Tiles + ty - 1]) {
            }
        }

    }

    void releasePermission(){
        tileComputed[tx*M_Tiles + ty] = true;
    }


    uint16 *D;
    int tx, ty;
    uint M, N;
    uint I, J;
    int M_Tiles, N_Tiles;
    const char *x, *y;
    std::atomic_bool* tileComputed;
    barrier *thread_barrier;
    ConcurrentQueue<ind>& indexQueue;
    int tid;
    static int nThreads;
};


#endif //LEVENSHTEIN_WORKER_H
