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
#include <mutex>
#include <condition_variable>
#include <atomic>
#include "BoostBarrier.h"
#include <boost/lockfree/queue.hpp>
#include "ind.h"
#include "uint.h"
#include "ConcurrentQueue.h"

class Worker {

public:
    void operator()() {
        int i = 0;
        //while(!indexQueue.empty()){// && !tileComputed[(M_Tiles-1)*M_Tiles + N_Tiles - 1]) {
        ind ind1;
        float cumulativeT = 0.0;
        float t2, t1;
        while(1) {
            /*while(indexQueue.empty()){
                if(tileComputed[(M_Tiles-1)*M_Tiles + N_Tiles - 1])
                    break;
            }*/

            //if(!tileComputed[(M_Tiles-1)*M_Tiles + N_Tiles - 1]) {
                t1 = omp_get_wtime();
                /*while(indexQueue.empty()){

                }*/
                indexQueue.pop(ind1);
                if(ind1.i == -1){
                    //printf("Thread %i is quitting...\n", tid);
                    break;
                }
                ++i;
                I = ind1.i * TILE_WIDTH + 1;
                J = ind1.j * TILE_WIDTH + 1;
                tx = ind1.i;
                ty = ind1.j;

                checkPermission();
                //t2 = omp_get_wtime();
                //cumulativeT += (t2 - t1);
                computeSubMatrix();
                //t1 = omp_get_wtime();
                releasePermission();
                //t2 = omp_get_wtime();
                //cumulativeT += (t2 - t1);
        }
        //printf("Total overhead by thread %i in %i iterations: %f\n", tid, i, cumulativeT);
        thread_barrier->count_down_and_wait();
    }

    Worker(const char *x, const char *y, uint16 *D, barrier *b, ConcurrentQueue<ind>& indexQueue,
            /*std::vector<std::mutex>& mutexVec, std::vector<std::condition_variable>& condVec,*/ std::atomic_bool* tileComputed) :
            D(D), x(x), y(y), /*mutexVec(mutexVec), condVec(condVec),*/ tileComputed(tileComputed), thread_barrier(b), indexQueue(indexQueue) {
        M = strlen(x);
        N = strlen(y);
        M_Tiles = ceil((float) M / Worker::TILE_WIDTH);
        N_Tiles = ceil((float) N / Worker::TILE_WIDTH);
        tid = nThreads++;
    }

    static const int TILE_WIDTH = 512;
    static const int MAX_THREAD_COUNT = 11;

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

    void checkPermission(){
        //std::unique_lock<std::mutex> lock(mutexVec[tx*M_Tiles + ty]);

        if(tx > 0 && ty > 0) {
            while (!(tileComputed[(tx - 1) * M_Tiles + ty] &&
                     tileComputed[tx * M_Tiles + ty - 1] &&
                     tileComputed[(tx - 1) * M_Tiles + ty - 1])) {
                //condVec[tx * M_Tiles + ty].wait(lock);
            }
        }else if(tx > 0){
            while (!tileComputed[(tx - 1) * M_Tiles + ty]) {
                //condVec[tx * M_Tiles + ty].wait(lock);
            }
        }else if(ty > 0){
            while (!tileComputed[tx * M_Tiles + ty - 1]) {
                //condVec[tx * M_Tiles + ty].wait(lock);
            }
        }

    }

    void releasePermission(){
        tileComputed[tx*M_Tiles + ty] = true;

        /*
        if(tx < N_Tiles - 1) {
            //condVec[(tx+1) * M_Tiles + ty].notify_one();
            indexQueue.push(ind(tx+1,ty));
        }
        if(ty < M_Tiles - 1) {
            //condVec[tx * M_Tiles + ty + 1].notify_one();
            indexQueue.push(ind(tx,ty+1));
        }
        if(tx < N_Tiles - 1 && ty < M_Tiles - 1){
            //condVec[(tx+1)*M_Tiles + ty + 1].notify_one();
            indexQueue.push(ind(tx+1,ty+1));
        }
        */

    }


    uint16 *D;
    int tx, ty;
    uint M, N;
    uint I, J;
    int M_Tiles, N_Tiles;
    const char *x, *y;
    //std::vector<std::mutex>& mutexVec;
    //std::vector<std::condition_variable>& condVec;
    std::atomic_bool* tileComputed;
    barrier *thread_barrier;
    //boost::lockfree::queue<ind>& indexQueue;
    ConcurrentQueue<ind>& indexQueue;
    int tid;
    static int nThreads;
};




#endif //LEVENSHTEIN_WORKER_H
