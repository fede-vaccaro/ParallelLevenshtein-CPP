//
// Created by fede on 22/01/19.
//

#ifndef LEVENSHTEIN_WORKER3_H
#define LEVENSHTEIN_WORKER3_H

#include "BoostBarrier.h"
#include <vector>
#include <boost/lockfree/queue.hpp>
#include "ind.h"
#include <string.h>

class Worker3 {
public:
    Worker3(const char *x, const char *y, int *D, barrier *b, boost::lockfree::queue<ind> &indexVector) : D(D), x(x),
                                                                                                          y(y), b(b),
                                                                                                          indexVector(
                                                                                                                  indexVector) {
        tid = threadCount++;
        M = strlen(x);
        N = strlen(y);
    }

    void operator()() {


        while (true) {
            if (!indexVector.empty()) {
                ind index;
                indexVector.pop(index);
                if (index.i == -1) { // poison pill
                    break;
                }
                //if (index.i == 3 && index.j == 3) {
                //    printf("ij=3,3 \n");
                //}
                computeSubMatrix(index.i * TILE_WIDTH + 1, index.j * TILE_WIDTH + 1);
                //printf("Computing submatrix %i%i\n", index.i, index.j);
            } else {
                b->count_down_and_wait();
            }
        }

        //b->count_down_and_wait();
    }


    static const int TILE_WIDTH = 256;
    static const int MAX_THREAD_COUNT = 12;

private:
    void computeSubMatrix(int I, int J) {
        int M_ = M + 1;
        int N_ = N + 1;
        for (int i = I; i < N_ && i < I + TILE_WIDTH; i++) {
            for (int j = J; j < M_ && j < J + TILE_WIDTH; j++) {
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

    int minimum(int i, int j, int k) {
        return std::min(std::min(i, j), k);
    }

    boost::lockfree::queue<ind> &indexVector;
    int *D;
    //int tx, ty;
    int M, N;
    int tid;
    //int M_Tiles, N_Tiles;
    barrier *b;
    const char *x, *y;
    static barrier sharedBarrier;
    static int threadCount;

};


#endif //LEVENSHTEIN_WORKER3_H
