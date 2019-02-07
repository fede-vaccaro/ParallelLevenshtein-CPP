//
// Created by fede on 21/01/19.
//

#include "Worker2.h"

int Worker2::nThreads = 0;
barrier Worker2::privateBarrier(Worker2::MAX_THREAD_COUNT);