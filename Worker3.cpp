//
// Created by fede on 21/01/19.
//

#include "Worker3.h"

int Worker3::nThreads = 0;
barrier Worker3::privateBarrier(MAX_THREAD_COUNT);