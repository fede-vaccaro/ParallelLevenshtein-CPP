//
// Created by fede on 21/01/19.
//

#include "BarrierBasedStatic.h"

int BarrierBasedStatic::nThreads = 0;
barrier BarrierBasedStatic::privateBarrier(MAX_THREAD_COUNT);