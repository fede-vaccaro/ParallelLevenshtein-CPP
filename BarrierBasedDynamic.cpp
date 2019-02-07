//
// Created by fede on 21/01/19.
//

#include "BarrierBasedDynamic.h"

int BarrierBasedDynamic::nThreads = 0;
barrier BarrierBasedDynamic::privateBarrier(BarrierBasedDynamic::MAX_THREAD_COUNT);