//
// Created by fede on 21/01/19.
//

#include "Worker2.h"

int Worker2::threadCount = 0;
barrier Worker2::sharedBarrier(Worker2::MAX_THREAD_COUNT);