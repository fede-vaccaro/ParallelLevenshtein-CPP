//
// Created by fede on 22/01/19.
//

#include "Worker3.h"

int Worker3::threadCount = 0;
barrier Worker3::sharedBarrier(Worker3::MAX_THREAD_COUNT);