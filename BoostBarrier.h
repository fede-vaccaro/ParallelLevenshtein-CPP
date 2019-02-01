//
// Created by fede on 19/01/19.
// From libboost
//

#ifndef LEVENSHTEIN_BARRIER_H
#define LEVENSHTEIN_BARRIER_H

#include <thread>
#include <mutex>
#include <condition_variable>
#
class barrier
{
public:
    barrier(const barrier&) = delete;
    barrier& operator=(const barrier&) = delete;
    explicit barrier(unsigned int count) :
            m_count(count), m_generation(0),
            m_count_reset_value(count)
    {
    }
    void count_down_and_wait()
    {
        std::unique_lock< std::mutex > lock(m_mutex);
        unsigned int gen = m_generation;
        if (--m_count == 0)
        {
            m_generation++;
            m_count = m_count_reset_value;
            m_cond.notify_all();
            return;
        }
        while (gen == m_generation)
            m_cond.wait(lock);
    }
private:
    std::mutex m_mutex;
    std::condition_variable m_cond;
    unsigned int m_count;
    unsigned int m_generation;
    unsigned int m_count_reset_value;
};

#endif //LEVENSHTEIN_BARRIER_H
