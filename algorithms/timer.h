//
// Created by cbl on 4/23/22.
//

#ifndef SIM_TIMER_H
#define SIM_TIMER_H

#include <cstdio>
#include <ctime>
#include <sys/time.h>
#include <boost/timer/timer.hpp>

class AutoTimer {
private:
    clock_t begin;
    clock_t end;
public:
    AutoTimer() {
        begin = clock();
        end = 0;
    }
    ~AutoTimer() {
        end = clock();
        float t = static_cast<float>(end - begin) * 1000 / CLOCKS_PER_SEC;
        std::cout << t << std::endl;
    }
};

class ClockTimer {
private:
    clock_t begin;
    clock_t end;
    clock_t checkpoint_start;
    clock_t checkpoint_end;
    clock_t cumulate_time;
public:
    ClockTimer() {
        begin = 0;
        end = 0;
        checkpoint_start = 0;
        checkpoint_end = 0;
        cumulate_time = 0;
    }
    void start() {begin = clock();}
    void stop() {
        end = clock();
        if (begin == 0 || end == 0) {
            printf("\ntimer: begin or stop method not called!\n");
            return;
        }
        // printf("\nElapsed time %.8f ms\n", (float) (end - begin) * 1000.0 / CLOCKS_PER_SEC);
        printf("%.17f\n", (float) (end - begin) * 1000.0 / CLOCKS_PER_SEC);
    }
    void cumulate_start_point() {
        checkpoint_start = clock();
    }
    void cumulate_end_point() {
        checkpoint_end = clock();
        cumulate_time += (checkpoint_end - checkpoint_start);
    }
    void stop_cumulation() const {
        printf("%.17f\n", (float) cumulate_time * 1000.0 / CLOCKS_PER_SEC);
    }
};

class WallTimer {
private:
    long cur_time_;
public:
    WallTimer() {
        cur_time_ = 0;
    }
    void start() {
        struct timeval time{};
        if(gettimeofday( &time, 0 )) return;
        cur_time_ = 1000000 * time.tv_sec + time.tv_usec;
    }
    void stop() {
        struct timeval time{};
        if(gettimeofday( &time, 0 )) return;

        long cur_time = 1000000 * time.tv_sec + time.tv_usec;
        double sec = (double) (cur_time - cur_time_) / 1000000.0;
        if(sec < 0) sec += 86400;
        cur_time_ = cur_time;
        printf("Elapsed time %.8f ms\n", 1000.0 * sec);
    }
};
#endif //SIM_TIMER_H
