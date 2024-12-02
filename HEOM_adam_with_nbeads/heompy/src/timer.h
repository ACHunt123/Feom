// File: timer.h
#ifndef TIMER_H
#define TIMER_H
#include <chrono>
#include <string>

class Timer
{
public:
    std::chrono::time_point<std::chrono::steady_clock> time_start, time_end;
    double time{0.0};
    // Constructors
    Timer() {start();}
    // Member functions
    void start();
    double end();
    double get_seconds();
    double get_minutes();
    double get_hours();
    double get_days();
    double get_months();
    double get_years();
    std::string get_string();
};
#endif
