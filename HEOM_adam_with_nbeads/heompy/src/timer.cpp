// File: timer.cpp
#include <string>
#include "timer.h"

void Timer::start()
{
    time_start = std::chrono::steady_clock::now();
}

double Timer::end()
{
    time_end = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::seconds> (time_end - time_start).count();
    return time;
}

double Timer::get_seconds()
{
    return time;
}

double Timer::get_minutes()
{
    return time/60.0;
}

double Timer::get_hours()
{
    return time/3600.0;
}


double Timer::get_days()
{
    return time/(24*3600.0);
}

double Timer::get_months()
{
    return time/(30*24*3600.0);
}

double Timer::get_years()
{
    return time/(365.25*24*3600.0);
}

std::string Timer::get_string()
{
    if (time<60)
    {
        return std::to_string(time) + " s";
    }
    else if (time<3600)
    {
        return std::to_string(time/60) + " min";
    }
    else if (time<86400)
    {
        return std::to_string(time/3600) + " h";
    }
    else 
    {
        return std::to_string(time/86400) + " days";
    }
}
