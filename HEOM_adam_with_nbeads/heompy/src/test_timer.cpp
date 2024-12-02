// File: test_timer.cpp
#include "timer.h"
#include <unistd.h> // for usleep
#include <iostream>
#include <thread>

int main()
{
    double time;
    Timer timer;
    std::cout << "The test_timer program" << std::endl;
    std::cout << "Please enter the time you want to wait in seconds" << std::endl;
    timer.start();
    std::cin >> time;
    timer.end();
    std::cout << "Time it took you to type the number : " 
    << timer.get_seconds()
    << " sec\n";

    time = 1.0e6*time;

    std::cout << "Waiting..." << std::endl;
    timer.start();
    usleep(time);
    timer.end();

    std::cout << "Elapsed time in seconds : " 
    << timer.get_seconds()
    << " sec\n";

    std::cout << "Elapsed time in minutes : " 
    << timer.get_minutes()
    << " min\n";

    std::cout << "Elapsed time in hours : " 
    << timer.get_hours()
    << " h\n";

    std::cout << "Elapsed time in days : " 
    << timer.get_days()
    << " days\n";

    std::cout << "Elapsed time in months : " 
    << timer.get_months()
    << " months\n";

    std::cout << "Elapsed time in years : " 
    << timer.get_years()
    << " years\n";

    return 0;
}
