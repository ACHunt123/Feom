// File: test_input.cpp
#include <string>
#include <iostream>
#include <cstdio>
#include "input.h"

int main(int argc, char *argv[])
{
    std::printf("Running program test_input by Adam Prada.\n");
    std::string filename{};
    if (argc < 2)
    {
        std::printf("Please give a file name as an argument.\n");
        std::exit(1);
    }
    else filename = argv[1];
    Input inp{filename};
    std::cout << "Opening file:\n";
    std::cout << filename << "\n";
    //inp.write_log(String("mylog"));
    inp.output_data(std::cout);

}
