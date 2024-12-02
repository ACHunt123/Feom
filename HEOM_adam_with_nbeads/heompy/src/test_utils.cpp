// File: test_utils.cpp
#include "general.h"
#include <iostream>
#include "utils.h"
int main()
{
    String testfile = "testfile";
    Real num1 = 11.11;
    Real num2 = 22.22;
    String message = "message\n";
    Vec vec = {101.1, 202.2, 303.3};

    std::cout << "writing num1\n message\n num1 num2\n numlist"
        << ", num1, num2, vec\n num1 vec\n vec\n" ; 
    write_out(testfile, num1);
    write_out(testfile, message);
    write_out(testfile, num1, num2);
    write_out(testfile, {num1,num2});
    write_out(testfile, num1, num2, vec);
    write_out(testfile, num1, vec);
    write_out(testfile, vec);

}
