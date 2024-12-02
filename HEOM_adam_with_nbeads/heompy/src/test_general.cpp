// File: test_general.cpp
#include <iostream>
#include "general.h"

void vec_fun(Vec vec)
{
    std::cout << "Original vector\n";
    std::cout << vec <<"\n";
    vec(1) = 99;
    std::cout << "new vector\n";
    std::cout << vec <<"\n";
}

int main()
{
    Int compute_until{};
    std::cout << "Enter the highest factorial you want to compute.\n";
    std::cin >> compute_until;
    for (Int i=0;i<=compute_until;i++)
    {
        std::cout << i << "           " << factorial(i) <<"\n";
    }
    for (Int i=0;i<=compute_until;i++)
    {
        for (Int j=0; j<=i; j++)
        {
        std::cout <<  "  " << binom(i,j) ;
        }
        std::cout <<"\n";
    }
}
