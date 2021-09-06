#include <iostream>

int main()
{
   int min = -10;

#pragma omp parallel for
   for (int i = min; i <= 10; ++i)
   {
      std::cout << i << std::endl;
   }

   return 0;
}