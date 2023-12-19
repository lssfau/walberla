#include <iostream>


int main()
{
   static_assert(std::is_floating_point_v<_Float16>);
}