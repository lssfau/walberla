#include <iostream>
#if defined(WALBERLA_USE_STD_ANY)
#include <any>
#elif defined(WALBERLA_USE_STD_EXPERIMENTAL_ANY)
#include <experimental/any>
#endif

int main() {
#if defined(WALBERLA_USE_STD_ANY)
   auto a = std::any(42);
#elif defined(WALBERLA_USE_STD_EXPERIMENTAL_ANY)
   auto a = std::experimental::any(42);
#endif
   std::cout << std::experimental::any_cast<int>(a) << std::endl;
   return 0;
}
