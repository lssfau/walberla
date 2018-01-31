#include <iostream>
#if defined(WALBERLA_USE_STD_FILESYSTEM)
#include <filesystem>
#elif defined(WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM)
#include <experimental/filesystem>
#endif

int main() {
#if defined(WALBERLA_USE_STD_FILESYSTEM)
   std::filesystem::path p("/tmp/test.txt");
#elif defined(WALBERLA_USE_STD_EXPERIMENTAL_FILESYSTEM)
   std::experimental::filesystem::path p("/tmp/test.txt");
#endif
   std::cout << p << std::endl;
   return 0;
}
