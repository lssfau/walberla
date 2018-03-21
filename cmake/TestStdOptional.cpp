#include <iostream>
#if defined(WALBERLA_USE_STD_OPTIONAL)
#include <optional>
#elif defined(WALBERLA_USE_STD_EXPERIMENTAL_OPTIONAL)
#include <experimental/optional>
#endif

int main() {
#if defined(WALBERLA_USE_STD_OPTIONAL)
   auto a = std::optional<int>();
   auto b = std::optional<int>(42);
#elif defined(WALBERLA_USE_STD_EXPERIMENTAL_OPTIONAL)
   auto a = std::experimental::optional<int>();
   auto b = std::experimental::optional<int>(42);
#endif
   if (b)
      std::cout << a.value_or(b.value()) << std::endl;
   return 0;
}
