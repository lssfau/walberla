#pragma once

#include "gpu/ErrorChecking.h"

namespace test_allocators {

template<typename T>
void checkVector(T* ptr, size_t size, bool* flag);

} // namespace test_allocators
