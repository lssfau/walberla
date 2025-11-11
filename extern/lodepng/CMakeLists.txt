add_library(
        lodepng
        EXCLUDE_FROM_ALL
        lodepng.cpp
        lodepng.h
)

set_target_properties(lodepng PROPERTIES CXX_CLANG_TIDY "")
target_include_directories(lodepng PUBLIC SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_options( lodepng PRIVATE -w )