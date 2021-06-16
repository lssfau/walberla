add_library(
      lodepng
      EXCLUDE_FROM_ALL
      lodepng.cpp
      lodepng.h
)

target_include_directories(lodepng PUBLIC SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}")
if (MSVC)
   target_compile_options( lodepng PRIVATE /w )
else ()
   target_compile_options( lodepng PRIVATE -w )
endif ()