add_library(
      sqlite3
      EXCLUDE_FROM_ALL
      sqlite3.c
      sqlite3.h
)

add_definitions( -DSQLITE_OMIT_LOAD_EXTENSION )

set_target_properties(sqlite3 PROPERTIES CXX_CLANG_TIDY "")
target_include_directories(sqlite3 PUBLIC SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_options( lodepng PRIVATE -w )