target_sources( walberla_core
        PUBLIC
        CheckFunctions.h
        CheckFunctions.impl.h
        Debug.h
        PrintStacktrace.h
        TestSubsystem.h
        all.h
        demangle.h
)

target_sources( walberla_core
        PRIVATE
        CMakeLists.txt
        CheckFunctions.cpp
        Debug.cpp
        PrintStacktrace.cpp
        TestSubsystem.cpp
)