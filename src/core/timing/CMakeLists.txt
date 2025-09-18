target_sources( walberla_core
        PUBLIC
        CpuPolicy.h
        ReduceType.h
        RemainingTimeLogger.h
        StaticPolicy.h
        Time.h
        Timer.h
        TimingJSON.h
        TimingNode.h
        TimingPool.h
        TimingTree.h
        WcPolicy.h
)

target_sources( walberla_core
        PRIVATE
        StaticPolicy.cpp
        Time.cpp
        TimingNode.cpp
        TimingPool.cpp
        TimingTree.cpp
)
