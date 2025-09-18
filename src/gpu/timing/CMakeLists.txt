target_sources( walberla_gpu
        PUBLIC
        DeviceTimer.h
        DeviceTimingNode.h
        DeviceTimingPool.h
        DeviceTimingTree.h
)

target_sources( walberla_gpu
        PRIVATE
        DeviceTimingNode.cpp
        DeviceTimingPool.cpp
        DeviceTimingTree.cpp
)
