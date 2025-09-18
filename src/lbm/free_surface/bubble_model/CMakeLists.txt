target_sources( walberla_lbm
        PUBLIC
        Bubble.h
        BubbleModel.impl.h
        BubbleDefinitions.h
        BubbleDistanceAdaptor.h
        BubbleIDFieldPackInfo.h
        BubbleModel.h
        BubbleModelFromConfig.h
        BubbleModelFromConfig.impl.h
        DisjoiningPressureBubbleModel.h
        DisjoiningPressureBubbleModel.impl.h
        DistanceInfo.h
        FloodFill.h
        FloodFill.impl.h
        Geometry.h
        Geometry.impl.h
        MergeInformation.h
        NewBubbleCommunication.h
        RegionalFloodFill.h
)

target_sources( walberla_lbm
        PRIVATE
        DistanceInfo.cpp
        MergeInformation.cpp
        NewBubbleCommunication.cpp
)
