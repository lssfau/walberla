target_sources( walberla_python_coupling
    PRIVATE
    ConfigFromDict.h
    MplHelpers.h
    PybindHelper.h
    SliceToCellInterval.h
    BlockStorageExportHelpers.h
    ModuleInit.cpp
    PythonIterableToStdVector.h
    OwningIterator.h
    ConfigFromDict.cpp
    BlockStorageExportHelpers.cpp
    ModuleInit.h     
    )
