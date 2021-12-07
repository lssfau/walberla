target_sources( pe_coupling
    PRIVATE
    LubricationCorrection.h
    TimeStep.h
    BodiesForceTorqueContainer.cpp
    LubricationCorrection.cpp
    TimeStep.cpp
    BodiesForceTorqueContainer.h
    ForceOnBodiesAdder.h
    ForceTorqueOnBodiesResetter.cpp
    ForceOnBodiesAdder.cpp
    BodySelectorFunctions.h
    ForceTorqueOnBodiesScaler.h
    BodySelectorFunctions.cpp
    ForceTorqueOnBodiesScaler.cpp
    ForceTorqueOnBodiesResetter.h     
    )
