target_sources( walberla_mesa_pd
        PUBLIC
        TensorGlyph.h
        OutputSelector.h
        WriteOutput.h
        ParticleVtkOutput.h
)

target_sources( walberla_mesa_pd
        PRIVATE
        TensorGlyph.cpp
        ParticleVtkOutput.cpp
)
