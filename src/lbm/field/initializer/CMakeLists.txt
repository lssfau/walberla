target_sources( walberla_lbm
        PUBLIC
        PdfFieldInitializer.h
        PdfFieldInitializer.impl.h
        ExprSystemInitFunction.h
)

target_sources( walberla_lbm
        PRIVATE
        ExprSystemInitFunction.cpp
)
