#include "core/all.h"
#include "blockforest/all.h"
#include "field/all.h"

#include "gen/Kernels1D2D.hpp"

using namespace walberla;

void test1D() {
    using Field_T = field::GhostLayerField< real_t, 3 >;

    auto blocks = blockforest::createUniformBlockGrid(
        1, 1, 1,
        32, 1, 1,
        1.0,
        true
    );

    auto fId_fzyx = field::addToStorage< Field_T >(blocks, "f", -1., field::fzyx);
    auto fId_zyxf = field::addToStorage< Field_T >(blocks, "f", -1., field::zyxf);

    Kernels1D::Set_fzyx set_fzyx{ blocks, fId_fzyx };
    Kernels1D::Set_zyxf set_zyxf{ blocks, fId_zyxf };

    for(auto& b: *blocks) {
        auto check = [&](BlockDataID fId)
         {
            Field_T & f = *b.getData< Field_T >(fId);

            CellInterval interior = f.xyzSize();

            for(auto fIt = f.beginWithGhostLayer(); fIt != f.end(); ++fIt) {
                Vector3< real_t > cc{ blocks->getBlockLocalCellCenter(b, fIt.cell()) };
                if(interior.contains(fIt.cell())){
                    WALBERLA_CHECK_FLOAT_EQUAL(*fIt, cc[uint_c(fIt.f())]);
                } else {
                    WALBERLA_CHECK_FLOAT_EQUAL(*fIt, -1.);
                }
            }
        };

        set_fzyx( &b );
        check(fId_fzyx);

        set_zyxf( &b );
        check(fId_zyxf);
    }    
}

void test2D() {
    using Field_T = field::GhostLayerField< real_t, 3 >;

    auto blocks = blockforest::createUniformBlockGrid(
        1, 1, 1,
        32, 32, 1,
        1.0,
        true
    );

    auto fId_fzyx = field::addToStorage< Field_T >(blocks, "f", -1., field::fzyx);
    auto fId_zyxf = field::addToStorage< Field_T >(blocks, "f", -1., field::zyxf);

    Kernels2D::Set_fzyx set_fzyx{ blocks, fId_fzyx };
    Kernels2D::Set_zyxf set_zyxf{ blocks, fId_zyxf };

    for(auto& b: *blocks) {
        auto check = [&](BlockDataID fId)
         {
            Field_T & f = *b.getData< Field_T >(fId);

            CellInterval interior = f.xyzSize();

            for(auto fIt = f.beginWithGhostLayer(); fIt != f.end(); ++fIt) {
                Vector3< real_t > cc{ blocks->getBlockLocalCellCenter(b, fIt.cell()) };
                if(interior.contains(fIt.cell())){
                    WALBERLA_CHECK_FLOAT_EQUAL(*fIt, cc[uint_c(fIt.f())]);
                } else {
                    WALBERLA_CHECK_FLOAT_EQUAL(*fIt, -1.);
                }
            }
        };

        set_fzyx( &b );
        check(fId_fzyx);

        set_zyxf( &b );
        check(fId_zyxf);
    }    
}

int main(int argc, char ** argv) {
    Environment env{ argc, argv };
    test1D();
    test2D();
}

