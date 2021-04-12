#include "core/math/Random.h"
#include "domain_decomposition/StructuredBlockStorage.h"


namespace walberla {


inline void initShearVelocity(const shared_ptr<StructuredBlockStorage> & blocks, BlockDataID velFieldID,
                              const real_t xMagnitude=real_t(0.005), const real_t fluctuationMagnitude=real_t(0.05) )
{
    math::seedRandomGenerator(0);
    auto halfZ = blocks->getDomainCellBB().zMax() / 2;
    for( auto & block: *blocks)
    {
        auto velField = block.getData<GhostLayerField<real_t, 3> >( velFieldID );
        WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(velField,
                                                         Cell globalCell;
        blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
        real_t randomReal = xMagnitude * math::realRandom<real_t>(-fluctuationMagnitude, fluctuationMagnitude);
        velField->get(x, y, z, 1) = real_t(0);
        velField->get(x, y, z, 2) = randomReal;

        if( globalCell[2] >= halfZ ) {
            velField->get(x, y, z, 0) = xMagnitude;
        } else {
            velField->get(x, y, z, 0) = -xMagnitude;
        }
        );
    }
}


}