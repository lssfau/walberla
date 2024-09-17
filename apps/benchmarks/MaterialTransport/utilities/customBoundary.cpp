#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/all.h"

#include "core/Environment.h"
#include "core/grid_generator/SCIterator.h"
#include "core/logging/all.h"
#include "core/math/all.h"
#include "core/waLBerlaBuildInfo.h"

#include "field/AddToStorage.h"
#include "field/vtk/all.h"

namespace walberla
{

using ScalarField_T = GhostLayerField< real_t, 1 >;

template< typename PdeField >
class CustomDirichletBoundaryFunctor
{
 public:
   CustomDirichletBoundaryFunctor(const shared_ptr< StructuredBlockStorage >& blocks, const math::AABB& domainAABB,
                                  const BlockDataID& pdeFieldID,
                                  const CellInterval& interval)
      : blocks_(blocks), domainAABB_(domainAABB), pdeFieldID_(pdeFieldID), interval_(interval)
   {}
   void operator()()
   {
      for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt)
      {
         auto* block = static_cast< blockforest::Block* >(&(*blockIt));

         PdeField* ScalarField = block->template getData< PdeField >(pdeFieldID_);

         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ(
            interval_, real_t boundaryCoord_x = 0.; real_t boundaryCoord_y = 0.; real_t boundaryCoord_z = 0.;

            const auto cellAABB = blocks_->getBlockLocalCellAABB(*block, Cell(x, y, z));
            auto cellCenter     = cellAABB.center();
            const real_t funcVal = 0.05*(1 - std::pow((y/real_c(160)),2));
            //WALBERLA_LOG_INFO_ON_ROOT("funcval "<< x << " " << y <<" " << z);
            ScalarField->get(x, y, z) = 0.05;//real_c(2) * funcVal ;//- ScalarField->get(x + 1, y + 0, z + 0);

         )
      }
   }

 private:
   shared_ptr< StructuredBlockStorage > blocks_;
   math::AABB domainAABB_;
   BlockDataID pdeFieldID_;
   CellInterval interval_;
};

} // namespace walberla
