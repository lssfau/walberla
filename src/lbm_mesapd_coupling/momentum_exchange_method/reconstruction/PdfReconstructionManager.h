//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file PdfReconstructionManager.h
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================


#pragma once

#include "core/debug/Debug.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/Field.h"
#include "field/FlagField.h"
#include "field/iterators/IteratorMacros.h"

#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/Reconstructor.h"
#include "lbm_mesapd_coupling/utility/ParticleFunctions.h"

#include "mesa_pd/data/IAccessor.h"

#include <functional>

namespace walberla {
namespace lbm_mesapd_coupling {


//**************************************************************************************************************************************
/*!\brief Class to manage the reconstruction of PDFs that is needed when cells are becoming uncovered by moving obstacles.
 *
 * Due to the explicit mapping of particles into the fluid domain via flags, the PDFs of cells that turned from obstacle to fluid
 * are missing and must be reconstructed in order to continue with the simulation.
 * This class is to be used as a sweep in a LBM simulation with moving obstacles and calls for each cell that is tagged as
 * 'formerObstacle' the specified reconstructor (see Reconstructor.h for the available variants).
 * After the successful reconstruction of all PDFs, the flags are updated to 'fluid'.
 * For small obstacle fractions, an optimized variant is available that only looks for 'formerObstacle' cells in the vicinity
 * of available particles. It is activated via the 'optimizeForSmallObstacleFraction' argument in the constructor.
 *
 *
 * Note:
 * A potential problem arises when the cell that is to be reconstructed, i.e. filled with fluid, is surrounded (given by the LBM stencil)
 * by only boundary cells. Consequently, the successive calls to LBM will simulate a local system containing a single fluid cell
 * surrounded by boundary cells with different velocity boundary conditions. This system is not divergence free (as required by the mass
 * conservation equation) and will thus lead to mass loss/gain depending on the boundary conditions.
 * This is a very rare case but if it happens, it will almost always lead to malicious oscillations in the density.
 * Note that this can also happen if several cells are locally "trapped" by surrounding boundary cells. Such a case can not be detected
 * easily and also not be curred easily.
 * It has to be inspected whether the behavior arising from this issue leads to problems in real simulations.
 */
//**************************************************************************************************************************************

template< typename PdfField_T, typename BoundaryHandling_T, typename ParticleAccessor_T, typename Reconstructor_T>
class PdfReconstructionManager
{
   static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

public:

   using FlagField_T = typename BoundaryHandling_T::FlagField;
   using flag_t = typename BoundaryHandling_T::flag_t;

   inline PdfReconstructionManager( const shared_ptr<StructuredBlockStorage> & blockStorage,
                                    const BlockDataID & pdfFieldID,
                                    const BlockDataID & boundaryHandlingID,
                                    const BlockDataID & particleFieldID,
                                    const shared_ptr<ParticleAccessor_T> & ac,
                                    const FlagUID & formerObstacle, const FlagUID & fluid,
                                    const shared_ptr<Reconstructor_T> & reconstructor,
                                    bool conserveMomentumWhenRestoring,
                                    const bool optimizeForSmallObstacleFraction = false ) :
      blockStorage_( blockStorage ), pdfFieldID_( pdfFieldID), boundaryHandlingID_( boundaryHandlingID ),
      particleFieldID_( particleFieldID ), ac_( ac ),
      reconstructor_ ( reconstructor ), formerObstacle_( formerObstacle ), fluid_( fluid ),
      conserveMomentumWhenRestoring_( conserveMomentumWhenRestoring ),
      optimizeForSmallObstacleFraction_( optimizeForSmallObstacleFraction )
   {}

   void inline operator()( IBlock * const block );

private:



   shared_ptr<StructuredBlockStorage> blockStorage_;

   const BlockDataID pdfFieldID_;
   const BlockDataID boundaryHandlingID_;
   const BlockDataID particleFieldID_;

   shared_ptr<ParticleAccessor_T> ac_;

   shared_ptr<Reconstructor_T> reconstructor_;

   const FlagUID formerObstacle_;
   const FlagUID fluid_;

   bool conserveMomentumWhenRestoring_;
   const bool optimizeForSmallObstacleFraction_;

};


template< typename PdfField_T, typename BoundaryHandling_T, typename ParticleAccessor_T, typename Reconstructor_T>
void PdfReconstructionManager< PdfField_T, BoundaryHandling_T, ParticleAccessor_T, Reconstructor_T >
::operator()( IBlock * const block )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   auto * pdfField         = block->getData< PdfField_T >( pdfFieldID_ );
   auto * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );
   auto * flagField        = boundaryHandling->getFlagField();
   auto * particleField    = block->getData< ParticleField_T >( particleFieldID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );
   WALBERLA_ASSERT_NOT_NULLPTR( flagField );
   WALBERLA_ASSERT_NOT_NULLPTR( particleField );

   WALBERLA_ASSERT_EQUAL( particleField->xyzSize(), flagField->xyzSize() );

   WALBERLA_ASSERT( flagField->flagExists( formerObstacle_ ) );
   WALBERLA_ASSERT( flagField->flagExists( fluid_ ) );

   const flag_t formerObstacle = flagField->getFlag( formerObstacle_ );
   const flag_t fluid          = flagField->getFlag( fluid_ );

   // reconstruct all missing PDFs (only inside the domain, ghost layer values get communicated)

   if( optimizeForSmallObstacleFraction_ )
   {
      const uint_t numberOfGhostLayersToInclude = uint_t(0);

      for( size_t idx = 0; idx < ac_->size(); ++idx )
      {
         auto cellIntervalForReconstruction = getParticleCellBB( idx, *ac_, *block, *blockStorage_, numberOfGhostLayersToInclude );

         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ(cellIntervalForReconstruction,
                                                size_t particleIdx = ac_->uidToIdx(particleField->get(x,y,z)); // to only reconstruct cells that belonged to this particle
                                                if (isFlagSet(flagField->get(x,y,z), formerObstacle) && particleIdx == idx) {
                                                   (*reconstructor_)(block, x, y, z, pdfField, idx, *ac_);
                                                   if( conserveMomentumWhenRestoring_ && pdfField->isInInnerPart(Cell(x,y,z)) ) {
                                                      // the (artificially) added momentum in the restored fluid cell has to be subtracted from the former particle to ensure momentum conservation
                                                      // force = momentum / dt, with dt = 1
                                                      Vector3<real_t> force = pdfField->getMomentumDensity(x, y, z);
                                                      Vector3< real_t > cellCenter = blockStorage_->getBlockLocalCellCenter(*block, Cell(x,y,z) );
                                                      lbm_mesapd_coupling::addHydrodynamicForceAtWFPosAtomic( particleIdx, *ac_, -force, cellCenter );
                                                   }
                                                }

         )
      }
   }
   else
   {
      WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ(flagField->xyzSize(),
                                             if (isFlagSet(flagField->get(x,y,z), formerObstacle)) {
                                                size_t particleIdx = ac_->uidToIdx(particleField->get(x,y,z));
                                                WALBERLA_ASSERT_UNEQUAL( particleIdx, ac_->getInvalidIdx(), "Index of particle is invalid!" );
                                                (*reconstructor_)(block, x, y, z, pdfField, particleIdx, *ac_);
                                                if( conserveMomentumWhenRestoring_ && pdfField->isInInnerPart(Cell(x,y,z)) ) {
                                                   // the (artificially) added momentum in the restored fluid cell has to be subtracted from the former particle to ensure momentum conservation
                                                   // force = momentum / dt, with dt = 1
                                                   Vector3<real_t> force = pdfField->getMomentumDensity(x, y, z);
                                                   Vector3< real_t > cellCenter = blockStorage_->getBlockLocalCellCenter(*block, Cell(x,y,z) );
                                                   lbm_mesapd_coupling::addHydrodynamicForceAtWFPosAtomic( particleIdx, *ac_, -force, cellCenter );
                                                }
                                             };
      )
   }

   // update the flags from formerObstacle to fluid (inside domain & in ghost layers), invalidate entry in particle field
   if( optimizeForSmallObstacleFraction_ )
   {
      const uint_t numberOfGhostLayersToInclude = flagField->nrOfGhostLayers();

      for( size_t idx = 0; idx < ac_->size(); ++idx )
      {
         auto cellIntervalToCleanUp = getParticleCellBB( idx, *ac_, *block, *blockStorage_, numberOfGhostLayersToInclude );
         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ(cellIntervalToCleanUp,
                                                if (isFlagSet(flagField->get(x,y,z), formerObstacle)) {
                                                   boundaryHandling->setDomain( fluid, x, y, z );
                                                   removeFlag( flagField->get(x,y,z), formerObstacle );
                                                   particleField->get(x,y,z) = ac_->getInvalidUid();
                                                }
         )
      }
   }
   else
   {
      WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ(flagField->xyzSizeWithGhostLayer(),
                                             if (isFlagSet(flagField->get(x,y,z), formerObstacle)) {
                                                boundaryHandling->setDomain( fluid, x, y, z );
                                                removeFlag( flagField->get(x,y,z), formerObstacle );
                                                particleField->get(x,y,z) = ac_->getInvalidUid();
                                             }
      )
   }
}


template< typename PdfField_T, typename BoundaryHandling_T, typename ParticleAccessor_T, typename Reconstructor_T >
shared_ptr< PdfReconstructionManager<PdfField_T,BoundaryHandling_T,ParticleAccessor_T,Reconstructor_T> >
makePdfReconstructionManager( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & pdfFieldID,
                              const BlockDataID & boundaryHandlingID, const BlockDataID & particleFieldID, const shared_ptr<ParticleAccessor_T> & ac,
                              const FlagUID & formerObstacle, const FlagUID & fluid,
                              const shared_ptr<Reconstructor_T> & reconstructor,
                              bool conserveMomentumWhenRestoring,
                              const bool optimizeForSmallObstacleFraction = false)
{
   using RM_T = PdfReconstructionManager<PdfField_T,BoundaryHandling_T,ParticleAccessor_T,Reconstructor_T>;
   return make_shared<RM_T>( RM_T(blockStorage, pdfFieldID, boundaryHandlingID, particleFieldID, ac,formerObstacle, fluid, reconstructor, conserveMomentumWhenRestoring, optimizeForSmallObstacleFraction) );
}


template< typename PdfField_T, typename BoundaryHandling_T, typename ParticleAccessor_T >
shared_ptr< PdfReconstructionManager<PdfField_T,BoundaryHandling_T,ParticleAccessor_T,EquilibriumReconstructor<BoundaryHandling_T> > >
makePdfReconstructionManager( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & pdfFieldID,
                              const BlockDataID & boundaryHandlingID, const BlockDataID & particleFieldID, const shared_ptr<ParticleAccessor_T> & ac,
                              const FlagUID & formerObstacle, const FlagUID & fluid,
                              bool conserveMomentumWhenRestoring,
                              const bool optimizeForSmallObstacleFraction = false)
{
   using RM_T = PdfReconstructionManager<PdfField_T,BoundaryHandling_T,ParticleAccessor_T,EquilibriumReconstructor<BoundaryHandling_T> >;
   return make_shared<RM_T>( RM_T(blockStorage, pdfFieldID, boundaryHandlingID, particleFieldID, ac,formerObstacle, fluid,
                                  makeEquilibriumReconstructor<BoundaryHandling_T>(blockStorage,boundaryHandlingID), conserveMomentumWhenRestoring, optimizeForSmallObstacleFraction) );
}


} // namespace lbm_mesapd_coupling
} // namespace walberla
