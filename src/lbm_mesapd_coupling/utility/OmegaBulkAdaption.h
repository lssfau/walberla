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
//! \file OmegaBulkAdaption.h
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "mesa_pd/common/AABBConversion.h"
#include "domain_decomposition/StructuredBlockStorage.h"

namespace walberla {
namespace lbm_mesapd_coupling {

// utility functions

real_t bulkViscosityFromOmegaBulk(real_t omegaBulk)
{
   return real_t(2) / real_t(9) * ( real_t(1) / omegaBulk - real_t(0.5) );
}


real_t omegaBulkFromBulkViscosity(real_t bulkViscosity)
{
   return real_t(2) / ( real_t(9) * bulkViscosity + real_t(1) );
}

// see Khirevich et al. - Coarse- and fine-grid numerical behavior of MRT/TRT lattice-Boltzmann schemes in regular and random sphere packings
// LambdaBulk is the "magic parameter" here, i.e. the ratio between Lambda_e and Lambda_nu, Eq. 19
real_t omegaBulkFromOmega(real_t omega, real_t LambdaBulk = real_t(1))
{
   return real_t(1) / (LambdaBulk * ( real_t(1) / omega - real_t(1)/ real_t(2) ) + real_t(1)/ real_t(2) );
}

/*
 * Adapts all cells that are inside a sphere of radius interactionRadius + adaptionLayerSize around the particle position (center)
 */
template< typename ParticleAccessor_T, typename ParticleSelector_T >
class OmegaBulkAdapter
{
   static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

public:
   OmegaBulkAdapter(const shared_ptr<StructuredBlockStorage> & blockStorage,
                    const BlockDataID & omegaBulkFieldID,
                    const shared_ptr<ParticleAccessor_T>& ac,
                    const real_t defaultOmegaBulk,
                    const real_t adaptedOmegaBulk,
                    const real_t adaptionLayerSize,
                    ParticleSelector_T particleSelector)
         : blockStorage_( blockStorage ),  omegaBulkFieldID_( omegaBulkFieldID ), ac_( ac ),
           defaultOmegaBulk_(defaultOmegaBulk),
           adaptedOmegaBulk_(adaptedOmegaBulk),
           adaptionLayerSize_(adaptionLayerSize),
           particleSelector_( particleSelector )
   {}


   void inline operator()( IBlock * const block )
   {

      using namespace lbm_mesapd_coupling;

      WALBERLA_ASSERT_NOT_NULLPTR( block );

      auto * omegaBulkField = block->getData< GhostLayerField< real_t, 1> >( omegaBulkFieldID_ );

      WALBERLA_ASSERT_NOT_NULLPTR( omegaBulkField );

      // reset whole field to default
      WALBERLA_FOR_ALL_CELLS_XYZ(omegaBulkField, omegaBulkField->get(x,y,z) = defaultOmegaBulk_; )

      const auto blockLevel = blockStorage_->getLevel(*block);
      const real_t dx = blockStorage_->dx( blockLevel );
      const real_t dy = blockStorage_->dy( blockLevel );
      const real_t dz = blockStorage_->dz( blockLevel );

      for( size_t particleIdx = 0; particleIdx < ac_->size(); ++particleIdx )
      {
         // currently only works for finite particles
         if (particleSelector_(particleIdx, *ac_) && !mesa_pd::data::particle_flags::isSet( ac_->getFlags(particleIdx), mesa_pd::data::particle_flags::INFINITE))
         {

            auto particlePosition = ac_->getPosition(particleIdx);
            auto particleInteractionRadius = ac_->getInteractionRadius(particleIdx);
            auto extendedParticleAABB = mesa_pd::getAABBFromInteractionRadius(particlePosition, particleInteractionRadius).getExtended(adaptionLayerSize_);
            auto cellBB = getCellBBFromAABB( extendedParticleAABB, false, *block, *blockStorage_, uint_t(0));
            // no ghost layers are included -> omega field has to be communicated afterwards if ghost layer values are needed in LBM sweep

            auto adaptionRadius = particleInteractionRadius + adaptionLayerSize_;

            Vector3<real_t> startCellCenter = blockStorage_->getBlockLocalCellCenter( *block, cellBB.min() );

            real_t cz = startCellCenter[2];
            for( cell_idx_t z = cellBB.zMin(); z <= cellBB.zMax(); ++z )
            {
               real_t cy = startCellCenter[1];
               for( cell_idx_t y = cellBB.yMin(); y <= cellBB.yMax(); ++y )
               {
                  real_t cx = startCellCenter[0];
                  for( cell_idx_t x = cellBB.xMin(); x <= cellBB.xMax(); ++x )
                  {
                     auto cellCenter = Vector3<real_t>(cx,cy,cz);
                     auto sqrDistanceToParticleCenter = (cellCenter-particlePosition).sqrLength();

                     if(sqrDistanceToParticleCenter < adaptionRadius * adaptionRadius)
                     {
                        // change omega to given adaption value
                        omegaBulkField->get(x,y,z) = adaptedOmegaBulk_;
                     }
                     cx += dx;
                  }
                  cy += dy;
               }
               cz += dz;
            }
         }
      }
   }

   void setDefaultOmegaBulk(real_t defaultOmegaBulk)
   {
      defaultOmegaBulk_ = defaultOmegaBulk;
   }
   real_t getDefaultOmegaBulk()
   {
      return defaultOmegaBulk_;
   }

   void setAdaptedOmegaBulk(real_t adaptedOmegaBulk)
   {
      adaptedOmegaBulk_ = adaptedOmegaBulk;
   }
   real_t getAdaptedOmegaBulk()
   {
      return adaptedOmegaBulk_;
   }

   void setAdaptionLayerSize( real_t adaptionLayerSize )
   {
      adaptionLayerSize_ = adaptionLayerSize;
   }
   real_t getAdaptionLayerSize()
   {
      return adaptionLayerSize_;
   }

private:

   shared_ptr<StructuredBlockStorage> blockStorage_;
   BlockDataID omegaBulkFieldID_;
   shared_ptr<ParticleAccessor_T> ac_;
   real_t defaultOmegaBulk_;
   real_t adaptedOmegaBulk_;
   real_t adaptionLayerSize_;
   ParticleSelector_T particleSelector_;

};

} // namespace lbm_mesapd_coupling
} // namespace walberla
