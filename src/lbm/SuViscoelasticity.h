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
//! \file SuViscoelasticity.h
//! \ingroup lbm
//! \author Cameron Stewart <cstewart@icp.uni-stuttgart.de>
//! \brief Oldroyd-B viscoelasticity extended from Su et al. Phys. Rev. E, 2013
//
//======================================================================================================================


#include "blockforest/communication/UniformBufferedScheme.h"
#include "core/math/Matrix3.h"
#include <field/AddToStorage.h>
#include "field/GhostLayerField.h"
#include "field/communication/PackInfo.h"
#include "lbm/field/PdfField.h"

namespace walberla {
namespace lbm {
namespace viscoelastic {

template < typename LatticeModel_T , typename BoundaryHandling_T>
class Su {

public:
            
   typedef GhostLayerField< Matrix3<real_t>, 1> StressField_T;
   typedef GhostLayerField< Vector3<real_t>, 1> VelocityField_T;
   typedef GhostLayerField< Vector3<real_t>, 1> ForceField_T;
   typedef PdfField< LatticeModel_T > PdfField_T;
   typedef shared_ptr< StructuredBlockForest > Blocks_T;


   Su( Blocks_T blocks, BlockDataID force, BlockDataID pdfId, BlockDataID boundaryHandlingId, BlockDataID stressId, BlockDataID stressOldId, real_t lambda_p, real_t eta_p,
       uint_t period = uint_c(1), bool compressibleFlag = false): blocks_( blocks ), forceId_( force ), pdfId_(pdfId), boundaryHandlingId_(boundaryHandlingId), stressId_(stressId), stressOldId_(stressOldId),
                                   inv_lambda_p_( real_c(1.0)/lambda_p ), eta_p_( eta_p ), delta_t_( real_c(period) ), executionCount_( 0 ), communicateStress_(blocks), 
                                   communicateVelocities_(blocks), compressibleFlag_(compressibleFlag)
   {
      // create velocity fields
      velocityId_ = walberla::field::addToStorage<VelocityField_T>( blocks_, "Velocity Field", Vector3<real_t>(0.0), field::zyxf, uint_c(1));

      // stress and velocity communication scheme
      communicateStress_.addPackInfo( make_shared< field::communication::PackInfo<StressField_T> >( stressId_ ));
      communicateVelocities_.addPackInfo( make_shared< field::communication::PackInfo<VelocityField_T> >( velocityId_ ));

      WALBERLA_ASSERT_GREATER_EQUAL(delta_t_, real_c(1.0));
   }

   Su( Blocks_T blocks, BlockDataID force, BlockDataID pdfId, BlockDataID boundaryHandlingId, BlockDataID stressId, BlockDataID stressOldId, BlockDataID velocityId, real_t lambda_p, real_t eta_p,
       uint_t period = uint_c(1), bool compressibleFlag = false): blocks_( blocks ), forceId_( force ), pdfId_(pdfId), boundaryHandlingId_(boundaryHandlingId), stressId_(stressId), stressOldId_(stressOldId),
                                   velocityId_(velocityId), inv_lambda_p_( real_c(1.0)/lambda_p ), eta_p_( eta_p ), delta_t_( real_c(period) ), executionCount_( 0 ), communicateStress_(blocks), 
                                   communicateVelocities_(blocks), compressibleFlag_(compressibleFlag)
   {
      // stress and velocity communication scheme
      communicateStress_.addPackInfo( make_shared< field::communication::PackInfo<StressField_T> >( stressId_ ));
      communicateVelocities_.addPackInfo( make_shared< field::communication::PackInfo<VelocityField_T> >( velocityId_ ));

      WALBERLA_ASSERT_GREATER_EQUAL(delta_t_, real_c(1.0));
   }

   void calculateStresses(IBlock * block) {
      StressField_T *stressNew = block->getData<StressField_T>(stressId_);
      StressField_T *stressOld = block->getData<StressField_T>(stressOldId_);
      VelocityField_T *velocity = block->getData<VelocityField_T>(velocityId_);
      PdfField_T *pdf = block->getData<PdfField_T>(pdfId_);
      BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingId_ );

      WALBERLA_ASSERT_GREATER_EQUAL(stressOld->nrOfGhostLayers(), 2);
      WALBERLA_ASSERT_GREATER_EQUAL(velocity->nrOfGhostLayers(), 1);

      WALBERLA_FOR_ALL_CELLS_XYZ(stressNew, {
         Cell cell(x,y,z);
         Matrix3<real_t> stress1 = Matrix3<real_t>(0.0);
         Matrix3<real_t> stress2 = Matrix3<real_t>(0.0);
         Matrix3<real_t> stress3 = Matrix3<real_t>(0.0);
         Matrix3<real_t> relstr  = Matrix3<real_t>(0.0);
         bool nearBoundaryFlag = false;

         // if cell is a fluid cell then calculate the stress tensor
         if (boundaryHandling->isDomain(cell)) {
            // check if near a wall
            if (boundaryHandling->isNearBoundary(cell)) {
               nearBoundaryFlag = true;
            }
            else {
               for (auto d = LatticeModel_T::Stencil::beginNoCenter(); d != LatticeModel_T::Stencil::end(); ++d) {
                  Cell cell1 = cell - *d;
                  if (boundaryHandling->isNearBoundary(cell1)) {
                     nearBoundaryFlag = true;
                  }
               }
            }
            for (auto d = LatticeModel_T::Stencil::beginNoCenter(); d != LatticeModel_T::Stencil::end(); ++d) {
               Cell cell1 = cell - *d;
               Cell cell2 = cell1 - *d;
               Cell cell3 = cell + *d;

               // check if using compressible of incompressible algorithm
               if (compressibleFlag_) {
                  if (nearBoundaryFlag) {
                     if (boundaryHandling->isDomain(cell1)) {
                        stress1 += (stressOld->get(cell1) * pdf->get(cell, d.toIdx())) * real_c(1 / 2.0);
                     }
                     if (boundaryHandling->isDomain(cell3)) {
                        stress2 += (stressOld->get(cell) * pdf->get(cell, d.toIdx())) * real_c(1 / 1.5);
                     }
                  }
                  else {
                     stress1 += stressOld->get(cell1) * pdf->get(cell, d.toIdx());
                     stress2 += stressOld->get(cell) * pdf->get(cell, d.toIdx());
                     stress3 += stressOld->get(cell2) * pdf->get(cell, d.toIdx());
                  }
               }
               else {
                  if (nearBoundaryFlag) {
                     if (boundaryHandling->isDomain(cell1)) {
                        stress1 += (stressOld->get(cell1) * pdf->get(cell1, d.toIdx())) * real_c(1 / 2.0);
                     }
                     if (boundaryHandling->isDomain(cell3)) {
                        stress2 += (stressOld->get(cell) * pdf->get(cell, d.toIdx())) * real_c(1 / 1.5);
                     }
                  }
                  else {
                     stress1 += stressOld->get(cell1) * pdf->get(cell1, d.toIdx());
                     stress2 += stressOld->get(cell) * pdf->get(cell, d.toIdx());
                     stress3 += stressOld->get(cell2) * pdf->get(cell2, d.toIdx());
                  }
               }
            }
            Matrix3<real_t> gradu = Matrix3<real_t>(0.0);

            // compute velocity gradient
            for (auto d = LatticeModel_T::Stencil::beginNoCenter(); d.direction() != stencil::NW; ++d) {
               for (uint_t a = 0; a < LatticeModel_T::Stencil::D; ++a) {
                  for (uint_t b = 0; b < LatticeModel_T::Stencil::D; ++b) {
                     if(boundaryHandling->isDomain(cell + *d) ) {
                        gradu(b, a) += velocity->get(cell + *d)[a] * real_c(d.c(b)) * real_c(0.5);
                     } else if(boundaryHandling->isDomain(cell - *d)){
                        gradu(b, a) += velocity->get(cell - *d)[a] * real_c(d.c(b)) * real_c(-0.5) + velocity->get(cell)[a] * real_c(d.c(b));
                     }
                  }
               }
            }
            auto graduT = gradu.getTranspose();

            // equation 16 from Su 2013
            relstr = stressOld->get(cell) * gradu + graduT * stressOld->get(cell);

            if(eta_p_ > real_t(0)) {
               relstr += ((gradu + graduT) * eta_p_ - stressOld->get(cell)) * inv_lambda_p_;
            }

            // equation 23 from Su 2013
            stressNew->get(cell) = stressOld->get(cell) + (stress1 * real_c(2.0) - stress2 * real_c(1.5) - stress3 * real_c(0.5)) * delta_t_ * (real_c(1.0) / pdf->getDensity(cell)) +
                                   relstr * delta_t_;
         }
      })
   }

   void swapStressBuffers( IBlock * block ){
      StressField_T * stressNew = block->getData< StressField_T >(stressId_);
      StressField_T * stressOld = block->getData< StressField_T >(stressOldId_);

      // swap pointers to old and new stress fields
      stressOld->swapDataPointers(stressNew);
   }

   void cacheVelocity( IBlock * block ){
      VelocityField_T * velocity = block->getData< VelocityField_T >(velocityId_);
      PdfField_T * pdf = block->getData< PdfField_T >(pdfId_);
      BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingId_ );

      WALBERLA_ASSERT_GREATER_EQUAL( velocity->nrOfGhostLayers(), 1 );

      // update velocity field for all fluid cells
      WALBERLA_FOR_ALL_CELLS_XYZ(velocity,{
         Cell cell(x, y, z);
         if( boundaryHandling->isDomain(cell) ) {
            velocity->get( cell ) = pdf->getVelocity(x, y, z);
         }
      })
   }

   void calculateForces( IBlock * block ) {
      using namespace stencil;
      StressField_T * stress = block->getData< StressField_T >(stressId_);
      ForceField_T * force = block->getData< ForceField_T >(forceId_);
      BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingId_ );

      WALBERLA_ASSERT_GREATER_EQUAL( stress->nrOfGhostLayers(), 1);

      WALBERLA_FOR_ALL_CELLS_XYZ(force, {
         Cell cell(x, y, z);
         uint_t k = 0;
         Vector3<real_t> f = Vector3<real_t>(0.0);

         // calculate force from finite difference divergence of extra stress in 3d or 2d
         if (boundaryHandling->isDomain(cell)) {
            for (auto d = LatticeModel_T::Stencil::beginNoCenter(); d.direction() != NW; ++d) {
               for (uint_t i = 0; i < LatticeModel_T::Stencil::D; ++i) {
                  if (d.direction() == E || d.direction() == W) {
                     k = 0;
                  } else if (d.direction() == N || d.direction() == S) {
                     k = 1;
                  } else {
                     k = 2;
                  }
                  if (boundaryHandling->isDomain(cell + *d)) {
                     f[i] += stress->get(cell + *d)(k, i) * real_c(d.c(k));
                  } else if(boundaryHandling->isDomain(cell - *d)){
                     f[i] += -stress->get(cell - *d)(k, i) * real_c(d.c(k)) + stress->get(cell)(k, i) * real_c(d.c(k)) * real_c(2.0);
                  }
               }
            }
            force->get(x, y, z) = f*real_c(0.5);
         }
      })
   }

   void operator()() {
      for( auto it = blocks_->begin(); it != blocks_->end(); ++it ) {
         auto block = it.get();
         if (executionCount_ % uint_c(delta_t_) == 0) {
            swapStressBuffers(block);
            cacheVelocity(block);
         }
      }
      communicateVelocities_();
      for( auto it = blocks_->begin(); it != blocks_->end(); ++it ){
         auto block = it.get();
         if (executionCount_ % uint_c(delta_t_) == 0)
         {
            calculateStresses(block);
         }
      }
      communicateStress_();
      for( auto it = blocks_->begin(); it != blocks_->end(); ++it ){
         auto block = it.get();
         calculateForces(block);
      }
      executionCount_++;
   }
private:
   Blocks_T blocks_;
   BlockDataID forceId_, pdfId_, boundaryHandlingId_, stressId_, stressOldId_, velocityId_;
   const real_t inv_lambda_p_, eta_p_, delta_t_;
   uint_t executionCount_;
   blockforest::communication::UniformBufferedScheme< typename LatticeModel_T::CommunicationStencil > communicateStress_, communicateVelocities_;
   bool compressibleFlag_;
};
        
} // namespace viscoelastic
} // namespace lbm
} // namespace walberla
