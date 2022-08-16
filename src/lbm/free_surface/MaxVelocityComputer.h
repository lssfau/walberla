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
//! \file MaxVelocityComputer.h
//! \ingroup free_surface
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Compute the maximum velocity of all liquid and interface cells (in each direction) in the system.
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"

namespace walberla
{
namespace free_surface
{
template< typename FreeSurfaceBoundaryHandling_T, typename PdfField_T, typename FlagField_T >
class MaxVelocityComputer
{
 public:
   MaxVelocityComputer(const std::weak_ptr< const StructuredBlockForest >& blockForest,
                       const std::weak_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling,
                       const ConstBlockDataID& pdfFieldID, uint_t frequency,
                       const std::shared_ptr< Vector3< real_t > >& maxVelocity)
      : blockForest_(blockForest), freeSurfaceBoundaryHandling_(freeSurfaceBoundaryHandling), pdfFieldID_(pdfFieldID),
        maxVelocity_(maxVelocity), frequency_(frequency), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      if (frequency_ == uint_c(0)) { return; }

      auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      auto freeSurfaceBoundaryHandling = freeSurfaceBoundaryHandling_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(freeSurfaceBoundaryHandling);

      if (executionCounter_ == uint_c(0)) { getMaxVelocity(blockForest, freeSurfaceBoundaryHandling); }
      else
      {
         // only evaluate in given frequencies
         if (executionCounter_ % frequency_ == uint_c(0)) { getMaxVelocity(blockForest, freeSurfaceBoundaryHandling); }
      }

      ++executionCounter_;
   }

   void getMaxVelocity(const std::shared_ptr< const StructuredBlockForest >& blockForest,
                       const std::shared_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling)
   {
      const BlockDataID flagFieldID = freeSurfaceBoundaryHandling->getFlagFieldID();
      const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

      real_t maxVelocityX = real_c(0);
      real_t maxVelocityY = real_c(0);
      real_t maxVelocityZ = real_c(0);

      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         const FlagField_T* const flagField = blockIt->template getData< const FlagField_T >(flagFieldID);
         const PdfField_T* const pdfField   = blockIt->template getData< const PdfField_T >(pdfFieldID_);

            WALBERLA_FOR_ALL_CELLS_OMP(flagFieldIt, flagField, pdfFieldIt, pdfField,
                                       omp parallel for schedule(static) reduction(max:maxVelocityX)
                                                                         reduction(max:maxVelocityY)
                                                                         reduction(max:maxVelocityZ),
            {
            if (flagInfo.isLiquid(flagFieldIt) || flagInfo.isInterface(flagFieldIt))
            {
               const Vector3< real_t > velocity = pdfField->getVelocity(pdfFieldIt.cell());

               if (velocity[0] > maxVelocityX) { maxVelocityX = velocity[0]; }
               if (velocity[1] > maxVelocityY) { maxVelocityY = velocity[1]; }
               if (velocity[2] > maxVelocityZ) { maxVelocityZ = velocity[2]; }
            }
            }) // WALBERLA_FOR_ALL_CELLS_OMP
      }

      Vector3< real_t > maxVelocity(maxVelocityX, maxVelocityY, maxVelocityZ);
      mpi::allReduceInplace< real_t >(maxVelocity, mpi::MAX);

      *maxVelocity_ = maxVelocity;
   };

 private:
   std::weak_ptr< const StructuredBlockForest > blockForest_;
   std::weak_ptr< const FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling_;

   const ConstBlockDataID pdfFieldID_;

   std::shared_ptr< Vector3< real_t > > maxVelocity_;

   uint_t frequency_;
   uint_t executionCounter_;
}; // class MaxVelocityComputer

} // namespace free_surface
} // namespace walberla