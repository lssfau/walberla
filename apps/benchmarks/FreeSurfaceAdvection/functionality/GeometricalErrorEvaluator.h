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
//! \file GeometricalErrorEvaluator.h
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Compute the relative geometrical error in free-surface advection test cases.
//======================================================================================================================

#include "blockforest/StructuredBlockForest.h"

#include "core/DataTypes.h"

#include "domain_decomposition/BlockDataID.h"

#include "field/iterators/IteratorMacros.h"

namespace walberla
{
namespace free_surface
{
template< typename FreeSurfaceBoundaryHandling_T, typename FlagField_T, typename ScalarField_T >
class GeometricalErrorEvaluator
{
 public:
   GeometricalErrorEvaluator(const std::weak_ptr< StructuredBlockForest >& blockForest,
                             const std::weak_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling,
                             const ConstBlockDataID& initialfillFieldID, const ConstBlockDataID& fillFieldID,
                             uint_t frequency, const std::shared_ptr< real_t >& geometricalError)
      : blockForest_(blockForest), freeSurfaceBoundaryHandling_(freeSurfaceBoundaryHandling),
        initialFillFieldID_(initialfillFieldID), fillFieldID_(fillFieldID), frequency_(frequency),
        geometricalError_(geometricalError), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      if (frequency_ == uint_c(0)) { return; }

      auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      auto freeSurfaceBoundaryHandling = freeSurfaceBoundaryHandling_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(freeSurfaceBoundaryHandling);

      if (executionCounter_ == uint_c(0))
      {
         computeInitialFillLevelSum(blockForest, freeSurfaceBoundaryHandling);
         computeError(blockForest, freeSurfaceBoundaryHandling);
      }
      else
      {
         // only evaluate in given frequencies
         if (executionCounter_ % frequency_ == uint_c(0)) { computeError(blockForest, freeSurfaceBoundaryHandling); }
      }

      ++executionCounter_;
   }

   void computeInitialFillLevelSum(
      const std::shared_ptr< const StructuredBlockForest >& blockForest,
      const std::shared_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling)
   {
      const BlockDataID flagFieldID = freeSurfaceBoundaryHandling->getFlagFieldID();
      const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         const ScalarField_T* const initialfillField = blockIt->getData< const ScalarField_T >(initialFillFieldID_);
         const FlagField_T* const flagField          = blockIt->getData< const FlagField_T >(flagFieldID);

         // avoid OpenMP here because initialFillLevelSum_ is a class member and not a regular variable
         WALBERLA_FOR_ALL_CELLS_OMP(initialfillFieldIt, initialfillField, flagFieldIt, flagField, omp critical, {
            if (flagInfo.isInterface(flagFieldIt) || flagInfo.isLiquid(flagFieldIt))
            {
               initialFillLevelSum_ += *initialfillFieldIt;
            }
         }) // WALBERLA_FOR_ALL_CELLS_OMP
      }

      mpi::allReduceInplace< real_t >(initialFillLevelSum_, mpi::SUM);
   }

   void computeError(const std::shared_ptr< const StructuredBlockForest >& blockForest,
                     const std::shared_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling)
   {
      const BlockDataID flagFieldID = freeSurfaceBoundaryHandling->getFlagFieldID();
      const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

      real_t geometricalError = real_c(0);

      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         const ScalarField_T* const initialfillField = blockIt->getData< const ScalarField_T >(initialFillFieldID_);
         const ScalarField_T* const fillField        = blockIt->getData< const ScalarField_T >(fillFieldID_);
         const FlagField_T* const flagField          = blockIt->getData< const FlagField_T >(flagFieldID);

         WALBERLA_FOR_ALL_CELLS_OMP(initialfillFieldIt, initialfillField, fillFieldIt, fillField, flagFieldIt,
                                    flagField, omp parallel for schedule(static) reduction(+:geometricalError), {
            if (flagInfo.isInterface(flagFieldIt) || flagInfo.isLiquid(flagFieldIt))
            {
               geometricalError += real_c(std::abs(*initialfillFieldIt - *fillFieldIt));
            }
         }) // WALBERLA_FOR_ALL_CELLS_OMP
      }

      mpi::allReduceInplace< real_t >(geometricalError, mpi::SUM);

      // compute L1 norms
      *geometricalError_ = geometricalError / initialFillLevelSum_;
   }

 private:
   std::weak_ptr< StructuredBlockForest > blockForest_;
   std::weak_ptr< const FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling_;
   ConstBlockDataID initialFillFieldID_;
   ConstBlockDataID fillFieldID_;
   uint_t frequency_;
   std::shared_ptr< real_t > geometricalError_;

   uint_t executionCounter_;
   real_t initialFillLevelSum_ = real_c(0);
}; // class GeometricalErrorEvaluator

} // namespace free_surface
} // namespace walberla