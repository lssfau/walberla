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
//! \file TotalMassComputer.h
//! \ingroup free_surface
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Compute the total mass of the system (including mass from the excessMassField).
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/StringUtility.h"

#include "lbm/free_surface/dynamics/SurfaceDynamicsHandler.h"
#include "lbm/free_surface/surface_geometry/SurfaceGeometryHandler.h"
#include "lbm/free_surface/surface_geometry/Utility.h"
#include "lbm/lattice_model/D3Q19.h"

namespace walberla
{
namespace free_surface
{
template< typename FreeSurfaceBoundaryHandling_T, typename PdfField_T, typename FlagField_T, typename ScalarField_T >
class TotalMassComputer
{
 public:
   TotalMassComputer(const std::weak_ptr< const StructuredBlockForest >& blockForest,
                     const std::weak_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling,
                     const ConstBlockDataID& pdfFieldID, const ConstBlockDataID& fillFieldID, uint_t frequency,
                     const std::shared_ptr< real_t >& totalMass)
      : blockForest_(blockForest), freeSurfaceBoundaryHandling_(freeSurfaceBoundaryHandling), pdfFieldID_(pdfFieldID),
        fillFieldID_(fillFieldID), totalMass_(totalMass), frequency_(frequency), executionCounter_(uint_c(0))
   {}

   TotalMassComputer(const std::weak_ptr< const StructuredBlockForest >& blockForest,
                     const std::weak_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling,
                     const ConstBlockDataID& pdfFieldID, const ConstBlockDataID& fillFieldID,
                     const ConstBlockDataID& excessMassFieldID, uint_t frequency,
                     const std::shared_ptr< real_t >& totalMass, const std::shared_ptr< real_t >& excessMass)
      : blockForest_(blockForest), freeSurfaceBoundaryHandling_(freeSurfaceBoundaryHandling), pdfFieldID_(pdfFieldID),
        fillFieldID_(fillFieldID), excessMassFieldID_(excessMassFieldID), totalMass_(totalMass),
        excessMass_(excessMass), frequency_(frequency), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      if (frequency_ == uint_c(0)) { return; }

      auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      auto freeSurfaceBoundaryHandling = freeSurfaceBoundaryHandling_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(freeSurfaceBoundaryHandling);

      // only evaluate in given frequencies
      if (executionCounter_ % frequency_ == uint_c(0) || executionCounter_ == uint_c(0))
      {
         computeMass(blockForest, freeSurfaceBoundaryHandling);
      }

      ++executionCounter_;
   }

   void computeMass(const std::shared_ptr< const StructuredBlockForest >& blockForest,
                    const std::shared_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling)
   {
      const BlockDataID flagFieldID = freeSurfaceBoundaryHandling->getFlagFieldID();
      const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

      real_t mass       = real_c(0);
      real_t excessMass = real_c(0);

      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         const FlagField_T* const flagField   = blockIt->template getData< const FlagField_T >(flagFieldID);
         const PdfField_T* const pdfField     = blockIt->template getData< const PdfField_T >(pdfFieldID_);
         const ScalarField_T* const fillField = blockIt->template getData< const ScalarField_T >(fillFieldID_);

         // if provided, also consider mass stored in excessMassField
         if (excessMassFieldID_ != ConstBlockDataID())
         {
            const ScalarField_T* const excessMassField =
               blockIt->template getData< const ScalarField_T >(excessMassFieldID_);

            WALBERLA_FOR_ALL_CELLS_OMP(flagFieldIt, flagField, pdfFieldIt, pdfField, fillFieldIt, fillField,
                                       excessMassFieldIt, excessMassField,
                                       omp parallel for schedule(static) reduction(+:mass),
            {
               if (flagInfo.isLiquid(flagFieldIt) || flagInfo.isInterface(flagFieldIt))
               {
                  const real_t density = pdfField->getDensity(pdfFieldIt.cell());
                  mass += *fillFieldIt * density + *excessMassFieldIt;

                  if (excessMass_ != nullptr) { excessMass += *excessMassFieldIt; }
               }
            }) // WALBERLA_FOR_ALL_CELLS_OMP
         }
         else
         {
            WALBERLA_FOR_ALL_CELLS_OMP(flagFieldIt, flagField, pdfFieldIt, pdfField, fillFieldIt, fillField,
                                       omp parallel for schedule(static) reduction(+:mass),
            {
               if (flagInfo.isLiquid(flagFieldIt) || flagInfo.isInterface(flagFieldIt))
               {
                  const real_t density = pdfField->getDensity(pdfFieldIt.cell());
                  mass += *fillFieldIt * density;
               }
            }) // WALBERLA_FOR_ALL_CELLS_OMP
         }
      }

      mpi::allReduceInplace< real_t >(mass, mpi::SUM);
      *totalMass_ = mass;

      if (excessMass_ != nullptr)
      {
         mpi::allReduceInplace< real_t >(excessMass, mpi::SUM);
         *excessMass_ = excessMass;
      }
   };

 private:
   std::weak_ptr< const StructuredBlockForest > blockForest_;
   std::weak_ptr< const FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling_;

   const ConstBlockDataID pdfFieldID_;
   const ConstBlockDataID fillFieldID_;
   const ConstBlockDataID excessMassFieldID_ = ConstBlockDataID();

   std::shared_ptr< real_t > totalMass_;
   std::shared_ptr< real_t > excessMass_ = nullptr; // mass stored in the excessMassField

   uint_t frequency_;
   uint_t executionCounter_;
}; // class TotalMassComputer

} // namespace free_surface
} // namespace walberla