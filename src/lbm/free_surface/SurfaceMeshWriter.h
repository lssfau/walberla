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
//! \file SurfaceMeshWriter.h
//! \ingroup free_surface
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Free surface-specific class for writing the free surface as triangle mesh.
//
//======================================================================================================================

#include "blockforest/StructuredBlockForest.h"

#include "core/Filesystem.h"

#include "field/AddToStorage.h"

#include "geometry/mesh/TriangleMeshIO.h"

#include "postprocessing/FieldToSurfaceMesh.h"

#include "timeloop/SweepTimeloop.h"

namespace walberla
{
namespace free_surface
{
namespace abortIfNullptr
{

// helper function to check validity of std::weak_ptr in constructors' initializer list; WALBERLA_CHECK_NOT_NULLPTR()
// does not work there, because the macro terminates with ";"
template< typename T >
void abortIfNullptr(const std::weak_ptr< T >& weakPointer)
{
   if (weakPointer.lock() == nullptr) { WALBERLA_ABORT("Weak pointer has expired."); }
}
} // namespace abortIfNullptr

/***********************************************************************************************************************
 * Write free surface as triangle mesh.
 *
 * Internally, a clone of the fill level field is stored and all cells not marked as liquidInterfaceGasFlagIDSet are set
 * to "obstacleFillLevel" in the cloned field. This is done to avoid writing e.g. obstacle cells that were possibly
 * assigned a fill level of 1 to not make them detect as gas cells in the bubble model.
 **********************************************************************************************************************/
template< typename ScalarField_T, typename FlagField_T >
class SurfaceMeshWriter
{
 public:
   SurfaceMeshWriter(const std::weak_ptr< StructuredBlockForest >& blockForest, const ConstBlockDataID& fillFieldID,
                     const ConstBlockDataID& flagFieldID, const Set< FlagUID >& liquidInterfaceGasFlagIDSet,
                     real_t obstacleFillLevel, uint_t writeFrequency, const std::string& baseFolder)
      : blockForest_(blockForest), fillFieldID_(fillFieldID), flagFieldID_(flagFieldID),
        liquidInterfaceGasFlagIDSet_(liquidInterfaceGasFlagIDSet), obstacleFillLevel_(obstacleFillLevel),
        writeFrequency_(writeFrequency), baseFolder_(baseFolder), executionCounter_(uint_c(0)),
        fillFieldCloneID_(
           (abortIfNullptr::abortIfNullptr(blockForest),
            field::addCloneToStorage< ScalarField_T >(blockForest_.lock(), fillFieldID_, "Fill level field clone")))
   {}

   // config block must be named "MeshOutputParameters"
   SurfaceMeshWriter(const std::weak_ptr< StructuredBlockForest >& blockForest, const ConstBlockDataID& fillFieldID,
                     const ConstBlockDataID& flagFieldID, const Set< FlagUID >& liquidInterfaceGasFlagIDSet,
                     real_t obstacleFillLevel, const std::weak_ptr< Config >& config)
      : SurfaceMeshWriter(
           blockForest, fillFieldID, flagFieldID, liquidInterfaceGasFlagIDSet, obstacleFillLevel,
           (abortIfNullptr::abortIfNullptr(config),
            config.lock()->getOneBlock("MeshOutputParameters").getParameter< uint_t >("writeFrequency")),
           (abortIfNullptr::abortIfNullptr(config),
            config.lock()->getOneBlock("MeshOutputParameters").getParameter< std::string >("baseFolder")))
   {}

   void operator()()
   {
      if (writeFrequency_ == uint_c(0)) { return; }

      if (executionCounter_ == uint_c(0))
      {
         createBaseFolder();
         writeMesh();
      }
      else { writeMesh(); }

      ++executionCounter_;
   }

 private:
   void createBaseFolder() const
   {
      WALBERLA_ROOT_SECTION()
      {
         const filesystem::path basePath(baseFolder_);
         if (filesystem::exists(basePath)) { filesystem::remove_all(basePath); }
         filesystem::create_directories(basePath);
      }
      WALBERLA_MPI_BARRIER();
   }

   void writeMesh()
   {
      // only write mesh in given frequency
      if (executionCounter_ % writeFrequency_ != uint_c(0)) { return; }

      // rank=0 is just an arbitrary choice here
      const int targetRank = 0;

      auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      // update clone of fill level field and set fill level of all non-liquid, -interface, or -gas cells to zero
      updateFillFieldClone(blockForest);

      const auto surfaceMesh = postprocessing::realFieldToSurfaceMesh< ScalarField_T >(
         blockForest, fillFieldCloneID_, real_c(0.5), uint_c(0), true, targetRank, MPI_COMM_WORLD);

      WALBERLA_EXCLUSIVE_WORLD_SECTION(targetRank)
      {
         geometry::writeMesh(baseFolder_ + "/" + "simulation_step_" + std::to_string(executionCounter_) + ".obj",
                             *surfaceMesh);
      }
   }

   // update clone of fill level field and set fill level of all non-liquid, -interface, or -gas cells to zero;
   // explicitly use shared_ptr instead of weak_ptr to avoid checking the latter's validity (is done in writeMesh()
   // already)
   void updateFillFieldClone(const shared_ptr< StructuredBlockForest >& blockForest)
   {
      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         ScalarField_T* const fillFieldClone  = blockIt->template getData< ScalarField_T >(fillFieldCloneID_);
         const ScalarField_T* const fillField = blockIt->template getData< const ScalarField_T >(fillFieldID_);
         const FlagField_T* const flagField   = blockIt->template getData< const FlagField_T >(flagFieldID_);

         const auto liquidInterfaceGasFlagMask = flagField->getMask(liquidInterfaceGasFlagIDSet_);

         WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(fillFieldClone, fillField->nrOfGhostLayers(), {
            const typename ScalarField_T::Ptr fillFieldClonePtr(*fillFieldClone, x, y, z);
            const typename ScalarField_T::ConstPtr fillFieldPtr(*fillField, x, y, z);
            const typename FlagField_T::ConstPtr flagFieldPtr(*flagField, x, y, z);

            // set fill level to zero in every non-liquid, -interface, or -gas cell
            if (!isPartOfMaskSet(flagFieldPtr, liquidInterfaceGasFlagMask)) { *fillFieldClonePtr = obstacleFillLevel_; }
            else
            {
               // copy fill level from fill level field
               *fillFieldClonePtr = *fillFieldPtr;
            }
         }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
      }
   }

   std::weak_ptr< StructuredBlockForest > blockForest_;
   ConstBlockDataID fillFieldID_;
   ConstBlockDataID flagFieldID_;
   Set< FlagUID > liquidInterfaceGasFlagIDSet_;
   real_t obstacleFillLevel_;
   uint_t writeFrequency_;
   std::string baseFolder_;
   uint_t executionCounter_;

   BlockDataID fillFieldCloneID_;
}; // class SurfaceMeshWriter
} // namespace free_surface
} // namespace walberla