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
//! \file VtkWriter.h
//! \ingroup free_surface
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Free surface-specific VTK writer function.
//
//======================================================================================================================

#include "blockforest/communication/UniformBufferedScheme.h"

#include "field/adaptors/AdaptorCreators.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/FlagFieldMapping.h"
#include "field/vtk/VTKWriter.h"

#include "lbm/field/Adaptors.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/Initialization.h"

#include "FlagInfo.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Add VTK output to time loop that includes all relevant free surface information. It must be configured via
 * config-file.
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FreeSurfaceBoundaryHandling_T, typename PdfField_T, typename FlagField_T,
          typename ScalarField_T, typename VectorField_T >
void addVTKOutput(const std::weak_ptr< StructuredBlockForest >& blockForestPtr, SweepTimeloop& timeloop,
                  const std::weak_ptr< Config >& configPtr,
                  const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo, const BlockDataID& pdfFieldID,
                  const BlockDataID& flagFieldID, const BlockDataID& fillFieldID, const BlockDataID& forceFieldID,
                  const BlockDataID& curvatureFieldID, const BlockDataID& normalFieldID,
                  const BlockDataID& obstacleNormalFieldID)
{
   const auto blockForest = blockForestPtr.lock();
   WALBERLA_CHECK_NOT_NULLPTR(blockForest);

   const auto config = configPtr.lock();
   WALBERLA_CHECK_NOT_NULLPTR(config);

   // add various adaptors (for simplified access to macroscopic quantities)
   const BlockDataID densityAdaptorID = field::addFieldAdaptor< typename lbm::Adaptor< LatticeModel_T >::Density >(
      blockForest, pdfFieldID, "DensityAdaptor");
   const BlockDataID velocityAdaptorID =
      field::addFieldAdaptor< typename lbm::Adaptor< LatticeModel_T >::VelocityVector >(blockForest, pdfFieldID,
                                                                                        "VelocityVectorAdaptor");
   // define VTK output (see src/vtk/Initialization.cpp, line 574 for usage)
   const auto vtkConfigFunc = [&](std::vector< std::shared_ptr< vtk::BlockCellDataWriterInterface > >& writers,
                                  std::map< std::string, vtk::VTKOutput::CellFilter >& filters,
                                  std::map< std::string, vtk::VTKOutput::BeforeFunction >& beforeFuncs) {
      using field::VTKWriter;

      // add fields to VTK output
      writers.push_back(std::make_shared< VTKWriter< typename lbm::Adaptor< LatticeModel_T >::VelocityVector > >(
         velocityAdaptorID, "velocity"));
      writers.push_back(std::make_shared< VTKWriter< typename lbm::Adaptor< LatticeModel_T >::Density > >(
         densityAdaptorID, "density"));
      writers.push_back(std::make_shared< VTKWriter< PdfField_T, float > >(pdfFieldID, "pdf"));
      writers.push_back(std::make_shared< VTKWriter< FlagField_T, float > >(flagFieldID, "flag"));
      writers.push_back(std::make_shared< VTKWriter< ScalarField_T, float > >(fillFieldID, "fill_level"));
      writers.push_back(std::make_shared< VTKWriter< ScalarField_T, float > >(curvatureFieldID, "curvature"));
      writers.push_back(std::make_shared< VTKWriter< VectorField_T, float > >(normalFieldID, "normal"));
      writers.push_back(
         std::make_shared< VTKWriter< VectorField_T, float > >(obstacleNormalFieldID, "obstacle_normal"));
      writers.push_back(std::make_shared< VTKWriter< VectorField_T, float > >(forceFieldID, "force"));

      // map flagIDs to integer values
      const auto flagMapper =
         std::make_shared< field::FlagFieldMapping< FlagField_T, walberla::uint_t > >(flagFieldID, "mapped_flag");
      flagMapper->addMapping(flagIDs::liquidFlagID, uint_c(1));
      flagMapper->addMapping(flagIDs::interfaceFlagID, uint_c(2));
      flagMapper->addMapping(flagIDs::gasFlagID, uint_c(3));
      flagMapper->addMapping(FreeSurfaceBoundaryHandling_T::noSlipFlagID, uint_c(4));
      flagMapper->addMapping(FreeSurfaceBoundaryHandling_T::freeSlipFlagID, uint_c(6));
      flagMapper->addMapping(FreeSurfaceBoundaryHandling_T::ubbFlagID, uint_c(6));
      flagMapper->addMapping(FreeSurfaceBoundaryHandling_T::ubbInflowFlagID, uint_c(7));
      flagMapper->addMapping(FreeSurfaceBoundaryHandling_T::pressureFlagID, uint_c(8));
      flagMapper->addMapping(FreeSurfaceBoundaryHandling_T::pressureOutflowFlagID, uint_c(9));
      flagMapper->addMapping(FreeSurfaceBoundaryHandling_T::outletFlagID, uint_c(10));

      writers.push_back(flagMapper);

      // filter for writing only liquid and interface cells to VTK
      auto liquidInterfaceFilter = field::FlagFieldCellFilter< FlagField_T >(flagFieldID);
      liquidInterfaceFilter.addFlag(flagIDs::liquidFlagID);
      liquidInterfaceFilter.addFlag(flagIDs::interfaceFlagID);
      filters["liquidInterfaceFilter"] = liquidInterfaceFilter;

      // communicate fields to update the ghost layer
      auto preVTKComm = blockforest::communication::UniformBufferedScheme< stencil::D3Q27 >(blockForest);
      preVTKComm.addPackInfo(std::make_shared< field::communication::PackInfo< PdfField_T > >(pdfFieldID));
      preVTKComm.addPackInfo(std::make_shared< field::communication::PackInfo< FlagField_T > >(flagFieldID));
      preVTKComm.addPackInfo(std::make_shared< field::communication::PackInfo< ScalarField_T > >(fillFieldID));
      preVTKComm.addPackInfo(std::make_shared< field::communication::PackInfo< ScalarField_T > >(curvatureFieldID));
      preVTKComm.addPackInfo(std::make_shared< field::communication::PackInfo< VectorField_T > >(normalFieldID));
      preVTKComm.addPackInfo(
         std::make_shared< field::communication::PackInfo< VectorField_T > >(obstacleNormalFieldID));

      beforeFuncs["ghost_layer_synchronization"] = preVTKComm;

      // set velocity and density to zero in obstacle and gas cells (only for visualization purposes); the PDF values in
      // these cells are not important and thus not set during the simulation;
      // only enable this functionality if the non-liquid and non-interface cells are not excluded anyway
      const auto vtkConfigBlock        = config->getOneBlock("VTK");
      const auto fluidFieldConfigBlock = vtkConfigBlock.getBlock("fluid_field");
      if (fluidFieldConfigBlock)
      {
         auto inclusionFiltersConfigBlock = fluidFieldConfigBlock.getBlock("inclusion_filters");

         // liquidInterfaceFilter limits VTK-output to only liquid and interface cells
         if (!inclusionFiltersConfigBlock.isDefined("liquidInterfaceFilter"))
         {
            class ZeroSetter
            {
             public:
               ZeroSetter(const weak_ptr< StructuredBlockForest >& blockForest, const BlockDataID& pdfFieldID,
                          const ConstBlockDataID& flagFieldID,
                          const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo)
                  : blockForest_(blockForest), pdfFieldID_(pdfFieldID), flagFieldID_(flagFieldID), flagInfo_(flagInfo)
               {}

               void operator()()
               {
                  auto blockForest = blockForest_.lock();
                  WALBERLA_CHECK_NOT_NULLPTR(blockForest);

                  for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
                  {
                     PdfField_T* const pdfField         = blockIt->template getData< PdfField_T >(pdfFieldID_);
                     const FlagField_T* const flagField = blockIt->template getData< const FlagField_T >(flagFieldID_);
                     WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(pdfField, uint_c(1), {
                        const typename PdfField_T::Ptr pdfFieldPtr(*pdfField, x, y, z);
                        const typename FlagField_T::ConstPtr flagFieldPtr(*flagField, x, y, z);

                        if (flagInfo_.isGas(*flagFieldPtr) || flagInfo_.isObstacle(*flagFieldPtr))
                        {
                           pdfField->setDensityAndVelocity(pdfFieldPtr.cell(), Vector3< real_t >(real_c(0)),
                                                           real_c(1.0));
                        }
                     }) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
                  }
               }

             private:
               weak_ptr< StructuredBlockForest > blockForest_;
               BlockDataID pdfFieldID_;
               ConstBlockDataID flagFieldID_;
               typename FreeSurfaceBoundaryHandling_T::FlagInfo_T flagInfo_;
            };

            beforeFuncs["gas_cell_zero_setter"] = ZeroSetter(blockForest, pdfFieldID, flagFieldID, flagInfo);
         }
      }
   };

   // add VTK output to timeloop
   std::map< std::string, vtk::SelectableOutputFunction > vtkOutputFunctions;
   vtk::initializeVTKOutput(vtkOutputFunctions, vtkConfigFunc, blockForest, config);
   for (auto output = vtkOutputFunctions.begin(); output != vtkOutputFunctions.end(); ++output)
   {
      timeloop.addFuncBeforeTimeStep(output->second.outputFunction, std::string("VTK: ") + output->first,
                                     output->second.requiredGlobalStates, output->second.incompatibleGlobalStates);
   }
}

} // namespace free_surface
} // namespace walberla