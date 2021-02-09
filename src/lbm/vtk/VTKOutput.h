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
//! \file VTKOutput.h
//! \ingroup lbm
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Default LBM VTK Output
//
//======================================================================================================================

#pragma once

#include "Density.h"
#include "Velocity.h"
#include "lbm/field/PdfField.h"

#include "blockforest/StructuredBlockForest.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/DataTypes.h"
#include "core/config/Config.h"

#include "domain_decomposition/BlockDataID.h"

#include "field/FlagField.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/FlagFieldMapping.h"
#include "field/vtk/VTKWriter.h"
#include "field/communication/PackInfo.h"

#include "stencil/D3Q27.h"

#include "timeloop/Timeloop.h"

#include "vtk/Initialization.h"
#include "vtk/VTKOutput.h"

#include <map>
#include <string>
#include <vector>


namespace walberla {
namespace lbm {


//**********************************************************************************************************************
/*!
* \brief Default setup for VTK in LBM simulations
*
* \ingroup lbm
*
* This functor is to be used with vtk::initializeVTKOutput.
* It adds five VTK writers: "Velocity", "VelocityMagnitude", and "Density" which write the respective quantities in
* lattice units, and "FlagField" and "MappedFlagField" which can be used to visualize the flag field.
* Additionally, it adds the cell filter "DomainFilter" which can be used to write only the domain
* (a.k.a. fluid) parts of the simulation.
* It is also possible to add a vtk::VTKOutput::BeforeFunction to synchronize the ghost layers before writing VTK output.
* This before function will be registered as "PDFGhostLayerSync".
*
* \tparam LatticeModel  The lattice model used for the simulation
* \tparam FlagFieldT    Type of the used flag field
*/
//**********************************************************************************************************************
template <typename LatticeModel, typename FlagFieldT>
class VTKOutput {

public:

   typedef typename FlagFieldT::flag_t flag_t;
   typedef std::map<FlagUID, flag_t> FlagMap;

   VTKOutput( const ConstBlockDataID & pdfField, const ConstBlockDataID & flagField,
              const FlagUID & domainFlagUID, const vtk::VTKOutput::BeforeFunction & ghostLayerSyncFunction,
              const FlagMap & flagFieldMapping = FlagMap() );

   VTKOutput( const ConstBlockDataID & pdfField, const ConstBlockDataID & flagField,
              const FlagUID & domainFlagUID,
              const FlagMap & flagFieldMapping = FlagMap() );

   virtual ~VTKOutput() = default;

   virtual void operator()( std::vector< shared_ptr<vtk::BlockCellDataWriterInterface> > & writers,
                            std::map< std::string, vtk::VTKOutput::CellFilter > & filters,
                            std::map< std::string, vtk::VTKOutput::BeforeFunction > & beforeFunctions ) const;

   static void addToTimeloop( Timeloop & timeloop, const shared_ptr< blockforest::StructuredBlockForest > & blocks,
                              const shared_ptr< Config > & config, const std::string & configBlockName,
                              const BlockDataID & pdfField, const ConstBlockDataID & flagField,
                              const FlagUID & domainFlagUID,
                              const FlagMap & flagFieldMapping = FlagMap() );

   static void addToTimeloop( Timeloop & timeloop, const shared_ptr< blockforest::StructuredBlockForest > & blocks,
                              const shared_ptr< Config > & config,
                              const BlockDataID & pdfField, const ConstBlockDataID & flagField,
                              const FlagUID & domainFlagUID,
                              const FlagMap & flagFieldMapping = FlagMap() );

protected:

   virtual void addWriters        ( std::vector< shared_ptr<vtk::BlockCellDataWriterInterface> >& writers     ) const;
   virtual void addFilters        ( std::map< std::string, vtk::VTKOutput::CellFilter >& filters              ) const;
   virtual void addBeforeFunctions( std::map< std::string, vtk::VTKOutput::BeforeFunction > & beforeFunctions ) const;

   BlockDataID pdfField_;
   ConstBlockDataID flagField_;

   FlagUID domainFlagUID_;

   vtk::VTKOutput::BeforeFunction ghostLayerSyncFunction_;

   FlagMap flagFieldMapping_;

}; // class VTKOutput



//**********************************************************************************************************************
/*!
* \ingroup lbm
*
* \param pdfField               BlockDataID of the PDF field used in the simulation
* \param flagField              BlockDataID of the flag field used in the simulation
* \param domainFlagUID          FlagUID of the domain (a.k.a fluid) flag
* \param ghostLayerSyncFunction A function used to synchronize the ghost layers before the VTK output. It will be called
*                               once every time step, before the output is written.
* \param flagFieldMapping       A mapping of FlagUIDs to flag_ts. The mapping is used to create the MappedFlagField
*                               writer. If omitted, the MappedFlagField writer will not be available.
*/
//**********************************************************************************************************************
template <typename LatticeModel, typename FlagFieldT>
VTKOutput<LatticeModel, FlagFieldT>::VTKOutput( const ConstBlockDataID & pdfField, const ConstBlockDataID & flagField,
                                                const FlagUID & domainFlagUID, const vtk::VTKOutput::BeforeFunction & ghostLayerSyncFunction,
                                                const FlagMap & flagFieldMapping /*= FlagMap()*/ )
   : pdfField_( pdfField ), flagField_( flagField ), domainFlagUID_( domainFlagUID ),
     ghostLayerSyncFunction_( ghostLayerSyncFunction ), flagFieldMapping_( flagFieldMapping )
{
}



//**********************************************************************************************************************
/*!
* \ingroup lbm
*
* \param pdfField               BlockDataID of the PDF field used in the simulation
* \param flagField              BlockDataID of the flag field used in the simulation
* \param domainFlagUID          FlagUID of the domain (a.k.a fluid) flag
* \param flagFieldMapping       A mapping of FlagUIDs to flag_ts. The mapping is used to create the MappedFlagField
*                               writer. If omitted, the MappedFlagField writer will not be available.
*/
//**********************************************************************************************************************
template <typename LatticeModel, typename FlagFieldT>
VTKOutput<LatticeModel, FlagFieldT>::VTKOutput( const ConstBlockDataID & pdfField, const ConstBlockDataID & flagField,
                                                const FlagUID & domainFlagUID,
                                                const FlagMap & flagFieldMapping /*= FlagMap*/ )
   : pdfField_( pdfField ), flagField_( flagField ), domainFlagUID_( domainFlagUID ), flagFieldMapping_( flagFieldMapping )
{
}

//**********************************************************************************************************************
/*!
* \brief see vtk::initializeVTKOutput for details
*
* \ingroup lbm
*/
//**********************************************************************************************************************
template <typename LatticeModel, typename FlagFieldT>
void VTKOutput<LatticeModel, FlagFieldT>::operator()( std::vector< shared_ptr<vtk::BlockCellDataWriterInterface> > & writers,
                                                      std::map< std::string, vtk::VTKOutput::CellFilter > & filters,
                                                      std::map< std::string, vtk::VTKOutput::BeforeFunction > &  beforeFunctions ) const
{
   addWriters( writers );
   addFilters( filters );
   addBeforeFunctions( beforeFunctions );
}



template <typename LatticeModel, typename FlagFieldT>
void VTKOutput<LatticeModel, FlagFieldT>::addToTimeloop( Timeloop & timeloop, const shared_ptr< blockforest::StructuredBlockForest > & blocks,
                                                         const shared_ptr< Config > & config, const std::string & configBlockName,
                                                         const BlockDataID & pdfField, const ConstBlockDataID & flagField,
                                                         const FlagUID & domainFlagUID,
                                                         const FlagMap & flagFieldMapping )
{
   typedef PdfField< LatticeModel > PdfField_t;

   blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > pdfGhostLayerSync( blocks );
   pdfGhostLayerSync.addPackInfo( make_shared< field::communication::PackInfo< PdfField_t > >( pdfField ) );

   VTKOutput vtkOutput( pdfField, flagField, domainFlagUID, pdfGhostLayerSync, flagFieldMapping );

   std::map< std::string, vtk::SelectableOutputFunction > vtkOutputFunctions;
   vtk::initializeVTKOutput( vtkOutputFunctions, vtkOutput, blocks, config, configBlockName );

   for( auto output = vtkOutputFunctions.begin(); output != vtkOutputFunctions.end(); ++output )
      timeloop.addFuncBeforeTimeStep( output->second.outputFunction, std::string("VTK (LBM): ") + output->first,
                                      output->second.requiredGlobalStates, output->second.incompatibleGlobalStates );
}



template <typename LatticeModel, typename FlagFieldT>
void VTKOutput<LatticeModel, FlagFieldT>::addToTimeloop( Timeloop & timeloop, const shared_ptr< blockforest::StructuredBlockForest > & blocks,
                                                         const shared_ptr< Config > & config,
                                                         const BlockDataID & pdfField, const ConstBlockDataID & flagField,
                                                         const FlagUID & domainFlagUID,
                                                         const FlagMap & flagFieldMapping )
{
   addToTimeloop( timeloop, blocks, config, "VTK", pdfField, flagField, domainFlagUID, flagFieldMapping );
}



template <typename LatticeModel, typename FlagFieldT>
void VTKOutput<LatticeModel, FlagFieldT>::addWriters( std::vector< shared_ptr<vtk::BlockCellDataWriterInterface> >& writers ) const
{
   writers.push_back( make_shared< VelocityVTKWriter< LatticeModel > >         ( pdfField_, "Velocity" ) );
   writers.push_back( make_shared< VelocityMagnitudeVTKWriter< LatticeModel > >( pdfField_, "VelocityMagnitude" ) );
   writers.push_back( make_shared< DensityVTKWriter< LatticeModel > >          ( pdfField_, "Density"  ) );

   if( !flagFieldMapping_.empty() )
      writers.push_back( walberla::make_shared< field::FlagFieldMapping< FlagFieldT, flag_t > >( flagField_, "MappedFlagField", flagFieldMapping_ ) );

   writers.push_back( make_shared< field::VTKWriter< FlagFieldT > >( flagField_, "FlagField" ) );
}



template <typename LatticeModel, typename FlagFieldT>
void VTKOutput<LatticeModel, FlagFieldT>::addFilters( std::map< std::string, vtk::VTKOutput::CellFilter >& filters ) const
{
   filters[ "DomainFilter" ] = field::FlagFieldCellFilter< FlagFieldT >( flagField_, domainFlagUID_ );
}



template <typename LatticeModel, typename FlagFieldT>
void VTKOutput<LatticeModel, FlagFieldT>::addBeforeFunctions( std::map< std::string, vtk::VTKOutput::BeforeFunction > & beforeFunctions ) const
{
   if( ghostLayerSyncFunction_ )
      beforeFunctions[ "PDFGhostLayerSync" ] = ghostLayerSyncFunction_;
}

} // namespace lbm
} // namespace walberla
