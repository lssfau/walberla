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
//! \file Initialization.h
//! \ingroup vtk
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "VTKOutput.h"
#include "core/DataTypes.h"
#include "core/config/Config.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include <functional>
#include <string>


namespace walberla {
namespace vtk {


// For documentation see the documentation of function "initializeVTKOutput" in Initialization.cpp


typedef std::function< void () > OutputFunction;

struct SelectableOutputFunction {

   SelectableOutputFunction() = default;
   SelectableOutputFunction( OutputFunction of, const Set<SUID>& rgs, const Set<SUID>& igs ) :
      outputFunction( of ), requiredGlobalStates( rgs ), incompatibleGlobalStates( igs ) {}

   OutputFunction outputFunction;
   Set<SUID> requiredGlobalStates;
   Set<SUID> incompatibleGlobalStates;
};

void initializeVTKOutput( std::map< std::string, SelectableOutputFunction > & outputFunctions,
                          const shared_ptr< const StructuredBlockStorage > & storage,
                          const shared_ptr< Config > & config, const std::string & configBlockName,
                          const std::vector< shared_ptr< BlockCellDataWriterInterface > > & writers,
                          const std::map< std::string, VTKOutput::CellFilter > & filters,
                          const std::map< std::string, VTKOutput::BeforeFunction > & beforeFunctions );

void initializeVTKOutput( std::map< std::string, SelectableOutputFunction > & outputFunctions,
                          const shared_ptr< const StructuredBlockStorage > & storage,
                          const Config::BlockHandle & parentBlockHandle, const std::string & configBlockName,
                          const std::vector< shared_ptr< BlockCellDataWriterInterface > > & writers,
                          const std::map< std::string, VTKOutput::CellFilter > & filters,
                          const std::map< std::string, VTKOutput::BeforeFunction > & beforeFunctions );

void initializeVTKOutput( std::map< std::string, SelectableOutputFunction > & outputFunctions,
                          const shared_ptr< const StructuredBlockStorage > & storage, const shared_ptr< Config > & config,
                          const std::vector< shared_ptr< BlockCellDataWriterInterface > > & writers,
                          const std::map< std::string, VTKOutput::CellFilter > & filters,
                          const std::map< std::string, VTKOutput::BeforeFunction > & beforeFunctions );

void initializeVTKOutput( std::map< std::string, SelectableOutputFunction > & outputFunctions,
                          const shared_ptr< const StructuredBlockStorage > & storage, const Config::BlockHandle & parentBlockHandle,
                          const std::vector< shared_ptr< BlockCellDataWriterInterface > > & writers,
                          const std::map< std::string, VTKOutput::CellFilter > & filters,
                          const std::map< std::string, VTKOutput::BeforeFunction > & beforeFunctions );

typedef std::function< void ( std::vector< shared_ptr< BlockCellDataWriterInterface > > & writers,
                                std::map< std::string, VTKOutput::CellFilter > &            filters,
                                std::map< std::string, VTKOutput::BeforeFunction > &        beforeFunctions ) > RegisterVTKOutputFunction;

void initializeVTKOutput( std::map< std::string, SelectableOutputFunction > & outputFunctions, const RegisterVTKOutputFunction& registerVTKOutputFunction,
                          const shared_ptr< const StructuredBlockStorage > & storage, const shared_ptr< Config > & config,
                          const std::string & configBlockName = std::string("VTK") );

void initializeVTKOutput( std::map< std::string, SelectableOutputFunction > & outputFunctions, const RegisterVTKOutputFunction& registerVTKOutputFunction,
                          const shared_ptr< const StructuredBlockStorage > & storage, const Config::BlockHandle & parentBlockHandle,
                          const std::string & configBlockName = std::string("VTK") );



} // namespace vtk
} // namespace walberla

