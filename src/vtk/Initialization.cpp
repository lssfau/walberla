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
//! \file Initialization.cpp
//! \ingroup vtk
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "ChainedFilter.h"
#include "Initialization.h"

#include "core/Abort.h"
#include "core/logging/Logging.h"
#include "core/StringUtility.h"

#include <functional>


namespace walberla {
namespace vtk {



template< typename T >
static void splitVector( T& x, T& y, T& z, const Config::BlockHandle& bb, const std::string& vertex, const std::string& errorMsg )
{
   std::vector< std::string > coordinates;
   std::string vector = bb.getParameter< std::string >( vertex );
   coordinates = string_split( vector, "<,> \t" );

   coordinates.erase( std::remove_if( coordinates.begin(), coordinates.end(), std::bind( &std::string::empty,  std::placeholders::_1 ) ), coordinates.end() );

   if( coordinates.size() != 3 )
      WALBERLA_ABORT( errorMsg );

   x = string_to_num< T >( coordinates[0] );
   y = string_to_num< T >( coordinates[1] );
   z = string_to_num< T >( coordinates[2] );
}



static std::vector< std::string > splitList( const std::string& string )
{
   std::vector< std::string > list;

   list = string_split( string, ", \t" );
   list.erase( std::remove_if( list.begin(), list.end(), std::bind( &std::string::empty,  std::placeholders::_1 ) ), list.end() );

   return list;
}



static void addStates( Set<SUID>& set, const std::string& string )
{
   std::vector< std::string > states;
   states = string_split( string, ", \t" );
   states.erase( std::remove_if( states.begin(), states.end(), std::bind( &std::string::empty,  std::placeholders::_1 ) ), states.end() );

   for( auto it = states.begin(); it != states.end(); ++it )
      set += SUID( *it );
}



void initializeVTKOutput( std::map< std::string, SelectableOutputFunction > & outputFunctions,
                          const shared_ptr< const StructuredBlockStorage > & storage,
                          const shared_ptr< Config > & config, const std::string & configBlockName,
                          const std::vector< shared_ptr< BlockCellDataWriterInterface > > & writers,
                          const std::map< std::string, VTKOutput::CellFilter > & filters,
                          const std::map< std::string, VTKOutput::BeforeFunction > & beforeFunctions )
{
   if( !!config )
      initializeVTKOutput( outputFunctions, storage, config->getGlobalBlock(), configBlockName, writers, filters, beforeFunctions );
}



struct CaseInsensitiveCompare {
   bool operator()( const std::string& lhs, const std::string& rhs ) const { return ( string_icompare(lhs, rhs) < 0 ); }
}; // only required in 'initializeVTKOutput' below

void initializeVTKOutput( std::map< std::string, SelectableOutputFunction > & outputFunctions,
                          const shared_ptr< const StructuredBlockStorage > & storage,
                          const Config::BlockHandle & parentBlockHandle, const std::string & configBlockName,
                          const std::vector< shared_ptr< BlockCellDataWriterInterface > > & _writers,
                          const std::map< std::string, VTKOutput::CellFilter > & _filters,
                          const std::map< std::string, VTKOutput::BeforeFunction > & _beforeFunctions )
{
   if( !parentBlockHandle )
      WALBERLA_ABORT("Invalid Argument: parentBlockHandle not valid!");

   std::map< std::string, shared_ptr< BlockCellDataWriterInterface >, CaseInsensitiveCompare > writers;
   std::map< std::string, VTKOutput::CellFilter,                      CaseInsensitiveCompare > globalFilters;
   std::map< std::string, VTKOutput::BeforeFunction,                  CaseInsensitiveCompare > beforeFunctions;

   for( auto writer = _writers.begin(); writer != _writers.end(); ++writer )
   {
      if( writers.find( (*writer)->identifier() ) != writers.end() )
         WALBERLA_ABORT( "There are at least two block data writers with identifier \"" << (*writer)->identifier() <<
                         "\" (test is case insensitive!).\nEvery writer must have a unique identifier!" );
      writers[ (*writer)->identifier() ] = *writer;
   }

   for( auto filter = _filters.begin(); filter != _filters.end(); ++filter )
   {
      if( globalFilters.find( filter->first ) != globalFilters.end() )
         WALBERLA_ABORT( "There are at least two cell filters with identifier \"" << filter->first <<
                         "\" (test is case insensitive!).\nEvery filter must have a unique identifier!" );
      globalFilters[ filter->first ] = filter->second;
   }

   for( auto function = _beforeFunctions.begin(); function != _beforeFunctions.end(); ++function )
   {
      if( beforeFunctions.find( function->first ) != beforeFunctions.end() )
         WALBERLA_ABORT( "There are at least two before functions with identifier \"" << function->first <<
                         "\" (test is case insensitive!).\nEvery function must have a unique identifier!" );
      beforeFunctions[ function->first ] = function->second;
   }

   Config::BlockHandle vtkBlock = parentBlockHandle.getBlock( configBlockName );

   if( !vtkBlock )
      return;

   Config::Blocks blocks;
   vtkBlock.getBlocks( blocks );

   for( auto block = blocks.begin(); block != blocks.end(); ++block )
   {
      Config::Blocks subBlocks;
      block->getBlocks( subBlocks );

      const std::string identifier( block->getKey() );

      const int    simultaneousIOOperations = block->getParameter< int >( "simultaneousIOOperations", 0 );
      const uint_t writeFrequency           = block->getParameter< uint_t >( "writeFrequency", uint_c(1) );
      const uint_t initialExecutionCount    = block->getParameter< uint_t >( "initialExecutionCount", uint_c(0) );
      const bool   forcePVTU                = block->getParameter< bool >( "forcePVTU", false );
      const uint_t ghostLayers              = block->getParameter< uint_t >( "ghostLayers", uint_c(0) );

      std::string baseFolder( block->getParameter< std::string >( "baseFolder", std::string( "vtk_out" ) ) );
      std::string executionFolder( block->getParameter< std::string >( "executionFolder", std::string( "simulation_step" ) ) );

      const bool outputDomainDecomposition = block->getParameter< bool >( "outputDomainDecomposition", false );
      const bool continuousNumbering       = block->getParameter< bool >( "continuousNumbering", false );
      const bool binary                    = block->getParameter< bool >( "binary", true );
      const bool littleEndian              = block->getParameter< bool >( "littleEndian", true );
      const bool useMPIIO                  = block->getParameter< bool >( "useMPIIO", true );

      Config::BlockHandle writersBlock = block->getBlock( "writers" );

      if( !writersBlock && !outputDomainDecomposition )
         WALBERLA_ABORT( "You declared a VTK output instance [\"" << identifier << "\"] without a \"writers\" block. "
                         "You have to specify at least on block data writer!" );

      if( useMPIIO && simultaneousIOOperations > 0 )
      {
         WALBERLA_LOG_WARNING( "In VTK output instance [\"" << identifier << "\"] you request the use of MPI I/O and "
                               "specified \"simultaneousIOOperations\". Those two settings are incompatible. "
                               "\"simultaneousIOOperations\" will be ignored!" );
      }

      std::vector< shared_ptr< BlockCellDataWriterInterface > > selectedWriters;

      if( writersBlock )
      {
         for( auto writerId = writersBlock.begin(); writerId != writersBlock.end(); ++writerId )
         {
            if( writers.find( writerId->first ) != writers.end() )
               selectedWriters.push_back( writers[ writerId->first ] );
            else
               WALBERLA_ABORT( "You have requested a block data writer \"" << writerId->first << "\". This writer is not available!" );
         }
      }

      if( selectedWriters.empty() && !outputDomainDecomposition )
         WALBERLA_ABORT( "No block data writers could be selected for your VTK output instance [\"" << identifier << "\"]. "
                         "Either you did not specify any writers or non of your specified writers are available." );

      shared_ptr< VTKOutput > vtkOutput;
      if( outputDomainDecomposition )
         vtkOutput = createVTKOutput_DomainDecomposition( *storage, identifier, writeFrequency, baseFolder, executionFolder,
                                                          continuousNumbering, binary, littleEndian, useMPIIO, initialExecutionCount );
      else
         vtkOutput = createVTKOutput_BlockData( *storage, identifier, writeFrequency, ghostLayers, forcePVTU,
                                                baseFolder, executionFolder, continuousNumbering, binary, littleEndian, useMPIIO, initialExecutionCount );

      const uint_t initialWriteCallsToSkip = block->getParameter< uint_t >( "initialWriteCallsToSkip", uint_t(0) );
      if( initialWriteCallsToSkip > uint_t(0) )
         vtkOutput->setInitialWriteCallsToSkip( initialWriteCallsToSkip );

      const real_t samplingResolution = block->getParameter< real_t >( "samplingResolution", real_c(-1) );
      vtkOutput->setSamplingResolution( samplingResolution );

      if( block->isDefined( "samplingDx" ) )
      {
         const real_t samplingDx = block->getParameter< real_t >( "samplingDx", real_c(-1) );
         const real_t samplingDy = block->getParameter< real_t >( "samplingDy", real_c(-1) );
         const real_t samplingDz = block->getParameter< real_t >( "samplingDz", real_c(-1) );

         vtkOutput->setSamplingResolution( samplingDx, samplingDy, samplingDz );
      }

      Config::BlockHandle beforeFunctionsBlock = block->getBlock( "before_functions" );
      if( beforeFunctionsBlock )
      {
         for( auto beforeFunctionId = beforeFunctionsBlock.begin(); beforeFunctionId != beforeFunctionsBlock.end(); ++beforeFunctionId )
         {
            if( beforeFunctions.find( beforeFunctionId->first ) != beforeFunctions.end() )
               vtkOutput->addBeforeFunction( beforeFunctions[ beforeFunctionId->first ] );
            else
               WALBERLA_ABORT( "You have requested a before function \"" << beforeFunctionId->first << "\". This function is not available!" );
         }
      }

      std::map< std::string, VTKOutput::CellFilter, CaseInsensitiveCompare > filters( globalFilters );

      Config::Blocks aabbBlocks;
      Config::Blocks cellBBBlocks;

      for( auto subBlock = subBlocks.begin(); subBlock != subBlocks.end(); ++subBlock )
      {
         if( string_icompare( std::string( subBlock->getKey(), 0, 11 ), std::string("AABB_filter") ) == 0 )
            aabbBlocks.push_back( *subBlock );

         if( string_icompare( std::string( subBlock->getKey(), 0, 13 ), std::string("CellBB_filter") ) == 0 )
            cellBBBlocks.push_back( *subBlock );
      }

      for( auto aabb = aabbBlocks.begin(); aabb != aabbBlocks.end(); ++aabb )
      {
         if( !aabb->isDefined("min") || !aabb->isDefined("max") )
            WALBERLA_ABORT( "You must specify a \"min\" and a \"max\" coordinate for AABB cell filter \"" << aabb->getKey() << "\"." );

         real_t xmin, ymin, zmin;
         splitVector< real_t >( xmin, ymin, zmin, *aabb, "min", std::string( "The \"min\" coordinate of AABB cell filter \"" ) + aabb->getKey() +
                                                                std::string( "\" must be a three-dimensional vector." ) );
         real_t xmax, ymax, zmax;
         splitVector< real_t >( xmax, ymax, zmax, *aabb, "max", std::string( "The \"max\" coordinate of AABB cell filter \"" ) + aabb->getKey() +
                                                                std::string( "\" must be a three-dimensional vector." ) );

         if( filters.find( aabb->getKey() ) != filters.end() )
            WALBERLA_ABORT( "There are at least two cell filters with identifier \"" << aabb->getKey() <<
                            "\" (test is case insensitive!).\nEvery filter must have a unique identifier!" );

         filters[ aabb->getKey() ] = AABBCellFilter( AABB( xmin, ymin, zmin, xmax, ymax, zmax ) );
      }

      for( auto bb = cellBBBlocks.begin(); bb != cellBBBlocks.end(); ++bb )
      {
         if( !bb->isDefined("min") || !bb->isDefined("max") )
            WALBERLA_ABORT( "You must specify a \"min\" and a \"max\" coordinate for CellBB cell filter \"" << bb->getKey() << "\"." );

         uint_t level = 0;
         if( bb->isDefined( "level" ) )
            level = bb->getParameter< uint_t >( "level" );

         cell_idx_t xmin, ymin, zmin;
         splitVector< cell_idx_t >( xmin, ymin, zmin, *bb, "min", std::string( "The \"min\" coordinate of CellBB cell filter \"" ) + bb->getKey() +
                                                                  std::string( "\" must be a three-dimensional vector." ) );
         cell_idx_t xmax, ymax, zmax;
         splitVector< cell_idx_t >( xmax, ymax, zmax, *bb, "max", std::string( "The \"max\" coordinate of CellBB cell filter \"" ) + bb->getKey() +
                                                                  std::string( "\" must be a three-dimensional vector." ) );

         if( filters.find( bb->getKey() ) != filters.end() )
            WALBERLA_ABORT( "There are at least two cell filters with identifier \"" << bb->getKey() <<
                            "\" (test is case insensitive!).\nEvery filter must have a unique identifier!" );

         filters[ bb->getKey() ] = CellBBCellFilter( CellInterval( xmin, ymin, zmin, xmax, ymax, zmax ), level );
      }

      for( auto selectedWriter = selectedWriters.begin(); selectedWriter != selectedWriters.end(); ++selectedWriter )
         vtkOutput->addCellDataWriter( *selectedWriter );

      Config::BlockHandle inclusionFiltersBlock = block->getBlock( "inclusion_filters" );
      if( inclusionFiltersBlock )
      {
         for( auto inclusionFilterId = inclusionFiltersBlock.begin(); inclusionFilterId != inclusionFiltersBlock.end(); ++inclusionFilterId )
         {
            if( inclusionFilterId->first == "combine" )
            {
               std::vector< std::string > filterList = splitList( inclusionFilterId->second );
               ChainedFilter combine;
               for( auto filter = filterList.begin(); filter != filterList.end(); ++filter )
               {
                  if( filters.find( *filter ) != filters.end() )
                     combine.addFilter( filters[ *filter ] );
                  else
                     WALBERLA_ABORT( "You have requested an inclusion cell filter \"" << *filter << "\". This filter is not available!" );
               }
               vtkOutput->addCellInclusionFilter( combine );
            }
            else
            {
               if( filters.find( inclusionFilterId->first ) != filters.end() )
                  vtkOutput->addCellInclusionFilter( filters[ inclusionFilterId->first ] );
               else
                  WALBERLA_ABORT( "You have requested an inclusion cell filter \"" << inclusionFilterId->first << "\". This filter is not available!" );
            }
         }
      }

      Config::BlockHandle exclusionFiltersBlock = block->getBlock( "exclusion_filters" );
      if( exclusionFiltersBlock )
      {
         for( auto exclusionFilterId = exclusionFiltersBlock.begin(); exclusionFilterId != exclusionFiltersBlock.end(); ++exclusionFilterId )
         {
            if( exclusionFilterId->first == "combine" )
            {
               std::vector< std::string > filterList = splitList( exclusionFilterId->second );
               ChainedFilter combine;
               for( auto filter = filterList.begin(); filter != filterList.end(); ++filter )
               {
                  if( filters.find( *filter ) != filters.end() )
                     combine.addFilter( filters[ *filter ] );
                  else
                     WALBERLA_ABORT( "You have requested an exclusion cell filter \"" << *filter << "\". This filter is not available!" );
               }
               vtkOutput->addCellExclusionFilter( combine );
            }
            else
            {
               if( filters.find( exclusionFilterId->first ) != filters.end() )
                  vtkOutput->addCellExclusionFilter( filters[ exclusionFilterId->first ] );
               else
                  WALBERLA_ABORT( "You have requested an exclusion cell filter \"" << exclusionFilterId->first << "\". This filter is not available!" );
            }
         }
      }

      Set<SUID> requiredGlobalStates;
      if( block->isDefined( "requiredGlobalStates" ) )
      {
         std::string states = block->getParameter< std::string >( "requiredGlobalStates" );
         addStates( requiredGlobalStates, states );
      }

      Set<SUID> incompatibleGlobalStates;
      if( block->isDefined( "incompatibleGlobalStates" ) )
      {
         std::string states = block->getParameter< std::string >( "incompatibleGlobalStates" );
         addStates( incompatibleGlobalStates, states );
      }

      Set<SUID> requiredBlockStates;
      if( block->isDefined( "requiredBlockStates" ) )
      {
         std::string states = block->getParameter< std::string >( "requiredBlockStates" );
         addStates( requiredBlockStates, states );
      }

      Set<SUID> incompatibleBlockStates;
      if( block->isDefined( "incompatibleBlockStates" ) )
      {
         std::string states = block->getParameter< std::string >( "incompatibleBlockStates" );
         addStates( incompatibleBlockStates, states );
      }

      outputFunctions[ identifier ] = SelectableOutputFunction( writeFiles( vtkOutput, true, simultaneousIOOperations,
                                                                                requiredGlobalStates +     requiredBlockStates,
                                                                            incompatibleGlobalStates + incompatibleBlockStates ),
                                                                requiredGlobalStates, incompatibleGlobalStates );
   }
}



void initializeVTKOutput( std::map< std::string, SelectableOutputFunction > & outputFunctions,
                          const shared_ptr< const StructuredBlockStorage > & storage, const shared_ptr< Config > & config,
                          const std::vector< shared_ptr< BlockCellDataWriterInterface > > & writers,
                          const std::map< std::string, VTKOutput::CellFilter > & filters,
                          const std::map< std::string, VTKOutput::BeforeFunction > & beforeFunctions )
{
   if( !!config )
      initializeVTKOutput( outputFunctions, storage, config->getGlobalBlock(), "VTK", writers, filters, beforeFunctions );
}



void initializeVTKOutput( std::map< std::string, SelectableOutputFunction > & outputFunctions,
                          const shared_ptr< const StructuredBlockStorage > & storage,
                          const Config::BlockHandle & parentBlockHandle,
                          const std::vector< shared_ptr< BlockCellDataWriterInterface > > & writers,
                          const std::map< std::string, VTKOutput::CellFilter > & filters,
                          const std::map< std::string, VTKOutput::BeforeFunction > & beforeFunctions )
{
   initializeVTKOutput( outputFunctions, storage, parentBlockHandle, "VTK", writers, filters, beforeFunctions );
}



//**********************************************************************************************************************
/*!
*   \brief Function for initializing VTKOutput objects from file and creating their corresponding output functions
*
*   This initialization function reads data stored in a configuration file, uses this data to create VTKOutput objects,
*   and finally constructs corresponding output functions which then are returned via the parameter "outputFunctions".
*   The returned output functions all have the signature "void (void)". Calling such an output function initiates one
*   output of the associated VTKOutput object.
*
*   \section docVTKConfigurationFile VTK via Configuration File
*
*   For reading VTK setup information from a configuration file, the structure of the configuration file must look
*   like as follows:
*
*   \code
*   VTK // every sub block of VTK corresponds to one VTKOutput object that is created
*   {
*      [name] // identifier for this VTKOutput object
*      {
*         simultaneousIOOperations [integer value]; // max. number of files that are written
*                                                   // in parallel (optional, default=0 [=disabled])
*
*         initialExecutionCount   [integer value]; // optional, default=0
*         initialWriteCallsToSkip [integer value]; // optional, default=0
*         writeFrequency          [integer value]; // output frequency = number of times "writeBlocks" must be
*                                                  // called for triggering an output (optional, default=1)
*                                                  // 0 disables output
*         forcePVTU               [boolean];       // if true, (P)VTU files are created, and not (P)VTI files
*                                                  // (optional, default=false)
*         ghostLayers             [integer value]; // number of ghost layers (optional, default=0)
*
*         baseFolder      [directory]; // base directory (optional, default=vtk_out)
*         executionFolder [directory]; // base directory for each time step, directory path is given relative
*                                      // to baseFolder (optional, default=simulation_step)
*
*         outputDomainDecomposition [boolean]; // if true, the domain decomposition is written to file
*                                              // specifying cell filters and block cell data writers is
*                                              // not allowed! (optional, default=false)
*         continuousNumbering       [boolean]; // if false, the actual time step is preserved
*                                              // (optional, default=false)
*         binary                    [boolean]; // if false, ascii files are written (optional, default=true)
*         littleEndian              [boolean]; // switch between little and big endianness (optional, default=true)
*
*         useMPIIO                  [boolean]; // use MPI I/O to write only one file per time step
*                                              // (optional, default=true)
*
*         // You can either specify "samplingResolution" or "samplingDx", "samplingDy", and "samplingDz"
*         samplingResolution [floating point value]; // "samplingResolution VALUE" has the same effect as
*                                                    // setting "samplingDx", "samplingDy", and "samplingDz" to VALUE
*         samplingDx         [floating point value]; // forces the output to use this dx (= cell x-spacing)
*         samplingDy         [floating point value]; // forces the output to use this dy (= cell y-spacing)
*         samplingDz         [floating point value]; // forces the output to use this dz (= cell z-spacing)
*
*         before_functions // (OPTIONAL, APPLICATION-DEPENDENT!)
*         {
*            [NameOfTheFirstFunction];   // the mapping of this name to a function pointer/functor is
*            [NameOfTheSecondFunction];  // done by the RegisterVTKOutputFunction "registerVTKOutputFunction"
*            [...]                       // which must be implemented by the user
*         }
*
*         // AABB filters are OPTIONAL. AABB filters are sub blocks that must start with "AABB_filter".
*         // AABB filters can be selected as either inclusion or exclusion filters. In order to select an
*         // AABB filter "AABB_filter_XXX" as inclusion/exclusion filter, its name/identifier (in this
*         // example: "AABB_filter_XXX") must be added to the list of inclusion/exclusion filters in
*         // the sub block inclusion_filters/inclusion_filters.
*
*         AABB_filter_0
*         {
*            min < [x: floating point value], [y: floating point value], [z: floating point value] >;
*            max < [x: floating point value], [y: floating point value], [z: floating point value] >;
*         }
*         // AABB_filter_1
*         // {
*         //    min < [x: floating point value], [y: floating point value], [z: floating point value] >;
*         //    max < [x: floating point value], [y: floating point value], [z: floating point value] >;
*         // }
*         // AABB_filter_2 { ... }
*         // AABB_filter_* { ... }
*
*         // CellBB filters are OPTIONAL. CellBB filters are sub blocks that must start with "CellBB_filter".
*         // CellBB filters can be selected as either inclusion or exclusion filters. In order to select an
*         // CellBB filter "CellBB_filter_XXX" as inclusion/exclusion filter, its name/identifier (in this
*         // example: "CellBB_filter_XXX") must be added to the list of inclusion/exclusion filters in
*         // the sub block inclusion_filters/inclusion_filters.
*         // CellBB filters are AABB filters that are defined using discrete (!) global cell coordinates.
*
*         CellBB_filter_0
*         {
*            level [integer value]; // the cell level to which the following coordinates correspond to
*                                   // (optional, default=0)
*            min < [x: integer value], [y: integer value], [z: integer value] >;
*            max < [x: integer value], [y: integer value], [z: integer value] >;
*         }
*         // CellBB_filter_1
*         // {
*         //    min < [x: integer value], [y: integer value], [z: integer value] >;
*         //    max < [x: integer value], [y: integer value], [z: integer value] >;
*         // }
*         // CellBB_filter_2 { ... }
*         // CellBB_filter_* { ... }
*
*         // In terms of set theory: all filters listed as inclusion filters are "added" together,
*         // resulting in a union of all filters.
*
*         inclusion_filters // (OPTIONAL, APPLICATION-DEPENDENT!)
*         {
*            [NameOfTheFirstFilter];   // the mapping of this name to an inclusion filter is
*            [NameOfTheSecondFilter];  // done by the RegisterVTKOutputFunction "registerVTKOutputFunction"
*            [...]                     // which must be implemented by the user
*
*            // In terms of set theory: combining filters results in an intersection of these filters.
*            combine [NameOfAFilter],[NameOfAnotherFilter],[NameOfYetAnotherFilter],[...];
*         }
*
*         // In terms of set theory: all filters listed as exclusion filters are "added" together,
*         // resulting in a union of all filters.
*
*         exclusion_filters // (OPTIONAL, APPLICATION-DEPENDENT!)
*         {
*            [NameOfTheFirstFilter];   // the mapping of this name to an exclusion filter is
*            [NameOfTheSecondFilter];  // done by the RegisterVTKOutputFunction "registerVTKOutputFunction"
*            [...]                     // which must be implemented by the user
*
*            // In terms of set theory: combining filters results in an intersection of these filters.
*            combine [NameOfAFilter],[NameOfAnotherFilter],[NameOfYetAnotherFilter],[...];
*         }
*
*         writers // (AT LEAST ONE IS MANDATORY [exception: if outputDomainDecomposition == true,
*                 // no writers are allowed!], APPLICATION-DEPENDENT!)
*         {
*            [NameOfTheFirstWriter];  // the mapping of this name to a writer is
*            [NameOfTheSecondWriter]; // done by the RegisterVTKOutputFunction "registerVTKOutputFunction"
*            [...]                    // which must be implemented by the user
*         }
*
*         requiredGlobalStates      [SUID identifier #1], [SUID identifier #2], ...; // (optional, default=[none])
*         incompatibleGlobalStates  [SUID identifier #1], [SUID identifier #2], ...; // (optional, default=[none])
*         requiredBlockStates       [SUID identifier #1], [SUID identifier #2], ...; // (optional, default=[none])
*         incompatibleBlockStates   [SUID identifier #1], [SUID identifier #2], ...; // (optional, default=[none])
*      }
*
*      [name of another VTK output object]
*      {
*         [...]
*      }
*
*   } // VTK
*   \endcode
*
*   \param outputFunctions           The output functions which correspond to the just created VTKOutput objects
*   \param registerVTKOutputFunction A function pointer that the user must provide and that is used for
*                                    registering cell filters and block data writers which then can be referenced in the
*                                    configuration file and which are used to assemble the VTKOutput objects
*   \param storage                   The structured block storage the VTKOutput object shall be associated with
*   \param config                    The configuration
*   \param configBlockName           Name of the block in the configuration that is used to setup the VTK output
*/
//**********************************************************************************************************************
void initializeVTKOutput( std::map< std::string, SelectableOutputFunction > & outputFunctions, const RegisterVTKOutputFunction& registerVTKOutputFunction,
                          const shared_ptr< const StructuredBlockStorage > & storage, const shared_ptr< Config > & config,
                          const std::string & configBlockName )
{
   if( !!config )
      initializeVTKOutput( outputFunctions, registerVTKOutputFunction, storage, config->getGlobalBlock(), configBlockName );
}



void initializeVTKOutput( std::map< std::string, SelectableOutputFunction > & outputFunctions, const RegisterVTKOutputFunction& registerVTKOutputFunction,
                          const shared_ptr< const StructuredBlockStorage > & storage, const Config::BlockHandle & parentBlockHandle,
                          const std::string & configBlockName )
{
   std::vector< shared_ptr< BlockCellDataWriterInterface > > _writers;
   std::map< std::string, VTKOutput::CellFilter >            _filters;
   std::map< std::string, VTKOutput::BeforeFunction >        _beforeFunctions;

   registerVTKOutputFunction( _writers, _filters, _beforeFunctions );

   initializeVTKOutput( outputFunctions, storage, parentBlockHandle, configBlockName, _writers, _filters, _beforeFunctions );
}



} // namespace vtk
} // namespace walberla
