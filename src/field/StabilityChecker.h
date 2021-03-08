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
//! \file StabilityCheckerer.h
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/Set.h"
#include "core/cell/CellSet.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "core/math/Vector3.h"
#include "core/mpi/Reduce.h"
#include "core/uid/SUID.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/EvaluationFilter.h"
#include "field/iterators/IteratorMacros.h"

#include "vtk/BlockCellDataWriter.h"
#include "vtk/DumpBlockStructureLevel.h"
#include "vtk/DumpBlockStructureProcess.h"
#include "vtk/VTKOutput.h"

#include <type_traits>

namespace walberla {
namespace field {

namespace internal {

const std::string stabilityCheckerVTKBase("vtk_out");
const std::string stabilityCheckerVTKFolder("output");
const std::string stabilityCheckerVTKIdentifier("error_field");

const bool stabilityCheckerVTKBinary( true );
const bool stabilityCheckerVTKLittleEndian( true );
const bool stabilityCheckerVTKMPIIO( true );
const bool stabilityCheckerVTKForcePVTU( false );

const std::string stabilityCheckerConfigBlock("StabilityChecker");

template< typename T >
inline bool stabilityCheckerIsFinite( const T & value ) { return math::finite( value ); }

template<>
inline bool stabilityCheckerIsFinite( const Vector3<real_t> & value ) { return math::finite(value[0]) && math::finite(value[1]) && math::finite(value[2]); }

}



//**********************************************************************************************************************
/*!
*   \brief Class/functor for checking a running simulation for non-finite values in a specific field
*
*   \section docStabilityChecker Stability Checker
*
*   If non-finite values are detected in the field that is checked, the simulation is aborted. Optionally, information
*   about all cells that contain non-finite vales can be logged via the Logging or saved as VTK output for further
*   investigation.
*
*   Do not create objects of class StabilityChecker directly, better use one of the various 'makeStabilityChecker'
*   functions below!
*
*   Template parameters:
*   - Field_T: the field storing the simulation values (also works if the field stores data of type Vector3)
*   - Filter_T: the type of the evaluation filter (see \ref docEvaluationFilter in 'EvaluationFilter.h')
*
*   For the parameters for setting up and controlling the stability checker see the documentation of the constructor
*   of this class.
*
*   You do not have to specify an evaluation filter! If you do not specify any filter, _all_ cells are processed and no
*   cell is excluded.
*
*   If you want to use a flag field as evaluation filter, fitting 'makeStabilityChecker' functions already exist.
*   These functions need an additional template parameter FlagField_T and you have to provide the block data ID of the
*   flag field together with a set of flag UIDs that specify which cells need to be processed.
*
*   There also exist 'makeStabilityChecker' functions that take configuration file data as an additional parameter in
*   order to parse the configuration file for setting up and controlling the stability checker. The configuration file
*   block looks like as follows:
*
*   \code
*   StabilityChecker
*   {
*      checkFrequency     [unsigned integer]; // check frequency [default:0]
*      streamOutput       [boolean]; // output to stream? [default: true]
*      vtkOutput          [boolean]; // output to VTK? [default:true]
*      vtkBaseFolder      [string]; // VTK base folder [default: vtk_out]
*      vtkExecutionFolder [string]; // VTK execution folder [default:output]
*      vtkIdentifier      [string]; // VTK identifier [default: error_field]
*      vtkBinary          [boolean]; // write VTK data in binary? [default: true]
*      vtkLittleEndian    [boolean]; // VTK binary file format [default: true (= little endian)]
*      vtkMPIIO           [boolean]; // use MPI IO for creating VTK output? [default: true]
*      vtkForcePVTU       [boolean]; // force VTK to generate a PVTU file? [default: false]
*   }
*   \endcode
*
*   Example:
*
*   \code
*   StabilityChecker
*   {
*      checkFrequency 100;
*      streamOutput   false;
*      vtkOutput      true;
*      vtkBaseFolder  /home/anonymous/vtk;
*   }
*   \endcode
*
*   For documentation of the VTK parameters see \ref docVTKConfigurationFile.
*
*   Note that the shared pointer returned by all 'makeStabilityChecker' functions can be captured by a SharedFunctor
*   for immediate registration at a time loop (see field::makeSharedFunctor).
*/

template< typename Field_T, typename Filter_T = DefaultEvaluationFilter >
class StabilityChecker
{
private:

   using BlockCellsMap = std::map<const IBlock *, std::map<Cell, std::set<cell_idx_t>>>;

   ////////////////////////////////
   //  VTK Output Helper Classes //
   ////////////////////////////////

   /// This cell filter selects only those cells in which at least one non-finite value (= infinite or NaN) was detected
   class VTKCellFilter
   {
   public:

      VTKCellFilter( BlockCellsMap & map ) : map_( map ) {}

      void operator()( CellSet & filteredCells, const IBlock & block, const StructuredBlockStorage & storage, const uint_t ghostLayers ) const
      {
         if( map_.find( &block ) != map_.end() )
         {
            const auto & cellMap = map_[ &block ];
            if( !cellMap.empty() )
            {
               const cell_idx_t gl    = cell_idx_c( ghostLayers );
               const cell_idx_t begin = cell_idx_c( -1 ) * gl;

               for( cell_idx_t z = begin; z < cell_idx_c( storage.getNumberOfZCells(block) ) + gl; ++z )
                  for( cell_idx_t y = begin; y < cell_idx_c( storage.getNumberOfYCells(block) ) + gl; ++y )
                     for( cell_idx_t x = begin; x < cell_idx_c( storage.getNumberOfXCells(block) ) + gl; ++x )
                        if( cellMap.find( Cell(x,y,z) ) != cellMap.end() ) filteredCells.insert( x, y, z );
            }
         }
      }

   private:

      BlockCellsMap & map_;
   };

   /// For each value of a cell, either '0' or '1' is stored in the VTK file.
   /// If the value is non-finite (= infinite or NaN), '1' is written to file - otherwise '0'.
   class FValueVTKWriter : public vtk::BlockCellDataWriter< uint8_t, Field_T::F_SIZE >
   {
   public:

      FValueVTKWriter( BlockCellsMap & map, const std::string & id ) : vtk::BlockCellDataWriter< uint8_t, Field_T::F_SIZE >( id ), map_( map ) {}

   protected:

      void configure() override {}

      uint8_t evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f ) override
      {
         WALBERLA_ASSERT( map_.find( this->block_ ) != map_.end() );
         WALBERLA_ASSERT( map_[ this->block_ ].find( Cell(x,y,z) ) != map_[ this->block_ ].end() );

         return ( map_[ this->block_ ][ Cell(x,y,z) ].find( f ) != map_[ this->block_ ][ Cell(x,y,z) ].end() ) ? uint8_t(1) : uint8_t(0);
      }

   private:

      BlockCellsMap & map_;
   };

   /// For each cell, the corresponding block local cell coordinates are stored in the VTK file.
   class LocalCoordVTKWriter : public vtk::BlockCellDataWriter< cell_idx_t, 3 >
   {
   public:
      LocalCoordVTKWriter( const std::string & id ) : vtk::BlockCellDataWriter< cell_idx_t, 3 >( id ) {}
   protected:
      void configure() override {}
      cell_idx_t evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f ) override
      {
         return ( f == cell_idx_t(0) ) ? x : ( ( f == cell_idx_t(1) ) ? y : z );
      }
   };

   /// For each cell, the corresponding global cell coordinates are stored in the VTK file.
   class GlobalCoordVTKWriter : public vtk::BlockCellDataWriter< cell_idx_t, 3 >
   {
   public:
      GlobalCoordVTKWriter( const std::string & id ) : vtk::BlockCellDataWriter< cell_idx_t, 3 >( id ) {}
   protected:
      void configure() override {}
      cell_idx_t evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f ) override
      {
         Cell cell(x,y,z);
         this->blockStorage_->transformBlockLocalToGlobalCell( cell, *(this->block_) );
         return ( f == cell_idx_t(0) ) ? cell.x(): ( ( f == cell_idx_t(1) ) ? cell.y() : cell.z() );
      }
   };

public:

   //*******************************************************************************************************************
   /*!
   *   \brief Constructor for class 'StabilityChecker'
   *
   *   \param blocks                Shared pointer to a structured block storage
   *   \param fieldId               Block data ID of the field that will be checked
   *   \param filter                The evaluation filter that indicates which cells are processed
   *   \param checkFrequency        If operator()() is called, the stability check is only performed every
   *                                'checkFrequency'-th time. Setting 'checkFrequency' to 1 means the stability check
   *                                is performed each time operator()() is called. Setting 'checkFrequency' to 0
   *                                disables the check entirely.
   *   \param outputToStream        If true, in case a non-finite value is detected in the field, information about the
   *                                corresponding cells is logged via WALBERLA_LOG_WARNING.
   *   \param outputVTK             If true, in case a non-finite value is detected in the field, VTK output is
   *                                generated and information about the corresponding cells is saved.
   *   \param requiredSelectors     Required selectors
   *   \param incompatibleSelectors Incompatible selectors
   */
   //*******************************************************************************************************************
   StabilityChecker( const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                     const Filter_T & filter, const uint_t checkFrequency,
                     const bool outputToStream = true, const bool outputVTK = true,
                     const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                     const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      blocks_( blocks ), filter_( filter ), executionCounter_( uint_c(0) ), checkFrequency_( checkFrequency ),
      fieldId_( fieldId ), outputToStream_( outputToStream ), outputVTK_( outputVTK ),
      vtkBaseFolder_( internal::stabilityCheckerVTKBase ),
      vtkExecutionFolder_( internal::stabilityCheckerVTKFolder ),
      vtkIdentifier_( internal::stabilityCheckerVTKIdentifier ),
      vtkBinary_( internal::stabilityCheckerVTKBinary ),
      vtkLittleEndian_( internal::stabilityCheckerVTKLittleEndian ),
      vtkMPIIO_( internal::stabilityCheckerVTKMPIIO ),
      vtkForcePVTU_( internal::stabilityCheckerVTKForcePVTU ),
      requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors ) {}

   StabilityChecker( const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                     const uint_t checkFrequency,
                     const bool outputToStream = true, const bool outputVTK = true,
                     const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                     const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      blocks_( blocks ), filter_( Filter_T() ), executionCounter_( uint_c(0) ), checkFrequency_( checkFrequency ),
      fieldId_( fieldId ), outputToStream_( outputToStream ), outputVTK_( outputVTK ),
      vtkBaseFolder_( internal::stabilityCheckerVTKBase ),
      vtkExecutionFolder_( internal::stabilityCheckerVTKFolder ),
      vtkIdentifier_( internal::stabilityCheckerVTKIdentifier ),
      vtkBinary_( internal::stabilityCheckerVTKBinary ),
      vtkLittleEndian_( internal::stabilityCheckerVTKLittleEndian ),
      vtkMPIIO_( internal::stabilityCheckerVTKMPIIO ),
      vtkForcePVTU_( internal::stabilityCheckerVTKForcePVTU ),
      requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors )
   {
      static_assert( (std::is_same< Filter_T, DefaultEvaluationFilter >::value),
                     "This constructor is only available if DefaultEvaluationFilter is set as filter type!" );
   }
   
   void setVTKBaseFolder     ( const std::string & vtkBaseFolder      ) { vtkBaseFolder_      = vtkBaseFolder; }
   void setVTKExecutionFolder( const std::string & vtkExecutionFolder ) { vtkExecutionFolder_ = vtkExecutionFolder; }
   void setVTKIdentifier     ( const std::string & vtkIdentifier      ) { vtkIdentifier_      = vtkIdentifier; }
   
   void setVTKBinary      ( const bool vtkBinary       ) { vtkBinary_       = vtkBinary; }
   void setVTKLittleEndian( const bool vtkLittleEndian ) { vtkLittleEndian_ = vtkLittleEndian; }
   void setVTKMPIIO       ( const bool vtkMPIIO        ) { vtkMPIIO_        = vtkMPIIO; }
   void setVTKForcePVTU   ( const bool vtkForcePVTU    ) { vtkForcePVTU_    = vtkForcePVTU; }
   
   void operator()();

private:

   void checkBlock( const IBlock * const block );



   weak_ptr< StructuredBlockStorage > blocks_;
   
   Filter_T filter_;

   uint_t executionCounter_;
   uint_t checkFrequency_;

   ConstBlockDataID fieldId_;

   BlockCellsMap failedCells_;

   bool outputToStream_;
   bool outputVTK_;
   
   std::string vtkBaseFolder_;
   std::string vtkExecutionFolder_;
   std::string vtkIdentifier_;

   bool vtkBinary_;
   bool vtkLittleEndian_;
   bool vtkMPIIO_;
   bool vtkForcePVTU_;

   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;

}; // StabilityChecker



template< typename Field_T, typename Filter_T >
void StabilityChecker< Field_T, Filter_T >::operator()()
{
   ++executionCounter_;
   if( checkFrequency_ == uint_t(0) || ( executionCounter_ - uint_c(1) ) % checkFrequency_ != 0 )
      return;

   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'StabilityChecker' for a block storage object that doesn't exist anymore" );

   for( auto block = blocks->begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks->end(); ++block )
      checkBlock( block.get() );

   if( outputToStream_ )
   {
      std::ostringstream oss;
      oss << "Check for finiteness failed on process " << MPIManager::instance()->rank();

      for( auto block = failedCells_.begin(); block != failedCells_.end(); ++block )
      {
         oss << "\n - on block " << block->first->getId() <<
                "\n - on level " << blocks->getLevel( *(block->first) )  <<
                "\n - for:";

         for( auto cell = block->second.begin(); cell != block->second.end(); ++cell )
         {
            const cell_idx_t x = cell->first[0];
            const cell_idx_t y = cell->first[1];
            const cell_idx_t z = cell->first[2];

            Vector3< real_t > center;
            blocks->getBlockLocalCellCenter( *(block->first), cell->first, center );

            Cell gCell(x,y,z);
            blocks->transformBlockLocalToGlobalCell( gCell, *(block->first) );

            for( auto f = cell->second.begin(); f != cell->second.end(); ++f )
            {
               oss << "\n   + block local cell( " << x << ", " << y << ", " << z << " ) at index " << *f <<
                      "\n     = global cell ( " << gCell.x() << ", " << gCell.y() << ", " << gCell.z() << " ) with "
                      "cell center ( " << center[0] << ", " << center[1] << ", " << center[2] << " )";
            }
         }
      }

      if( !failedCells_.empty() )
         WALBERLA_LOG_WARNING( oss.str() );
   }

   bool abort = !failedCells_.empty();
   mpi::allReduceInplace( abort, mpi::LOGICAL_OR );

   if( abort )
   {
      if( outputVTK_ )
      {
         auto vtkWriter = vtk::createVTKOutput_BlockData( blocks, vtkIdentifier_, uint_t(1), uint_t(0), vtkForcePVTU_,
                                                          vtkBaseFolder_, vtkExecutionFolder_, false, vtkBinary_, vtkLittleEndian_, vtkMPIIO_ );                                                         

         vtkWriter->addCellInclusionFilter( VTKCellFilter( failedCells_ ) );

         vtkWriter->addCellDataWriter( walberla::make_shared< vtk::DumpBlockStructureProcess >( "process" ) );
         vtkWriter->addCellDataWriter( walberla::make_shared< vtk::DumpBlockStructureLevel >( "level" ) );
         vtkWriter->addCellDataWriter( walberla::make_shared< FValueVTKWriter >( std::ref( failedCells_ ), "F" ) );
         vtkWriter->addCellDataWriter( walberla::make_shared< LocalCoordVTKWriter >( "blockLocalCoordinate" ) );
         vtkWriter->addCellDataWriter( walberla::make_shared< GlobalCoordVTKWriter >( "globalCoordinate" ) );

         vtkWriter->write();
      }

      WALBERLA_LOG_WARNING_ON_ROOT( "Field stability check failed - aborting program ..." );
      WALBERLA_MPI_WORLD_BARRIER();

      WALBERLA_ABORT_NO_DEBUG_INFO("");
   }
}



template< typename Field_T, typename Filter_T >
void StabilityChecker< Field_T, Filter_T >::checkBlock( const IBlock * const block )
{
   const Field_T * field = block->getData< Field_T >(  fieldId_ );
   
   filter_( *block );
   
#ifndef _OPENMP

      WALBERLA_FOR_ALL_CELLS_XYZ( field,

         if( filter_(x,y,z) )
         {
            for( uint_t f = uint_t(0); f < Field_T::F_SIZE; ++f )
            {
               if( !internal::stabilityCheckerIsFinite( field->get( x, y, z, cell_idx_c(f) ) ) )
                  failedCells_[ block ][ Cell(x,y,z) ].insert( cell_idx_c(f) );
            }
         }
      )

#else

   // WALBERLA_FOR_ALL_CELLS macros cannot be used since they do not support additional omp pragmas inside the kernel.
   // The additional omp critical section, however, is required.

   const CellInterval & size = field->xyzSize();

   if( size.zSize() >= size.ySize() )
   {
      const int izSize = int_c( size.zSize() );
      #pragma omp parallel for schedule(static)
      for( int iz = 0; iz < izSize; ++iz )
      {
         cell_idx_t z = cell_idx_c( iz );
         for( cell_idx_t y = size.yMin(); y <= size.yMax(); ++y ) {
            for( cell_idx_t x = size.xMin(); x <= size.xMax(); ++x )
            {
               if( filter_(x,y,z) )
               {
                  for( uint_t f = uint_t(0); f < Field_T::F_SIZE; ++f )
                  {
                     if( !internal::stabilityCheckerIsFinite( field->get( x, y, z, cell_idx_c(f) ) ) )
                     {
                        #pragma omp critical (StabilityChecker)
                        {
                           failedCells_[ block ][ Cell(x,y,z) ].insert( cell_idx_c(f) );
                        }
                     }
                  }
               }
            }
         }
      }
   }
   else
   {
      const int iySize = int_c( size.ySize() );
      #pragma omp parallel for schedule(static)
      for( int iy = 0; iy < iySize; ++iy )
      {
         cell_idx_t y = cell_idx_c( iy );
         for( cell_idx_t z = size.zMin(); z <= size.zMax(); ++z ) {
            for( cell_idx_t x = size.xMin(); x <= size.xMax(); ++x )
            {
               if( filter_(x,y,z) )
               {
                  for( uint_t f = uint_t(0); f < Field_T::F_SIZE; ++f )
                  {
                     if( !internal::stabilityCheckerIsFinite( field->get( x, y, z, cell_idx_c(f) ) ) )
                     {
                        #pragma omp critical (StabilityChecker)
                        {
                           failedCells_[ block ][ Cell(x,y,z) ].insert( cell_idx_c(f) );
                        }
                     }
                  }
               }
            }
         }
      }
   }

#endif

}



///////////////////////////////////////////////////////////////
// makeStabilityChecker functions without configuration file //
///////////////////////////////////////////////////////////////

template< typename Field_T >
shared_ptr< StabilityChecker< Field_T > > makeStabilityChecker( const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                                                                const uint_t checkFrequency,
                                                                const bool outputToStream = true, const bool outputVTK = true,
                                                                const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                                                const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using SC_T = StabilityChecker<Field_T>;
   return shared_ptr< SC_T >( new SC_T( blocks, fieldId, checkFrequency, outputToStream, outputVTK, requiredSelectors, incompatibleSelectors ) );
}

template< typename Field_T, typename FlagField_T >
shared_ptr< StabilityChecker< Field_T, FlagFieldEvaluationFilter<FlagField_T> > >
makeStabilityChecker( const weak_ptr< StructuredBlockStorage > & blocks,
                      const ConstBlockDataID & fieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                      const uint_t checkFrequency,
                      const bool outputToStream = true, const bool outputVTK = true,
                      const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                      const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using SC_T = StabilityChecker<Field_T, FlagFieldEvaluationFilter<FlagField_T>>;
   return shared_ptr< SC_T >( new SC_T( blocks, fieldId, FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                        checkFrequency, outputToStream, outputVTK, requiredSelectors, incompatibleSelectors ) );
}

template< typename Field_T, typename Filter_T >
shared_ptr< StabilityChecker< Field_T, Filter_T > >
makeStabilityChecker( const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                      const Filter_T & filter, const uint_t checkFrequency,
                      const bool outputToStream = true, const bool outputVTK = true,
                      const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                      const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   using SC_T = StabilityChecker<Field_T, Filter_T>;
   return shared_ptr< SC_T >( new SC_T( blocks, fieldId, filter, checkFrequency, outputToStream, outputVTK, requiredSelectors, incompatibleSelectors ) );
}



///////////////////////////////////////////////////////////
// makeStabilityChecker functions + configuration file //
///////////////////////////////////////////////////////////

namespace internal {

inline void stabilityCheckerConfigParser( const Config::BlockHandle & parentBlockHandle, const std::string & configBlockName,
                                          uint_t & defaultCheckFrequency, bool & defaultOutputToStream, bool & defaultOutputVTK,
                                          std::string & defaultVTKBaseFolder, std::string & defaultVTKExecutionFolder, std::string & defaultVTKIdentifier,
                                          bool & defaultVTKBinary, bool & defaultVTKLittleEndian, bool & defaultVTKMPIIO, bool & defaultVTKForcePVTU )
{
   if( parentBlockHandle )
   {
      Config::BlockHandle block = parentBlockHandle.getBlock( configBlockName );
      if( block )
      {
         defaultCheckFrequency = block.getParameter< uint_t >( "checkFrequency", defaultCheckFrequency );
         defaultOutputToStream = block.getParameter< bool >( "streamOutput", defaultOutputToStream );
         defaultOutputVTK = block.getParameter< bool >( "vtkOutput", defaultOutputVTK );
         defaultVTKBaseFolder = block.getParameter< std::string >( "vtkBaseFolder", defaultVTKBaseFolder );
         defaultVTKExecutionFolder = block.getParameter< std::string >( "vtkExecutionFolder", defaultVTKExecutionFolder );
         defaultVTKIdentifier = block.getParameter< std::string >( "vtkIdentifier", defaultVTKIdentifier );
         defaultVTKBinary = block.getParameter< bool >( "vtkBinary", defaultVTKBinary );
         defaultVTKLittleEndian = block.getParameter< bool >( "vtkLittleEndian", defaultVTKLittleEndian );
         defaultVTKMPIIO = block.getParameter< bool >( "vtkMPIIO", defaultVTKMPIIO );
         defaultVTKForcePVTU = block.getParameter< bool >( "vtkForcePVTU", defaultVTKForcePVTU );
      }
   }
}

inline void stabilityCheckerConfigParser( const shared_ptr< Config > & config, const std::string & configBlockName,
                                          uint_t & defaultCheckFrequency, bool & defaultOutputToStream, bool & defaultOutputVTK,
                                          std::string & defaultVTKBaseFolder, std::string & defaultVTKExecutionFolder, std::string & defaultVTKIdentifier,
                                          bool & defaultVTKBinary, bool & defaultVTKLittleEndian, bool & defaultVTKMPIIO, bool & defaultVTKForcePVTU )
{
   if( !!config )
      stabilityCheckerConfigParser( config->getGlobalBlock(), configBlockName, defaultCheckFrequency, defaultOutputToStream, defaultOutputVTK,
                                    defaultVTKBaseFolder, defaultVTKExecutionFolder, defaultVTKIdentifier,
                                    defaultVTKBinary, defaultVTKLittleEndian, defaultVTKMPIIO, defaultVTKForcePVTU );
}

} // namespace internal

#define WALBERLA_FIELD_MAKE_STABILITY_CHECKER_CONFIG_PARSER( config ) \
   uint_t defaultCheckFrequency = uint_t(0); \
   bool defaultOutputToStream = true; \
   bool defaultOutputVTK = true; \
   std::string defaultVTKBaseFolder = internal::stabilityCheckerVTKBase; \
   std::string defaultVTKExecutionFolder = internal::stabilityCheckerVTKFolder; \
   std::string defaultVTKIdentifier = internal::stabilityCheckerVTKIdentifier; \
   bool defaultVTKBinary = internal::stabilityCheckerVTKBinary; \
   bool defaultVTKLittleEndian = internal::stabilityCheckerVTKLittleEndian; \
   bool defaultVTKMPIIO = internal::stabilityCheckerVTKMPIIO; \
   bool defaultVTKForcePVTU = internal::stabilityCheckerVTKForcePVTU; \
   internal::stabilityCheckerConfigParser( config, configBlockName, defaultCheckFrequency, defaultOutputToStream, defaultOutputVTK, \
                                           defaultVTKBaseFolder, defaultVTKExecutionFolder, defaultVTKIdentifier, \
                                           defaultVTKBinary, defaultVTKLittleEndian, defaultVTKMPIIO, defaultVTKForcePVTU );

#define WALBERLA_FIELD_MAKE_STABILITY_CHECKER_SET_AND_RETURN() \
   checker->setVTKBaseFolder( defaultVTKBaseFolder ); \
   checker->setVTKExecutionFolder( defaultVTKExecutionFolder ); \
   checker->setVTKIdentifier( defaultVTKIdentifier ); \
   checker->setVTKBinary( defaultVTKBinary ); \
   checker->setVTKLittleEndian( defaultVTKLittleEndian ); \
   checker->setVTKMPIIO( defaultVTKMPIIO ); \
   checker->setVTKForcePVTU( defaultVTKForcePVTU ); \
   return checker;

template< typename Field_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< StabilityChecker< Field_T > > makeStabilityChecker( const Config_T & config,
                                                                const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                                                                const std::string & configBlockName = internal::stabilityCheckerConfigBlock,
                                                                const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                                                const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_STABILITY_CHECKER_CONFIG_PARSER( config )
   using SC_T = StabilityChecker<Field_T>;
   auto checker = shared_ptr< SC_T >( new SC_T( blocks, fieldId, defaultCheckFrequency, defaultOutputToStream, defaultOutputVTK, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_STABILITY_CHECKER_SET_AND_RETURN()
}

template< typename Field_T, typename FlagField_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< StabilityChecker< Field_T, FlagFieldEvaluationFilter<FlagField_T> > >
makeStabilityChecker( const Config_T & config,
                      const weak_ptr< StructuredBlockStorage > & blocks,
                      const ConstBlockDataID & fieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                      const std::string & configBlockName = internal::stabilityCheckerConfigBlock,
                      const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                      const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_STABILITY_CHECKER_CONFIG_PARSER( config )
   using SC_T = StabilityChecker<Field_T, FlagFieldEvaluationFilter<FlagField_T>>;
   auto checker = shared_ptr< SC_T >( new SC_T( blocks, fieldId, FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                                defaultCheckFrequency, defaultOutputToStream, defaultOutputVTK, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_STABILITY_CHECKER_SET_AND_RETURN()
}

template< typename Field_T, typename Filter_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< StabilityChecker< Field_T, Filter_T > >
makeStabilityChecker( const Config_T & config,
                      const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId, const Filter_T & filter,
                      const std::string & configBlockName = internal::stabilityCheckerConfigBlock,
                      const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                      const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_STABILITY_CHECKER_CONFIG_PARSER( config )
   using SC_T = StabilityChecker<Field_T, Filter_T>;
   auto checker = shared_ptr< SC_T >( new SC_T( blocks, fieldId, filter, defaultCheckFrequency, defaultOutputToStream, defaultOutputVTK,
                                                requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_STABILITY_CHECKER_SET_AND_RETURN()
}



#undef WALBERLA_FIELD_MAKE_STABILITY_CHECKER_CONFIG_PARSER
#undef WALBERLA_FIELD_MAKE_STABILITY_CHECKER_SET_AND_RETURN

} // namespace field
} // namespace walberla
