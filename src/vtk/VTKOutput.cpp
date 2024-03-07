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
//! \file VTKOutput.cpp
//! \ingroup vtk
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "VTKOutput.h"
#include "core/debug/Debug.h"
#include "core/math/Vector3.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/MPITextFile.h"
#include "core/mpi/Reduce.h"
#include "core/mpi/Tokenizing.h"
#include "core/mpi/MPIWrapper.h"
#include "core/uid/GlobalState.h"
#include "core/selectable/IsSetSelected.h"

#include <algorithm>
#include <numeric>


namespace walberla {
namespace vtk {



VTKOutput::VTKOutput( const BlockStorage & bs, const std::string & identifier, const uint_t writeFrequency,
                      const std::string & baseFolder, const std::string & executionFolder,
                      const bool continuousNumbering, const bool binary, const bool littleEndian, const bool useMPIIO,
                      const uint_t initialExecutionCount ) :

   unstructuredBlockStorage_( &bs ),
   blockStorage_( nullptr ),
   baseFolder_( baseFolder ), executionFolder_( executionFolder ),
   executionCounter_(initialExecutionCount), initialWriteCallsToSkip_(0), writeFrequency_( writeFrequency ), continuousNumbering_( continuousNumbering ),
   pvdEnd_(-2), binary_( binary ), format_( binary ? std::string("binary") : std::string("ascii") ),
   endianness_( littleEndian ? std::string("LittleEndian") : std::string("BigEndian") ),
   useMPIIO_( useMPIIO ),
   outputDomainDecomposition_( true ),
   samplingDx_( real_c(-1) ), samplingDy_( real_c(-1) ), samplingDz_( real_c(-1) ),
   forcePVTU_( false ), configured_( false ), uniformGrid_( false ), ghostLayers_( uint_c(0) ), writeNextStep_( false )
{
   init( identifier );
}



VTKOutput::VTKOutput( const StructuredBlockStorage & sbs, const std::string & identifier, const uint_t writeFrequency,
                      const std::string & baseFolder, const std::string & executionFolder,
                      const bool continuousNumbering, const bool binary, const bool littleEndian, const bool useMPIIO,
                      const uint_t ghostLayers, const bool forcePVTU, const uint_t initialExecutionCount, const bool amrFileFormat, const bool oneFilePerProcess ) :

   unstructuredBlockStorage_( &sbs.getBlockStorage() ),
   blockStorage_( &sbs ),
   baseFolder_( baseFolder ), executionFolder_( executionFolder ),
   executionCounter_(initialExecutionCount), initialWriteCallsToSkip_(0), writeFrequency_( writeFrequency ), continuousNumbering_( continuousNumbering ),
   pvdEnd_(-2), binary_( binary ), format_( binary ? std::string("binary") : std::string("ascii") ),
   endianness_( littleEndian ? std::string("LittleEndian") : std::string("BigEndian") ),
   useMPIIO_( useMPIIO ),
   outputDomainDecomposition_( false ),
   samplingDx_( real_c(-1) ), samplingDy_( real_c(-1) ), samplingDz_( real_c(-1) ),
   forcePVTU_( forcePVTU ), configured_( false ), uniformGrid_( false ), amrFileFormat_(amrFileFormat), oneFilePerProcess_(oneFilePerProcess), ghostLayers_( ghostLayers ), writeNextStep_( false )
{
   if(ghostLayers > 0 && oneFilePerProcess_)
      WALBERLA_LOG_WARNING_ON_ROOT("Writing out ghostlayers is not supported with oneFilePerProcess. The ghostlayers are just dropped. Alternatively MPI-IO could be used to achieve a similar task")

   if (amrFileFormat && oneFilePerProcess)
   {
      WALBERLA_LOG_WARNING_ON_ROOT("Choose either oneFilePerProcess or amrFileFormat. amrFileFormat is set to false in this combination")
      amrFileFormat_ = false;
   }

   if (useMPIIO_ && amrFileFormat_)
   {
      WALBERLA_LOG_WARNING_ON_ROOT("Choose either MPI-I0 or amrFileFormat. amrFileFormat is set to false in this combination")
      amrFileFormat_ = false;
   }

   if (useMPIIO_ && oneFilePerProcess_)
   {
      WALBERLA_LOG_WARNING_ON_ROOT("Choose either MPI-I0 or oneFilePerProcess. oneFilePerProcess is set to false in this combination")
      oneFilePerProcess_ = false;
   }

   init( identifier );
}



VTKOutput::VTKOutput( const shared_ptr< PointDataSource >& pds, const std::string & identifier, const uint_t writeFrequency,
                      const std::string & baseFolder, const std::string & executionFolder,
                      const bool continuousNumbering, const bool binary, const bool littleEndian, const bool useMPIIO,
                      const uint_t initialExecutionCount ) :

   unstructuredBlockStorage_(nullptr), blockStorage_( nullptr ), pointDataSource_( pds ),
   baseFolder_( baseFolder ), executionFolder_( executionFolder ),
   executionCounter_(initialExecutionCount), initialWriteCallsToSkip_(0), writeFrequency_( writeFrequency ), continuousNumbering_( continuousNumbering ),
   pvdEnd_(-2), binary_( binary ), format_( binary ? std::string("binary") : std::string("ascii") ),
   endianness_( littleEndian ? std::string("LittleEndian") : std::string("BigEndian") ),
   useMPIIO_( useMPIIO ),
   outputDomainDecomposition_( false ),
   samplingDx_( real_c(-1) ), samplingDy_( real_c(-1) ), samplingDz_( real_c(-1) ),
   forcePVTU_( false ), configured_( false ), uniformGrid_( false ), ghostLayers_( uint_c(0) ), writeNextStep_( false )
{
   init( identifier );
}




VTKOutput::VTKOutput( const shared_ptr< PolylineDataSource >& pds, const std::string & identifier, const uint_t writeFrequency,
                      const std::string & baseFolder, const std::string & executionFolder,
                      const bool continuousNumbering, const bool binary, const bool littleEndian, const bool useMPIIO,
                      const uint_t initialExecutionCount ) :

   unstructuredBlockStorage_(nullptr), blockStorage_( nullptr ), polylineDataSource_( pds ),
   baseFolder_( baseFolder ), executionFolder_( executionFolder ),
   executionCounter_(initialExecutionCount), initialWriteCallsToSkip_(0), writeFrequency_( writeFrequency ), continuousNumbering_( continuousNumbering ),
   pvdEnd_(-2), binary_( binary ), format_( binary ? std::string("binary") : std::string("ascii") ),
   endianness_( littleEndian ? std::string("LittleEndian") : std::string("BigEndian") ),
   useMPIIO_( useMPIIO ),
   outputDomainDecomposition_( false ),
   samplingDx_( real_c(-1) ), samplingDy_( real_c(-1) ), samplingDz_( real_c(-1) ),
   forcePVTU_( false ), configured_( false ), uniformGrid_( false ), ghostLayers_( uint_c(0) ), writeNextStep_( false )
{
   init( identifier );
}




void VTKOutput::init( const std::string & identifier )
{
   identifier_ = identifier;

   std::stringstream uniqueIdentifierStrStream;
   uniqueIdentifierStrStream << baseFolder_ << "/" << identifier;
   std::string uniqueIdentifier = uniqueIdentifierStrStream.str();

   if( VTKUID::exists( uniqueIdentifier ) )
      WALBERLA_ABORT( "VTKOutput with identifier \"" << uniqueIdentifier << "\" already exists. VTKOutput objects must have unique identifiers!" );

   VTKUID( uniqueIdentifier, true );

   WALBERLA_ROOT_SECTION()
   {
      filesystem::path path( baseFolder_ + "/" + identifier_ );
      if( filesystem::exists( path ) && executionCounter_ == 0 )
         filesystem::remove_all( path );

      filesystem::path pvd( baseFolder_ + "/" + identifier_ + ".pvd" );
      if( filesystem::exists( pvd ) && executionCounter_ == 0 )
         std::remove( pvd.string().c_str() );

      filesystem::path vthbSeries( baseFolder_ + "/" + identifier_ + ".vthb.series" );
      if( filesystem::exists( vthbSeries ) && executionCounter_ == 0 )
         std::remove( vthbSeries.string().c_str() );

      filesystem::path basePath( baseFolder_ );
      if( !filesystem::exists( basePath ) )
         filesystem::create_directories( basePath );

      if ( !filesystem::exists( path ) )
         filesystem::create_directories( path );
   }

   WALBERLA_MPI_WORLD_BARRIER();
}



VTKOutput::~VTKOutput()
{
   WALBERLA_ROOT_SECTION()
   {
      if( !collectorsToWrite_.empty() )
         writeCollectors( false );
   }
}

void VTKOutput::forceWrite( uint_t number, const bool immediatelyWriteCollectors, const int simultaneousIOOperations,
                            const Set<SUID>& requiredStates, const Set<SUID>& incompatibleStates )
{
   for( auto func = beforeFunctions_.begin(); func != beforeFunctions_.end(); ++func )
      ( *func )( );

   std::ostringstream path;
   path << baseFolder_ << "/" << identifier_ << "/" << executionFolder_ << "_" << number;

   WALBERLA_ROOT_SECTION()
   {
      if( !useMPIIO_ )
      {
         filesystem::path tpath( path.str() );
         if( !filesystem::exists( tpath ) )
            filesystem::create_directory( tpath );
      }
   }
   WALBERLA_MPI_WORLD_BARRIER();


   if( useMPIIO_ )
   {
      bool fileWritten = false;
      if( outputDomainDecomposition_ )
      {
         std::ostringstream oss;
         writeDomainDecompositionPieces( oss, requiredStates, incompatibleStates );
         fileWritten = writeCombinedVTU( oss.str(), number );
      }
      else if( pointDataSource_ )
      {
         std::ostringstream oss;
         writePointDataPieces( oss );
         fileWritten = writeCombinedVTU( oss.str(), number );

      }
      else if( polylineDataSource_ )
      {
         std::ostringstream oss;
         writePolylineDataPieces( oss );
         fileWritten = writeCombinedVTU( oss.str(), number );
      }
      else
      {
         std::ostringstream oss;
         writeBlockPieces( oss, requiredStates, incompatibleStates );

         if( uniformGrid_ )
         {
            if( samplingDx_ <= real_c( 0 ) || samplingDy_ <= real_c( 0 ) || samplingDz_ <= real_c( 0 ) )
               fileWritten = writeCombinedVTI( oss.str(), number );
            else
               fileWritten = writeCombinedVTI_sampled( oss.str(), number );
         }
         else
         {
            fileWritten = writeCombinedVTU( oss.str(), number );
         }
      }
      WALBERLA_ROOT_SECTION()
      {
         if( fileWritten )
            allCollectors_.push_back( number );
      }
   }
   else // no MPIIO
   {
      mpi::Tokenizing tokenizing( simultaneousIOOperations > 0 ? uint_c( simultaneousIOOperations ) : uint_t( 0 ) );
      tokenizing.pre();

      if( outputDomainDecomposition_ )
         writeDomainDecomposition( path.str(), requiredStates, incompatibleStates );
      else if( pointDataSource_ )
         writePointData( path.str() );
      else if( polylineDataSource_ )
         writePolylineData( path.str() );
      else
         writeBlocks( path.str(), requiredStates, incompatibleStates );

      WALBERLA_ROOT_SECTION()
      {
         allCollectors_.push_back( number );
         collectorsToWrite_.push_back( number );
      }

      tokenizing.post();
   }

   if( immediatelyWriteCollectors )
      writeCollectors( true );
}

void VTKOutput::write( const bool immediatelyWriteCollectors, const int simultaneousIOOperations,
                       const Set<SUID>& requiredStates, const Set<SUID>& incompatibleStates )
{
   ++executionCounter_;

   if( executionCounter_ <= initialWriteCallsToSkip_ )
      return;

   if( ( writeFrequency_ == uint_t(0) || ( executionCounter_ - initialWriteCallsToSkip_ - uint_c( 1 ) ) % writeFrequency_ != 0 ) && !writeNextStep_ )
      return;

   uint_t number = continuousNumbering_ ? ( ( executionCounter_ - initialWriteCallsToSkip_ - uint_c( 1 ) ) / writeFrequency_ ) :
                   ( executionCounter_ - uint_c( 1 ) );

   forceWrite(number, immediatelyWriteCollectors, simultaneousIOOperations, requiredStates, incompatibleStates);
   writeNextStep_ = false;
}



void VTKOutput::writeDomainDecomposition( const std::string& path, const Set<SUID>& requiredStates, const Set<SUID>& incompatibleStates ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( unstructuredBlockStorage_ );

   const uint_t numberOfBlocks = unstructuredBlockStorage_->getNumberOfBlocks();

   if( numberOfBlocks == 0 )
      return;

   const int process = MPIManager::instance()->rank();

   std::ostringstream file;
   file << path << "/process_" << process << ".vtu";
   std::ofstream ofs( file.str().c_str() );

   ofs << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
       << " <UnstructuredGrid>\n";

   writeDomainDecompositionPieces( ofs, requiredStates, incompatibleStates );

   ofs << " </UnstructuredGrid>\n"
       << "</VTKFile>" << std::endl;

   ofs.close();
}


void VTKOutput::writeDomainDecompositionPieces( std::ostream& ofs, const Set<SUID>& requiredStates, const Set<SUID>& incompatibleStates ) const
{
   std::vector< const IBlock * > blocks;
   for( auto block = unstructuredBlockStorage_->begin(); block != unstructuredBlockStorage_->end(); ++block )
   {
      if( selectable::isSetSelected( uid::globalState() + block->getState(), requiredStates, incompatibleStates ) )
         blocks.push_back( block.get() );
   }

   const uint_t numberOfBlocks = blocks.size();
   const uint_t points = uint_c( 8 ) * numberOfBlocks;
   const int process = MPIManager::instance()->rank();

   ofs << "  <Piece NumberOfPoints=\"" << points << "\" NumberOfCells=\"" << numberOfBlocks << "\">\n"
       << "   <Points>\n"
       << "    <DataArray type=\"" << vtk::typeToString< float >() << "\" NumberOfComponents=\"3\" format=\"" << format_ << "\">\n";

   std::vector< float > vertex;
   for( auto block = blocks.begin(); block != blocks.end(); ++block )
   {
      const AABB& aabb = (*block)->getAABB();
      vertex.push_back( numeric_cast< float >( aabb.xMin() ) ); vertex.push_back( numeric_cast< float >( aabb.yMin() ) );
      vertex.push_back( numeric_cast< float >( aabb.zMin() ) );
      vertex.push_back( numeric_cast< float >( aabb.xMax() ) ); vertex.push_back( numeric_cast< float >( aabb.yMin() ) );
      vertex.push_back( numeric_cast< float >( aabb.zMin() ) );
      vertex.push_back( numeric_cast< float >( aabb.xMin() ) ); vertex.push_back( numeric_cast< float >( aabb.yMax() ) );
      vertex.push_back( numeric_cast< float >( aabb.zMin() ) );
      vertex.push_back( numeric_cast< float >( aabb.xMax() ) ); vertex.push_back( numeric_cast< float >( aabb.yMax() ) );
      vertex.push_back( numeric_cast< float >( aabb.zMin() ) );
      vertex.push_back( numeric_cast< float >( aabb.xMin() ) ); vertex.push_back( numeric_cast< float >( aabb.yMin() ) );
      vertex.push_back( numeric_cast< float >( aabb.zMax() ) );
      vertex.push_back( numeric_cast< float >( aabb.xMax() ) ); vertex.push_back( numeric_cast< float >( aabb.yMin() ) );
      vertex.push_back( numeric_cast< float >( aabb.zMax() ) );
      vertex.push_back( numeric_cast< float >( aabb.xMin() ) ); vertex.push_back( numeric_cast< float >( aabb.yMax() ) );
      vertex.push_back( numeric_cast< float >( aabb.zMax() ) );
      vertex.push_back( numeric_cast< float >( aabb.xMax() ) ); vertex.push_back( numeric_cast< float >( aabb.yMax() ) );
      vertex.push_back( numeric_cast< float >( aabb.zMax() ) );
   }

   if( binary_ )
   {
      Base64Writer base64;
      for( auto v = vertex.begin(); v != vertex.end(); ++v )
         base64 << *v;
      ofs << "     "; base64.toStream( ofs );
   }
   else for( uint_t i = 0; i != vertex.size(); i += 3 )
      ofs << "     " << vertex[ i ] << " " << vertex[ i + 1 ] << " " << vertex[ i + 2 ] << "\n";

   ofs << "    </DataArray>\n"
      << "   </Points>\n"
      << "   <Cells>\n"
      << "    <DataArray type=\"" << vtk::typeToString< Index >() << "\" Name=\"connectivity\" format=\"" << format_ << "\">\n";

   if( binary_ )
   {
      Base64Writer base64;
      for( int32_t i = 0; i != int32_c( points ); ++i )
         base64 << i;
      ofs << "     "; base64.toStream( ofs );
   }
   else for( int32_t i = 0; i != int32_c( points ); i += 8 )
      ofs << "     " << ( i ) << " " << ( i + 1 ) << " " << ( i + 2 ) << " " << ( i + 3 ) << " "
      << ( i + 4 ) << " " << ( i + 5 ) << " " << ( i + 6 ) << " " << ( i + 7 ) << "\n";

   ofs << "    </DataArray>\n"
      << "    <DataArray type=\"" << vtk::typeToString< Index >() << "\" Name=\"offsets\" format=\"" << format_ << "\">\n";

   if( binary_ )
   {
      Base64Writer base64;
      for( int32_t i = 0; i != int32_c( points ); i += 8 )
         base64 << i + int32_c( 8 );
      ofs << "     "; base64.toStream( ofs );
   }
   else for( int32_t i = 0; i != int32_c( points ); i += 8 )
      ofs << "     " << i + int32_c( 8 ) << "\n";

   ofs << "    </DataArray>\n"
      << "    <DataArray type=\"" << vtk::typeToString< uint8_t >() << "\" Name=\"types\" format=\"" << format_ << "\">\n";

   if( binary_ )
   {
      Base64Writer base64;
      for( uint_t i = 0; i != points; i += 8 )
         base64 << uint8_c( 11 );
      ofs << "     "; base64.toStream( ofs );
   }
   else for( uint_t i = 0; i != points; i += 8 )
      ofs << "     11" << "\n";

   ofs << "    </DataArray>\n"
      << "   </Cells>\n"
      << "   <CellData>\n"
      << "    <DataArray type=\"" << vtk::typeToString< uint8_t >()
      << "\" Name=\"Level\" NumberOfComponents=\"1\" format=\"" << format_ << "\">\n"
      << "     ";

   if( binary_ )
   {
      Base64Writer base64;
      for( auto block = blocks.begin(); block != blocks.end(); ++block )
         base64 << uint8_c( unstructuredBlockStorage_->getLevel( **block ) );
      base64.toStream( ofs );
   }
   else
   {
      for( auto block = blocks.begin(); block != blocks.end(); ++block )
         ofs << unstructuredBlockStorage_->getLevel( **block ) << " ";
      ofs << "\n";
   }

   ofs << "    </DataArray>\n"
      << "    <DataArray type=\"" << vtk::typeToString< int >()
      << "\" Name=\"Process\" NumberOfComponents=\"1\" format=\"" << format_ << "\">\n"
      << "     ";

   if( binary_ )
   {
      Base64Writer base64;
      for( uint_t i = 0; i != numberOfBlocks; ++i )
         base64 << process;
      base64.toStream( ofs );
   }
   else
   {
      for( uint_t i = 0; i != numberOfBlocks; ++i )
         ofs << process << " ";
      ofs << "\n";
   }

   ofs << "    </DataArray>\n"
       << "   </CellData>\n"
       << "  </Piece>\n";
}


void VTKOutput::computeOutputPoints( std::vector<Vector3<real_t> > & points, std::vector<bool> & outputPoint,
                                     uint_t & numberOfPoints ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( pointDataSource_ );

   pointDataSource_->configure();
   points = pointDataSource_->getPoints();

   numberOfPoints = points.size();
   outputPoint.assign( points.size(), true );
   for( uint_t i = 0; i != points.size(); ++i )
   {
      bool included = aabbInclusionFilters_.empty();
      for( auto aabb = aabbInclusionFilters_.begin(); aabb != aabbInclusionFilters_.end() && !included; ++aabb )
         if( aabb->contains( points[ i ][ 0 ], points[ i ][ 1 ], points[ i ][ 2 ] ) ) included = true;
      if( !included )
      {
         outputPoint[ i ] = false;
         --numberOfPoints;
         continue;
      }
      bool excluded = false;
      for( auto aabb = aabbExclusionFilters_.begin(); aabb != aabbExclusionFilters_.end() && !excluded; ++aabb )
         if( aabb->contains( points[ i ][ 0 ], points[ i ][ 1 ], points[ i ][ 2 ] ) ) excluded = true;
      if( excluded )
      {
         outputPoint[ i ] = false;
         --numberOfPoints;
      }
   }
}

void VTKOutput::writePointDataPieceHelper( const std::vector<Vector3<real_t> > & points,
         const std::vector<bool> & outputPoint, const uint_t numberOfPoints, std::ostream & ofs ) const
{
   ofs << "  <Piece NumberOfPoints=\"" << numberOfPoints << "\" NumberOfCells=\"" << numberOfPoints << "\">\n"
       << "   <Points>\n"
       << "    <DataArray type=\"" << vtk::typeToString< float >() << "\" NumberOfComponents=\"3\" format=\"" << format_ << "\">\n";

   if( binary_ )
   {
      Base64Writer base64;
      for( uint_t i = 0; i != points.size(); ++i )
         if( outputPoint[ i ] ) base64 << numeric_cast<float>( points[ i ][ 0 ] ) << numeric_cast<float>( points[ i ][ 1 ] )
            << numeric_cast<float>( points[ i ][ 2 ] );
      ofs << "     "; base64.toStream( ofs );
   }
   else for( uint_t i = 0; i != points.size(); ++i )
      if( outputPoint[ i ] ) ofs << "     " << numeric_cast<float>( points[ i ][ 0 ] ) << " " << numeric_cast<float>( points[ i ][ 1 ] )
         << " " << numeric_cast<float>( points[ i ][ 2 ] ) << "\n";

   ofs << "    </DataArray>\n"
       << "   </Points>\n"
       << "   <Cells>\n"
       << "    <DataArray type=\"" << vtk::typeToString< Index >() << "\" Name=\"connectivity\" format=\"" << format_ << "\">\n";

   Index j = 0;
   if( binary_ )
   {
      Base64Writer base64;
      for( uint_t i = 0; i != points.size(); ++i )
         if( outputPoint[ i ] ) base64 << j++;
      ofs << "     "; base64.toStream( ofs );
   }
   else for( uint_t i = 0; i != points.size(); ++i )
      if( outputPoint[ i ] ) ofs << "     " << j++ << "\n";

   ofs << "    </DataArray>\n"
       << "    <DataArray type=\"" << vtk::typeToString< Index >() << "\" Name=\"offsets\" format=\"" << format_ << "\">\n";

   j = 0;
   if( binary_ )
   {
      Base64Writer base64;
      for( uint_t i = 0; i != points.size(); ++i )
         if( outputPoint[ i ] ) base64 << ++j;
      ofs << "     "; base64.toStream( ofs );
   }
   else for( uint_t i = 0; i != points.size(); ++i )
      if( outputPoint[ i ] ) ofs << "     " << ++j << "\n";

   ofs << "    </DataArray>\n"
       << "    <DataArray type=\"" << vtk::typeToString< uint8_t >() << "\" Name=\"types\" format=\"" << format_ << "\">\n";

   if( binary_ )
   {
      Base64Writer base64;
      for( uint_t i = 0; i != points.size(); ++i )
         if( outputPoint[ i ] ) base64 << uint8_c( 1 );
      ofs << "     "; base64.toStream( ofs );
   }
   else for( uint_t i = 0; i != points.size(); ++i )
      if( outputPoint[ i ] ) ofs << "     1" << "\n";

   ofs << "    </DataArray>\n"
       << "   </Cells>\n"
       << "   <PointData>\n";

   uint_t dataIndex = 0;
   auto attributes = pointDataSource_->getAttributes();
   for( auto dataArray = attributes.begin(); dataArray != attributes.end(); ++dataArray, ++dataIndex )
   {
      const uint_t components = dataArray->components;
      WALBERLA_ASSERT_GREATER( components, 0 );

      ofs << "    <DataArray type=\"" << dataArray->type << "\" Name=\"" << dataArray->name << "\" NumberOfComponents=\""
         << components << "\" format=\"" << format_ << "\">\n";

      if( binary_ )
      {
         Base64Writer base64;
         for( uint_t i = 0; i != points.size(); ++i )
            if( outputPoint[ i ] )
               for( uint_t c = 0; c != components; ++c )
                  pointDataSource_->push( base64, dataIndex, i, c );
         ofs << "     "; base64.toStream( ofs );
      }
      else for( uint_t i = 0; i != points.size(); ++i )
         if( outputPoint[ i ] )
         {
         ofs << "    ";
         for( uint_t c = 0; c != components; ++c )
         {
            ofs << " ";
            pointDataSource_->push( ofs, dataIndex, i, c );
         }
         ofs << "\n";
         }

      ofs << "    </DataArray>\n";
   }

   ofs << "   </PointData>\n"
       << "  </Piece>\n";
}






void VTKOutput::writePointData( const std::string& path )
{
   configured_ = true;

   std::vector<Vector3<real_t> > points;
   std::vector<bool> outputPoint;
   uint_t numberOfPoints;

   computeOutputPoints( points, outputPoint, numberOfPoints );

   if( numberOfPoints == 0 )
      return;

   std::ostringstream file;
   file << path << "/process_" << MPIManager::instance()->rank() << ".vtu";
   std::ofstream ofs( file.str().c_str() );

   ofs << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
       << " <UnstructuredGrid>\n";


   writePointDataPieceHelper( points, outputPoint, numberOfPoints, ofs );

   ofs << " </UnstructuredGrid>\n"
       << "</VTKFile>" << std::endl;

   ofs.close();
}


void VTKOutput::writePointDataPieces(std::ostream& ofs)
{
   configured_ = true;

   std::vector<Vector3<real_t> > points;
   std::vector<bool> outputPoint;
   uint_t numberOfPoints;

   computeOutputPoints( points, outputPoint, numberOfPoints );

   //if( numberOfPoints == 0 )
   //   return;

   writePointDataPieceHelper( points, outputPoint, numberOfPoints, ofs );
}


void VTKOutput::computeOutputPolylines( std::vector< std::vector< Vector3< real_t > > > & lines,
   std::vector< std::vector< bool > > & outputPolylinePoint, std::vector< size_t > & polylineSize,
   uint_t & numberOfPolylines, uint_t & numberOfPolylinePoints ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( polylineDataSource_ );

   polylineDataSource_->configure();
   lines = polylineDataSource_->getPolyLines();

   numberOfPolylines = 0u;
   numberOfPolylinePoints = 0;
   outputPolylinePoint.resize( lines.size() );
   polylineSize.resize( lines.size() );

   auto outputPolylinePointIt = outputPolylinePoint.begin();
   auto polylineSizeIt = polylineSize.begin();
   for( auto polylineIt = lines.begin(); polylineIt != lines.end(); ++polylineIt, ++outputPolylinePointIt, ++polylineSizeIt )
   {
      numberOfPolylinePoints += polylineIt->size();
      outputPolylinePointIt->assign( polylineIt->size(), true );
      *polylineSizeIt = polylineIt->size();
   }

   outputPolylinePointIt = outputPolylinePoint.begin();
   polylineSizeIt = polylineSize.begin();
   for( auto polylineIt = lines.begin(); polylineIt != lines.end(); ++polylineIt, ++outputPolylinePointIt, ++polylineSizeIt )
   {
      auto outputPointIt = outputPolylinePointIt->begin();
      for( auto pointIt = polylineIt->begin(); pointIt != polylineIt->end(); ++pointIt, ++outputPointIt )
      {
         bool included = aabbInclusionFilters_.empty();
         for( auto aabb = aabbInclusionFilters_.begin(); aabb != aabbInclusionFilters_.end() && !included; ++aabb )
            if( aabb->contains( ( *pointIt )[ 0 ], ( *pointIt )[ 1 ], ( *pointIt )[ 2 ] ) ) included = true;
         if( !included )
         {
            *outputPointIt = false;
            --numberOfPolylinePoints;
            --( *polylineSizeIt );
            continue;
         }
         bool excluded = false;
         for( auto aabb = aabbExclusionFilters_.begin(); aabb != aabbExclusionFilters_.end() && !excluded; ++aabb )
            if( aabb->contains( ( *pointIt )[ 0 ], ( *pointIt )[ 1 ], ( *pointIt )[ 2 ] ) ) excluded = true;
         if( excluded )
         {
            *outputPointIt = false;
            --numberOfPolylinePoints;
            --( *polylineSizeIt );
         }
      }

      if( *polylineSizeIt > 0u )
         ++numberOfPolylines;
   }
}


void VTKOutput::writePolylineDataPieceHelper( const std::vector< std::vector< Vector3< real_t > > > & lines,
   const std::vector< std::vector< bool > > & outputPolylinePoint, const std::vector< size_t > & polylineSize,
   const uint_t numberOfPolylines, const uint_t numberOfPolylinePoints, std::ostream & ofs ) const
{
   ofs << "  <Piece NumberOfPoints=\"" << numberOfPolylinePoints << "\" NumberOfCells=\"" << numberOfPolylines << "\">\n"
       << "   <Points>\n"
       << "    <DataArray type=\"" << vtk::typeToString< float >() << "\" NumberOfComponents=\"3\" format=\"" << format_ << "\">\n";

   if( binary_ )
   {
      Base64Writer base64;

      auto outputPolylinePointIt = outputPolylinePoint.begin();
      for( auto polylineIt = lines.begin(); polylineIt != lines.end(); ++polylineIt, ++outputPolylinePointIt )
      {
         auto outputPointIt = outputPolylinePointIt->begin();
         for( auto pointIt = polylineIt->begin(); pointIt != polylineIt->end(); ++pointIt, ++outputPointIt )
         {
            if( *outputPointIt ) base64 << numeric_cast< float >( ( *pointIt )[ 0 ] ) << numeric_cast< float >( ( *pointIt )[ 1 ] )
               << numeric_cast< float >( ( *pointIt )[ 2 ] );
         }
      }

      ofs << "     "; base64.toStream( ofs );
   }
   else
   {
      auto outputPolylinePointIt = outputPolylinePoint.begin();
      for( auto polylineIt = lines.begin(); polylineIt != lines.end(); ++polylineIt, ++outputPolylinePointIt )
      {
         auto outputPointIt = outputPolylinePointIt->begin();
         for( auto pointIt = polylineIt->begin(); pointIt != polylineIt->end(); ++pointIt, ++outputPointIt )
         {
            if( *outputPointIt ) ofs << "     " << numeric_cast< float >( ( *pointIt )[ 0 ] ) << " " << numeric_cast< float >( ( *pointIt )[ 1 ] )
               << " " << numeric_cast< float >( ( *pointIt )[ 2 ] ) << "\n";
         }
      }
   }

   ofs << "    </DataArray>\n"
       << "   </Points>\n"
       << "   <Cells>\n"
       << "    <DataArray type=\"" << vtk::typeToString< Index >() << "\" Name=\"connectivity\" format=\"" << format_ << "\">\n";

   Index j = 0;
   if( binary_ )
   {
      Base64Writer base64;

      auto outputPolylinePointIt = outputPolylinePoint.begin();
      for( auto polylineIt = lines.begin(); polylineIt != lines.end(); ++polylineIt, ++outputPolylinePointIt )
      {
         auto outputPointIt = outputPolylinePointIt->begin();
         for( auto pointIt = polylineIt->begin(); pointIt != polylineIt->end(); ++pointIt, ++outputPointIt )
         {
            if( *outputPointIt ) base64 << j++;
         }
      }

      ofs << "     "; base64.toStream( ofs );

   }
   else
   {
      auto outputPolylinePointIt = outputPolylinePoint.begin();
      for( auto polylineIt = lines.begin(); polylineIt != lines.end(); ++polylineIt, ++outputPolylinePointIt )
      {
         auto outputPointIt = outputPolylinePointIt->begin();
         for( auto pointIt = polylineIt->begin(); pointIt != polylineIt->end(); ++pointIt, ++outputPointIt )
         {
            if( *outputPointIt ) ofs << "     " << j++ << "\n";
         }
      }
   }

   ofs << "    </DataArray>\n"
       << "    <DataArray type=\"" << vtk::typeToString< Index >() << "\" Name=\"offsets\" format=\"" << format_ << "\">\n";

   j = 0;
   if( binary_ )
   {
      Base64Writer base64;

      auto polylineSizeIt = polylineSize.begin();
      for( auto polylineIt = lines.begin(); polylineIt != lines.end(); ++polylineIt, ++polylineSizeIt )
      {
         if( *polylineSizeIt > 0u )
         {
            j += numeric_cast< Index > ( *polylineSizeIt );
            base64 << j;
         }
      }

      ofs << "     "; base64.toStream( ofs );
   }
   else
   {
      auto polylineSizeIt = polylineSize.begin();
      for( auto polylineIt = lines.begin(); polylineIt != lines.end(); ++polylineIt, ++polylineSizeIt )
      {
         if( *polylineSizeIt > 0u )
         {
            j += numeric_cast< Index > ( *polylineSizeIt );
            ofs << "     " << j << "\n";
         }
      }

   }

   ofs << "    </DataArray>\n"
       << "    <DataArray type=\"" << vtk::typeToString< uint8_t >() << "\" Name=\"types\" format=\"" << format_ << "\">\n";

   if( binary_ )
   {
      Base64Writer base64;

      auto polylineSizeIt = polylineSize.begin();
      for( auto polylineIt = lines.begin(); polylineIt != lines.end(); ++polylineIt, ++polylineSizeIt )
         if( *polylineSizeIt > 0u )
         {
         base64 << uint8_c( 4 );
         }

      ofs << "     "; base64.toStream( ofs );
   }
   else
   {
      auto polylineSizeIt = polylineSize.begin();
      for( auto polylineIt = lines.begin(); polylineIt != lines.end(); ++polylineIt )
         if( *polylineSizeIt > 0u )
         {
         ofs << "     4" << "\n";
         }
   }

   ofs << "    </DataArray>\n"
       << "   </Cells>\n"
       << "   <PointData>\n";

   uint_t dataIndex = 0;
   auto attributes = polylineDataSource_->getAttributes();
   for( auto dataArray = attributes.begin(); dataArray != attributes.end(); ++dataArray, ++dataIndex )
   {
      const uint_t components = dataArray->components;
      WALBERLA_ASSERT_GREATER( components, 0 );

      ofs << "    <DataArray type=\"" << dataArray->type << "\" Name=\"" << dataArray->name << "\" NumberOfComponents=\""
         << components << "\" format=\"" << format_ << "\">\n";

      if( binary_ )
      {
         Base64Writer base64;

         for( uint_t polylineIdx = 0; polylineIdx < lines.size(); ++polylineIdx )
            for( uint_t polylinePointIdx = 0; polylinePointIdx < lines[ polylineIdx ].size(); ++polylinePointIdx )
            {
            if( outputPolylinePoint[ polylineIdx ][ polylinePointIdx ] )
               for( uint_t c = 0; c != components; ++c )
                  polylineDataSource_->push( base64, dataIndex, polylineIdx, polylinePointIdx, c );
            }

         ofs << "     "; base64.toStream( ofs );
      }
      else
      {
         for( uint_t polylineIdx = 0; polylineIdx < lines.size(); ++polylineIdx )
            for( uint_t polylinePointIdx = 0; polylinePointIdx < lines[ polylineIdx ].size(); ++polylinePointIdx )
            {
            if( outputPolylinePoint[ polylineIdx ][ polylinePointIdx ] )
            {
               ofs << "    ";
               for( uint_t c = 0; c != components; ++c )
               {
                  ofs << " ";
                  polylineDataSource_->push( ofs, dataIndex, polylineIdx, polylinePointIdx, c );
               }
               ofs << "\n";
            }
            }
      }

      ofs << "    </DataArray>\n";
   }

   ofs << "   </PointData>\n"
       << "  </Piece>\n";
}


void VTKOutput::writePolylineData( const std::string& path )
{
   configured_ = true;

   std::vector< std::vector< Vector3< real_t > > > lines;
   uint_t numberOfPolylines;
   uint_t numberOfPolylinePoints;
   std::vector< std::vector< bool > > outputPolylinePoint;
   std::vector< size_t > polylineSize;

   computeOutputPolylines( lines, outputPolylinePoint, polylineSize, numberOfPolylines, numberOfPolylinePoints );

   if( numberOfPolylinePoints == 0 )
      return;

   std::ostringstream file;
   file << path << "/process_" << MPIManager::instance()->rank() << ".vtu";
   std::ofstream ofs( file.str().c_str() );

   ofs << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
      << " <UnstructuredGrid>\n";

   writePolylineDataPieceHelper( lines, outputPolylinePoint, polylineSize, numberOfPolylines, numberOfPolylinePoints, ofs );

   ofs << " </UnstructuredGrid>\n"
       << "</VTKFile>" << std::endl;

   ofs.close();
}




void VTKOutput::writePolylineDataPieces(std::ostream& ofs)
{
   configured_ = true;

   std::vector< std::vector< Vector3< real_t > > > lines;
   uint_t numberOfPolylines;
   uint_t numberOfPolylinePoints;
   std::vector< std::vector< bool > > outputPolylinePoint;
   std::vector< size_t > polylineSize;

   computeOutputPolylines( lines, outputPolylinePoint, polylineSize, numberOfPolylines, numberOfPolylinePoints );

   if (numberOfPolylinePoints == 0)
      return;

   writePolylineDataPieceHelper( lines, outputPolylinePoint, polylineSize, numberOfPolylines, numberOfPolylinePoints, ofs );
}


void VTKOutput::computeVTUCells( const IBlock& block, CellVector & cellsOut ) const
{
   const cell_idx_t gl = cell_idx_c( ghostLayers_ );
   const cell_idx_t begin = cell_idx_c( -1 ) * gl;

   CellSet includedCells;

   if( cellInclusionFunctions_.empty() )
   {
      for( cell_idx_t z = begin; z < cell_idx_c( blockStorage_->getNumberOfZCells( block ) ) + gl; ++z )
         for( cell_idx_t y = begin; y < cell_idx_c( blockStorage_->getNumberOfYCells( block ) ) + gl; ++y )
            for( cell_idx_t x = begin; x < cell_idx_c( blockStorage_->getNumberOfXCells( block ) ) + gl; ++x )
               includedCells.insert( x, y, z );
   }
   else if( cellInclusionFunctions_.size() == 1 ) {
      cellInclusionFunctions_[ 0 ]( includedCells, block, *blockStorage_, ghostLayers_ );
   }
   else {
      for( auto func = cellInclusionFunctions_.begin(); func != cellInclusionFunctions_.end(); ++func ) {
         CellSet cellSet;
         ( *func )( cellSet, block, *blockStorage_, ghostLayers_ );
         includedCells += cellSet;
      }
   }

   CellSet excludedCells;

   if( cellExclusionFunctions_.size() == 1 ) {
      cellExclusionFunctions_[ 0 ]( excludedCells, block, *blockStorage_, ghostLayers_ );
   }
   else {
      for( auto func = cellExclusionFunctions_.begin(); func != cellExclusionFunctions_.end(); ++func ) {
         CellSet cellSet;
         ( *func )( cellSet, block, *blockStorage_, ghostLayers_ );
         excludedCells += cellSet;
      }
   }


   for( cell_idx_t z = begin; z < cell_idx_c( blockStorage_->getNumberOfZCells( block ) ) + gl; ++z )
      for( cell_idx_t y = begin; y < cell_idx_c( blockStorage_->getNumberOfYCells( block ) ) + gl; ++y )
         for( cell_idx_t x = begin; x < cell_idx_c( blockStorage_->getNumberOfXCells( block ) ) + gl; ++x )
            if( includedCells.contains( x, y, z ) && !excludedCells.contains( x, y, z ) ) cellsOut.push_back( x, y, z );
}


void VTKOutput::writeBlocks( const std::string& path, const Set<SUID>& requiredStates, const Set<SUID>& incompatibleStates )
{
   WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ );

   if( !configured_ ) {
      if( !forcePVTU_ && cellInclusionFunctions_.empty() && cellExclusionFunctions_.empty() &&
          blockStorage_->getNumberOfLevels() == 1 && ghostLayers_ == 0 ) // uniform data -> vti
      {
            uniformGrid_ = true;
      }
      configured_ = true;
   }

   if(!uniformGrid_ && oneFilePerProcess_)
   {
      const int rank         = MPIManager::instance()->rank();
      std::ostringstream file;
      file << path << "/dataRank[" << rank << "].vtu";
      std::ofstream ofs(file.str().c_str());
      writeParallelVTU( ofs, requiredStates, incompatibleStates );
      ofs.close();
   }
   else
   {
      std::vector< const IBlock* > blocks;
      for( auto block = blockStorage_->begin(); block != blockStorage_->end(); ++block )
      {
         if( selectable::isSetSelected( uid::globalState() + block->getState(), requiredStates, incompatibleStates ) )
            blocks.push_back( block.get() );
      }
      for( auto it = blocks.begin(); it != blocks.end(); ++it )
      {
         WALBERLA_ASSERT_NOT_NULLPTR( *it );
         const IBlock& block = **it;
         const uint_t level = blockStorage_->getLevel(block);

         std::ostringstream file;
         file << path << "/block [" << block.getId() << "] level[" << level << "].";

         if( uniformGrid_ || amrFileFormat_ ) // uniform data -> vti  amr data -> vti
         {
            file << "vti";
            std::ofstream ofs( file.str().c_str()  );
            if( samplingDx_ <= real_c(0) || samplingDy_ <= real_c(0) || samplingDz_ <= real_c(0) )
               writeVTI( ofs, block );
            else
               writeVTI_sampling( ofs, block );
            ofs.close();
         }
         else // unstructured data -> vtu
         {
            CellVector cells; // cells to be written to file
            computeVTUCells( block, cells );

            if( !cells.empty() )
            {
               file << "vtu";
               std::ofstream ofs( file.str().c_str()  );
               if( samplingDx_ <= real_c(0) || samplingDy_ <= real_c(0) || samplingDz_ <= real_c(0) )
                  writeVTU( ofs, block, cells );
               else
                  writeVTU_sampling( ofs, block, cells );
               ofs.close();
            }
         }
      }
   }
}



void VTKOutput::writeBlockPieces( std::ostream & oss, const Set<SUID>& requiredStates, const Set<SUID>& incompatibleStates )
{
   WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ );

   std::vector< const IBlock* > blocks;
   for( auto block = blockStorage_->begin(); block != blockStorage_->end(); ++block )
   {
      if( selectable::isSetSelected( uid::globalState() + block->getState(), requiredStates, incompatibleStates ) )
         blocks.push_back( block.get() );
   }

   if( !configured_ ) {
      if( !forcePVTU_ && cellInclusionFunctions_.empty() && cellExclusionFunctions_.empty() &&
         blockStorage_->getNumberOfLevels() == 1 && ghostLayers_ == 0 ) // uniform data -> vti
      {
         uniformGrid_ = true;
      }
      configured_ = true;
   }

   for( auto it = blocks.begin(); it != blocks.end(); ++it )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( *it );
      const IBlock& block = **it;

      if( uniformGrid_ ) // uniform data -> vti
      {
         if( samplingDx_ <= real_c(0) || samplingDy_ <= real_c(0) || samplingDz_ <= real_c(0) )
            writeVTIPiece( oss, block );
         else
            writeVTIPiece_sampling( oss, block );
      }
      else // unstructured data -> vtu
      {
         CellVector cells; // cells to be written to file
         computeVTUCells( block, cells );

         if( !cells.empty() )
         {
            if( samplingDx_ <= real_c(0) || samplingDy_ <= real_c(0) || samplingDz_ <= real_c(0) )
               writeVTUPiece( oss, block, cells );
            else
               writeVTUPiece_sampling( oss, block, cells );
         }
      }
   }
}



void VTKOutput::writeVTI( std::ostream& ofs, const IBlock& block ) const
{
   const CellInterval& cellBB = blockStorage_->getBlockCellBB( block );
   const AABB&         domain = blockStorage_->getDomain();
   const uint_t        level  = blockStorage_->getLevel( block );

   ofs << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
       << " <ImageData WholeExtent=\"" << cellBB.xMin() << " " << ( cellBB.xMax() + 1 ) << " "
       << cellBB.yMin() << " " << ( cellBB.yMax() + 1 ) << " "
       << cellBB.zMin() << " " << ( cellBB.zMax() + 1 ) << "\""
       << " Origin=\"" << domain.xMin() << " " << domain.yMin() << " " << domain.zMin() << "\""
       << " Spacing=\"" << blockStorage_->dx(level) << " " << blockStorage_->dy(level) << " " << blockStorage_->dz(level) << "\">\n";

   writeVTIPiece( ofs, block );

   ofs << " </ImageData>\n"
       << "</VTKFile>" << std::endl;
}

void VTKOutput::writeVTIPiece( std::ostream& ofs, const IBlock& block ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ );
   WALBERLA_ASSERT_EQUAL( ghostLayers_, 0 );

   const CellInterval& cellBB = blockStorage_->getBlockCellBB( block );

   ofs << "  <Piece Extent=\"" << cellBB.xMin() << " " << ( cellBB.xMax() + 1 ) << " "
       << cellBB.yMin() << " " << ( cellBB.yMax() + 1 ) << " "
       << cellBB.zMin() << " " << ( cellBB.zMax() + 1 ) << "\">\n"
       << "   <CellData>\n";

   CellVector cells;

   for( uint_t z = 0; z < blockStorage_->getNumberOfZCells( block ); ++z )
      for( uint_t y = 0; y < blockStorage_->getNumberOfYCells( block ); ++y )
         for( uint_t x = 0; x < blockStorage_->getNumberOfXCells( block ); ++x )
            cells.push_back( x, y, z );

   writeCellData( ofs, block, cells );

   ofs << "   </CellData>\n"
      << "  </Piece>\n";
}

void VTKOutput::writeVTI_sampling( std::ostream& ofs, const IBlock& block ) const
{
   const AABB&  blockBB = block.getAABB();
   const AABB&  domain  = blockStorage_->getDomain();

   CellInterval cellBB = getSampledCellInterval( blockBB );

   ofs << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
       << " <ImageData WholeExtent=\"" << cellBB.xMin() << " " << ( cellBB.xMax() + 1 ) << " "
       << cellBB.yMin() << " " << ( cellBB.yMax() + 1 ) << " "
       << cellBB.zMin() << " " << ( cellBB.zMax() + 1 ) << "\""
       << " Origin=\"" << domain.xMin() << " " << domain.yMin() << " " << domain.zMin() << "\""
       << " Spacing=\"" << samplingDx_ << " " << samplingDy_ << " " << samplingDz_ << "\">\n";

   writeVTIPiece_sampling( ofs, block );

   ofs << " </ImageData>\n"
       << "</VTKFile>" << std::endl;
}

void VTKOutput::writeVTIPiece_sampling( std::ostream& ofs, const IBlock& block ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ );
   WALBERLA_ASSERT_EQUAL( ghostLayers_, 0 );

   WALBERLA_ASSERT_GREATER( samplingDx_, real_t(0) );
   WALBERLA_ASSERT_GREATER( samplingDy_, real_t(0) );
   WALBERLA_ASSERT_GREATER( samplingDz_, real_t(0) );

   const AABB&  blockBB = block.getAABB();
   const AABB&  domain = blockStorage_->getDomain();
   const uint_t level = blockStorage_->getLevel( block );

   CellInterval cellBB = getSampledCellInterval( blockBB );

   ofs << "  <Piece Extent=\"" << cellBB.xMin() << " " << (cellBB.xMax()+1) << " "
       << cellBB.yMin() << " " << (cellBB.yMax()+1) << " "
       << cellBB.zMin() << " " << (cellBB.zMax()+1) << "\">\n"
       << "   <CellData>\n";

   CellVector cells;
   for( auto it = cellBB.begin(); it != cellBB.end(); ++it )
   {
      Vector3<real_t> world( domain.xMin() + ( real_c( it->x() ) + real_t(0.5) ) * samplingDx_,
         domain.yMin() + ( real_c( it->y() ) + real_t(0.5) ) * samplingDy_,
         domain.zMin() + ( real_c( it->z() ) + real_t(0.5) ) * samplingDz_ );

      Cell cell;
      blockStorage_->getCell( cell, world[0], world[1], world[2], level );

      blockStorage_->transformGlobalToBlockLocalCell( cell, block );
      cells.push_back( cell );
   }

   writeCellData( ofs, block, cells );

   ofs << "   </CellData>\n"
       << "  </Piece>\n";
}



void VTKOutput::writeVTU( std::ostream& ofs, const IBlock& block, const CellVector& cells ) const
{
   ofs << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
       << " <UnstructuredGrid>\n";

   writeVTUPiece( ofs, block, cells );

   ofs << " </UnstructuredGrid>\n"
       << "</VTKFile>" << std::endl;
}

void VTKOutput::writeParallelVTU( std::ostream& ofs, const Set<SUID>& requiredStates, const Set<SUID>& incompatibleStates  ) const
{
   ofs << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
       << " <UnstructuredGrid>\n";

   WALBERLA_ASSERT_NOT_NULLPTR(blockStorage_);
   const uint_t finestLevel = blockStorage_->getNumberOfLevels() - 1;

   std::map< Vertex, Index, VertexCompare > vimap; // vertex<->index map
   std::vector< VertexCoord >               vc;    // vertex coordinates
   std::vector< Index >                     ci;    // ci[0] to ci[7]: indices for cell number one, ci[8] to ci[15]: ...
   uint_t numberOfCells = 0;

   for( auto block = blockStorage_->begin(); block != blockStorage_->end(); ++block )
   {
      if( !selectable::isSetSelected( uid::globalState() + block->getState(), requiredStates, incompatibleStates ) )
         continue;

      // CellVector cells; // cells to be written to file
      // computeVTUCells( *block, cells );

      const uint_t level = blockStorage_->getLevel(*block);
      const cell_idx_t factorToFinest = 1 << (finestLevel - level);
      const CellInterval cells = blockStorage_->getBlockCellBB(*block); //  These are global cells

      for (auto cell = cells.begin(); cell != cells.end(); ++cell)
      {
         numberOfCells++;
         const AABB aabb = blockStorage_->getCellAABB(*cell, level);
         for (cell_idx_t z = 0; z != 2; ++z) {
            for (cell_idx_t y = 0; y != 2; ++y) {
               for (cell_idx_t x = 0; x != 2; ++x)
               {
                  const Vertex v((cell->x() + x) * factorToFinest, (cell->y() + y) * factorToFinest, (cell->z() + z) * factorToFinest);
                  auto mapping = vimap.find(v);
                  if (mapping != vimap.end()) // vertex already exists
                  {
                     ci.push_back(mapping->second);
                  }
                  else // new vertex
                  {
                     vimap[v] = numeric_cast< Index >(vc.size());
                     ci.push_back(numeric_cast< Index >(vc.size()));
                     vc.emplace_back((x == 0) ? aabb.xMin() : aabb.xMax(),
                                     (y == 0) ? aabb.yMin() : aabb.yMax(),
                                     (z == 0) ? aabb.zMin() : aabb.zMax());
                  }
               }
            }
         }
      }
   }
   // <--- setting up vertex-index mapping
   writeVTUHeaderPiece(ofs, numberOfCells, vc, ci);

   ofs << "   <CellData>\n";

   writeCellData(ofs, requiredStates, incompatibleStates);

   ofs << "   </CellData>\n"
       << "  </Piece>\n";

   ofs << " </UnstructuredGrid>\n"
       << "</VTKFile>" << std::endl;

}



void VTKOutput::writeVTUPiece( std::ostream& ofs, const IBlock& block, const CellVector& cells ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR(blockStorage_);

   // setting up vertex-index mapping --->

   std::map< Vertex, Index, VertexCompare > vimap; // vertex<->index map
   std::vector< VertexCoord >               vc;    // vertex coordinates
   std::vector< Index >                     ci;    // ci[0] to ci[7]: indices for cell number one, ci[8] to ci[15]: ...

   for (auto cell = cells.begin(); cell != cells.end(); ++cell)
   {
      AABB aabb;
      blockStorage_->getBlockLocalCellAABB(block, *cell, aabb);

      for (cell_idx_t z = 0; z != 2; ++z) {
         for (cell_idx_t y = 0; y != 2; ++y) {
            for (cell_idx_t x = 0; x != 2; ++x)
            {
               Vertex v(cell->x() + x, cell->y() + y, cell->z() + z);

               auto mapping = vimap.find(v);
               if (mapping != vimap.end()) // vertex already exists
               {
                  ci.push_back(mapping->second);
               }
               else // new vertex
               {
                  vimap[v] = numeric_cast< Index >(vc.size());
                  ci.push_back(numeric_cast< Index >(vc.size()));
                  vc.emplace_back((x == 0) ? aabb.xMin() : aabb.xMax(),
                     (y == 0) ? aabb.yMin() : aabb.yMax(),
                     (z == 0) ? aabb.zMin() : aabb.zMax());
               }
            }
         }
      }
   }

   WALBERLA_ASSERT_EQUAL(ci.size(), 8 * cells.size());

   // <--- setting up vertex-index mapping

   writeVTUHeaderPiece(ofs, uint_c(cells.size()), vc, ci);

   ofs << "   <CellData>\n";

   if (ghostLayers_ > 0)
   {
      ofs << "    <DataArray type=\"" << vtk::typeToString< uint8_t >()
          << "\" Name=\"vtkGhostLevels\" NumberOfComponents=\"1\" format=\"" << format_ << "\">\n"
          << "     ";

      if (binary_)
      {
         Base64Writer base64;
         for (auto cell = cells.begin(); cell != cells.end(); ++cell)
            base64 << ghostLayerNr(block, cell->x(), cell->y(), cell->z());
         base64.toStream(ofs);
      }
      else
      {
         for (auto cell = cells.begin(); cell != cells.end(); ++cell)
            ofs << uint_c(ghostLayerNr(block, cell->x(), cell->y(), cell->z())) << " ";
         ofs << "\n";
      }

      ofs << "    </DataArray>\n";
   }

   writeCellData(ofs, block, cells);

   ofs << "   </CellData>\n"
       << "  </Piece>\n";
}



void VTKOutput::writeVTU_sampling( std::ostream& ofs, const IBlock& block, const CellVector& localCells ) const
{
   ofs << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
       << " <UnstructuredGrid>\n";

   writeVTUPiece_sampling( ofs, block, localCells );

   ofs << " </UnstructuredGrid>\n"
       << "</VTKFile>" << std::endl;

}



void VTKOutput::writeVTUPiece_sampling(std::ostream& ofs, const IBlock& block, const CellVector& localCells) const
{
   std::vector< SamplingCell > cells = getSamplingCells(block, localCells);

   // setting up vertex-index mapping --->

   std::map< Vertex, Index, VertexCompare > vimap; // vertex<->index map
   std::vector< VertexCoord >               vc;    // vertex coordinates
   std::vector< Index >                     ci;    // ci[0] to ci[7]: indices for cell number one, ci[8] to ci[15]: ...

   for (auto cell = cells.begin(); cell != cells.end(); ++cell) {
      for (cell_idx_t z = 0; z != 2; ++z) {
         for (cell_idx_t y = 0; y != 2; ++y) {
            for (cell_idx_t x = 0; x != 2; ++x)
            {
               Vertex v(cell->coordinates_.x() + x, cell->coordinates_.y() + y, cell->coordinates_.z() + z);

               auto mapping = vimap.find(v);
               if (mapping != vimap.end()) // vertex already exists
               {
                  ci.push_back(mapping->second);
               }
               else // new vertex
               {
                  vimap[v] = numeric_cast< Index >(vc.size());
                  ci.push_back(numeric_cast< Index >(vc.size()));
                  vc.emplace_back((x == 0) ? cell->aabb_.xMin() : cell->aabb_.xMax(),
                     (y == 0) ? cell->aabb_.yMin() : cell->aabb_.yMax(),
                     (z == 0) ? cell->aabb_.zMin() : cell->aabb_.zMax());
               }
            }
         }
      }
   }

   WALBERLA_ASSERT_EQUAL(ci.size(), 8 * cells.size());

   // <--- setting up vertex-index mapping

   writeVTUHeaderPiece(ofs, uint_c(cells.size()), vc, ci);

   ofs << "   <CellData>\n";

   writeCellData(ofs, block, cells);

   ofs << "   </CellData>\n"
      << "  </Piece>\n";
}



void VTKOutput::writeVTUHeader( std::ofstream& ofs, const uint_t numberOfCells,
                                const std::vector< VertexCoord > & vc, const std::vector< Index > & ci ) const
{
   ofs << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
       << " <UnstructuredGrid>\n";

   writeVTUHeaderPiece( ofs, numberOfCells, vc, ci );
}



void VTKOutput::writeVTUHeaderPiece( std::ostream& ofs, const uint_t numberOfCells,
                                     const std::vector< VertexCoord > & vc, const std::vector< Index > & ci ) const
{
   ofs << "  <Piece NumberOfPoints=\"" << vc.size() << "\" NumberOfCells=\"" << numberOfCells << "\">\n"
       << "   <Points>\n"
       << "    <DataArray type=\"" << vtk::typeToString< float >() << "\" NumberOfComponents=\"3\" format=\"" << format_ << "\">\n";

   if( binary_ )
   {
      Base64Writer base64;
      for( auto vertex = vc.begin(); vertex != vc.end(); ++vertex )
         base64 << numeric_cast<float>( std::get<0>( *vertex ) ) << numeric_cast<float>( std::get<1>( *vertex ) )
                << numeric_cast<float>( std::get<2>( *vertex ) );
      ofs << "     "; base64.toStream( ofs );
   }
   else for( auto vertex = vc.begin(); vertex != vc.end(); ++vertex )
      ofs << "     " << numeric_cast<float>( std::get<0>( *vertex ) ) << " " << numeric_cast<float>( std::get<1>( *vertex ) )
          << " " << numeric_cast<float>( std::get<2>( *vertex ) ) << "\n";

   ofs << "    </DataArray>\n"
       << "   </Points>\n"
       << "   <Cells>\n"
       << "    <DataArray type=\"" << vtk::typeToString< Index >() << "\" Name=\"connectivity\" format=\"" << format_ << "\">\n";

   if( binary_ )
   {
      Base64Writer base64;
      for( uint_t i = 0; i != ci.size(); i += 8 )
         base64 << ci[ i ] << ci[ i + 1 ] << ci[ i + 2 ] << ci[ i + 3 ] << ci[ i + 4 ] << ci[ i + 5 ] << ci[ i + 6 ] << ci[ i + 7 ];
      ofs << "     "; base64.toStream( ofs );
   }
   else for( uint_t i = 0; i != ci.size(); i += 8 )
      ofs << "     " << ci[ i ] << " " << ci[ i + 1 ] << " " << ci[ i + 2 ] << " " << ci[ i + 3 ] << " "
          << ci[ i + 4 ] << " " << ci[ i + 5 ] << " " << ci[ i + 6 ] << " " << ci[ i + 7 ] << "\n";

   ofs << "    </DataArray>\n"
       << "    <DataArray type=\"" << vtk::typeToString< Index >() << "\" Name=\"offsets\" format=\"" << format_ << "\">\n";

   if( binary_ )
   {
      Base64Writer base64;
      for( uint_t i = 0; i != ci.size(); i += 8 )
         base64 << numeric_cast<Index>( i + uint_c( 8 ) );
      ofs << "     "; base64.toStream( ofs );
   }
   else for( uint_t i = 0; i != ci.size(); i += 8 )
      ofs << "     " << numeric_cast<Index>( i + uint_c( 8 ) ) << "\n";

   ofs << "    </DataArray>\n"
       << "    <DataArray type=\"" << vtk::typeToString< uint8_t >() << "\" Name=\"types\" format=\"" << format_ << "\">\n";

   if( binary_ )
   {
      Base64Writer base64;
      for( uint_t i = 0; i != ci.size(); i += 8 )
         base64 << uint8_c( 11 );
      ofs << "     "; base64.toStream( ofs );
   }
   else for( uint_t i = 0; i != ci.size(); i += 8 )
      ofs << "     " << "11" << "\n";

   ofs << "    </DataArray>\n"
       << "   </Cells>\n";
}



uint8_t VTKOutput::ghostLayerNr( const IBlock& block, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ );

   const cell_idx_t xMax = cell_idx_c( blockStorage_->getNumberOfXCells(block) ) - cell_idx_c(1);
   const cell_idx_t yMax = cell_idx_c( blockStorage_->getNumberOfYCells(block) ) - cell_idx_c(1);
   const cell_idx_t zMax = cell_idx_c( blockStorage_->getNumberOfZCells(block) ) - cell_idx_c(1);

   cell_idx_t xDistance = 0;
   if( x < 0 ) xDistance = -x;
   else if( x > xMax ) xDistance = x - xMax;
   WALBERLA_ASSERT_GREATER_EQUAL( xDistance, 0 );

   cell_idx_t yDistance = 0;
   if( y < 0 ) yDistance = -y;
   else if( y > yMax ) yDistance = y - yMax;
   WALBERLA_ASSERT_GREATER_EQUAL( xDistance, 0 );

   cell_idx_t zDistance = 0;
   if( z < 0 ) zDistance = -z;
   else if( z > zMax ) zDistance = z - zMax;
   WALBERLA_ASSERT_GREATER_EQUAL( xDistance, 0 );

   return uint8_c( std::max( xDistance, std::max( yDistance, zDistance )) );
}



std::vector< VTKOutput::SamplingCell > VTKOutput::getSamplingCells( const IBlock& block, const CellVector& cells ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ );

   std::vector< SamplingCell > samplingCells;
   CellSet cellSet( cells );

   const AABB& domainBB  = blockStorage_->getDomain();
   const AABB& blockAABB = block.getAABB();

   const real_t xEnd = std::floor( ( blockAABB.xMax() - domainBB.xMin() ) / samplingDx_ );
   const real_t yEnd = std::floor( ( blockAABB.yMax() - domainBB.yMin() ) / samplingDy_ );
   const real_t zEnd = std::floor( ( blockAABB.zMax() - domainBB.zMin() ) / samplingDz_ );

   for( real_t z = std::floor( ( blockAABB.zMin() - domainBB.zMin() ) / samplingDz_ ); z <= zEnd; ++z )
   {
      real_t zMin = domainBB.zMin() +   z   * samplingDz_;
      real_t zMax = domainBB.zMin() + (z+1) * samplingDz_;

      for( real_t y = std::floor( ( blockAABB.yMin() - domainBB.yMin() ) / samplingDy_ ); y <= yEnd; ++y )
      {
         real_t yMin = domainBB.yMin() +   y   * samplingDy_;
         real_t yMax = domainBB.yMin() + (y+1) * samplingDy_;

         for( real_t x = std::floor( ( blockAABB.xMin() - domainBB.xMin() ) / samplingDx_ ); x <= xEnd; ++x )
         {
            real_t xMin = domainBB.xMin() +   x   * samplingDx_;
            real_t xMax = domainBB.xMin() + (x+1) * samplingDx_;

            SamplingCell cell;
            cell.aabb_.init( xMin, yMin, zMin, xMax, yMax, zMax );
            cell.globalX_ = ( xMin + xMax ) / real_c(2);
            cell.globalY_ = ( yMin + yMax ) / real_c(2);
            cell.globalZ_ = ( zMin + zMax ) / real_c(2);

            blockStorage_->getBlockLocalCell( cell.localCell_, block, cell.globalX_, cell.globalY_, cell.globalZ_ );

            if( cellSet.contains( cell.localCell_ ) )
            {
               cell.coordinates_.x() = cell_idx_c( x );
               cell.coordinates_.y() = cell_idx_c( y );
               cell.coordinates_.z() = cell_idx_c( z );

               AABB localCellAABB;
               blockStorage_->getBlockLocalCellAABB( block, cell.localCell_, localCellAABB );

               cell.localCellX_ = real_c( cell.localCell_.x() ) +
                                  ( ( cell.globalX_ - localCellAABB.xMin() ) / ( localCellAABB.xMax() - localCellAABB.xMin() ) );
               cell.localCellY_ = real_c( cell.localCell_.y() ) +
                                  ( ( cell.globalY_ - localCellAABB.yMin() ) / ( localCellAABB.yMax() - localCellAABB.yMin() ) );
               cell.localCellZ_ = real_c( cell.localCell_.z() ) +
                                  ( ( cell.globalZ_ - localCellAABB.zMin() ) / ( localCellAABB.zMax() - localCellAABB.zMin() ) );

               samplingCells.push_back( cell );
            }
         }
      }
   }

   return samplingCells;
}



void VTKOutput::writeCellData( std::ostream& ofs, const IBlock& block, const CellVector& cells ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ );

   for( auto writer = cellDataWriter_.begin(); writer != cellDataWriter_.end(); ++writer )
   {
      (*writer)->configure( block, *blockStorage_ );

      ofs << "    <DataArray type=\"" << (*writer)->typeString() << "\" Name=\"" << (*writer)->identifier()
                                      << "\" NumberOfComponents=\"" << (*writer)->fSize() << "\" format=\"" << format_ << "\">\n";

      if( binary_ )
      {
         Base64Writer base64;
         for( auto cell = cells.begin(); cell != cells.end(); ++cell )
            for( uint_t f = 0; f != (*writer)->fSize(); ++f )
               (*writer)->push( base64, cell->x(), cell->y(), cell->z(), cell_idx_c(f) );
         ofs << "     "; base64.toStream( ofs );
      }
      else
      {
         for( auto cell = cells.begin(); cell != cells.end(); ++cell ) {
            ofs << "     ";
            for( uint_t f = 0; f != (*writer)->fSize(); ++f )
            {
               (*writer)->push( ofs, cell->x(), cell->y(), cell->z(), cell_idx_c(f) );
               ofs << ( ( f == (*writer)->fSize() - 1 ) ? "\n" : " " );
            }
         }
      }

      ofs << "    </DataArray>\n";
   }
}


void VTKOutput::writeCellData( std::ostream& ofs, const Set<SUID>& requiredStates, const Set<SUID>& incompatibleStates ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ );

   for( auto writer = cellDataWriter_.begin(); writer != cellDataWriter_.end(); ++writer )
   {
      ofs << "    <DataArray type=\"" << (*writer)->typeString() << "\" Name=\"" << (*writer)->identifier()
          << "\" NumberOfComponents=\"" << (*writer)->fSize() << "\" format=\"" << format_ << "\">\n";

      for( auto block = blockStorage_->begin(); block != blockStorage_->end(); ++block )
      {
         if (!selectable::isSetSelected(uid::globalState() + block->getState(), requiredStates, incompatibleStates))
            continue;

         CellVector cells; // cells to be written to file
         computeVTUCells(*block, cells);
         (*writer)->configure( *block, *blockStorage_ );

         if( binary_ )
         {
            Base64Writer base64;
            for( auto cell = cells.begin(); cell != cells.end(); ++cell )
               for( uint_t f = 0; f != (*writer)->fSize(); ++f )
                  (*writer)->push( base64, cell->x(), cell->y(), cell->z(), cell_idx_c(f) );
            ofs << "     "; base64.toStream( ofs );
         }
         else
         {
            for( auto cell = cells.begin(); cell != cells.end(); ++cell ) {
               ofs << "     ";
               for( uint_t f = 0; f != (*writer)->fSize(); ++f )
               {
                  (*writer)->push( ofs, cell->x(), cell->y(), cell->z(), cell_idx_c(f) );
                  ofs << ( ( f == (*writer)->fSize() - 1 ) ? "\n" : " " );
               }
            }
         }
      }
      ofs << "    </DataArray>\n";
   }
}


void VTKOutput::writeCellData( std::ostream& ofs, const IBlock& block, const std::vector< SamplingCell >& cells ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ );

   for( auto writer = cellDataWriter_.begin(); writer != cellDataWriter_.end(); ++writer )
   {
      (*writer)->configure( block, *blockStorage_ );

      ofs << "    <DataArray type=\"" << (*writer)->typeString() << "\" Name=\"" << (*writer)->identifier()
                                      << "\" NumberOfComponents=\"" << (*writer)->fSize() << "\" format=\"" << format_ << "\">\n";

      if( binary_ )
      {
         Base64Writer base64;
         for( auto cell = cells.begin(); cell != cells.end(); ++cell )
            for( uint_t f = 0; f != (*writer)->fSize(); ++f )
               (*writer)->push( base64, cell->localCell_.x(), cell->localCell_.y(), cell->localCell_.z(), cell_idx_c(f),
                                        cell->localCellX_,    cell->localCellY_,    cell->localCellZ_,
                                        cell->globalX_   ,    cell->globalY_,       cell->globalZ_,
                                        samplingDx_,          samplingDy_,          samplingDz_ );
         ofs << "     "; base64.toStream( ofs );
      }
      else
      {
         for( auto cell = cells.begin(); cell != cells.end(); ++cell ) {
            ofs << "     ";
            for( uint_t f = 0; f != (*writer)->fSize(); ++f )
            {
               (*writer)->push( ofs, cell->localCell_.x(), cell->localCell_.y(), cell->localCell_.z(), cell_idx_c(f),
                                     cell->localCellX_,    cell->localCellY_,    cell->localCellZ_,
                                     cell->globalX_   ,    cell->globalY_,       cell->globalZ_,
                                     samplingDx_,          samplingDy_,          samplingDz_ );
               ofs << ( ( f == (*writer)->fSize() - 1 ) ? "\n" : " " );
            }
         }
      }

      ofs << "    </DataArray>\n";
   }
}



void VTKOutput::writeCollectors( const bool barrier )
{
   if( barrier )
      WALBERLA_MPI_WORLD_BARRIER();

   WALBERLA_NON_ROOT_SECTION() { return; }

   WALBERLA_ASSERT_EQUAL( MPIManager::instance()->worldRank(), 0 );

   if(!amrFileFormat_)
      writePVD();




   for( auto collector = collectorsToWrite_.begin(); collector != collectorsToWrite_.end(); ++collector )
   {
      if( uniformGrid_ )
      {
         if( samplingDx_ <= real_c(0) || samplingDy_ <= real_c(0) || samplingDz_ <= real_c(0) )
            writePVTI( *collector );
         else
            writePVTI_sampled( *collector );
      }
      else if (amrFileFormat_)
      {
         writeVTHBSeries();
         writeVTHB( *collector );
      }
      else
      {
         writePVTU( *collector ); // also applies for outputDomainDecomposition_ == true and pointDataSource_ != NULL
                                  // and polylineDataSource_ != NULL (uniformGrid_ will be false)
      }
   }

   collectorsToWrite_.clear();
}



void VTKOutput::writePVD()
{
   if ( !useMPIIO_  && collectorsToWrite_.empty() )
      return;

   std::string file( baseFolder_ + "/" + identifier_ + ".pvd" );
   std::fstream ofs( file.c_str() );

   if( !ofs ) // failed because file does not yet exist
   {
      ofs.open( file.c_str(), std::ios::out );

      ofs << "<?xml version=\"1.0\"?>\n"
          << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
          << " <Collection>\n";
   }
   else if( pvdEnd_ == std::streampos(-2) )
   {
      for( std::string line; std::getline(ofs, line); )
      {
         if( line.find("</Collection>") != std::string::npos )
         {
            WALBERLA_ASSERT_GREATER( ofs.tellg(), 0 );
            pvdEnd_ = ofs.tellg();
            pvdEnd_ -= int_c(line.size()) + 1;
            break;
         }
      }
      WALBERLA_ASSERT_GREATER( pvdEnd_, 0 );
      ofs.seekp(pvdEnd_);
   }
   else
   {
      ofs.seekp(pvdEnd_);
   }
   WALBERLA_ASSERT_GREATER(ofs.tellp(), 0);

   std::string ending;
   if( useMPIIO_ )
   {
      ending = ".vtu";
      if( uniformGrid_ )
         ending = ".vti";
   }
   else
   {
      ending = amrFileFormat_? ".vthb" : ".pvtu";
      if( uniformGrid_ )
         ending = ".pvti";
   }

   for( auto collector = allCollectors_.begin(); collector != allCollectors_.end(); ++collector )
   {
      std::ostringstream collection;
      collection << identifier_ << "/" << executionFolder_ << "_" << *collector << ending;
      ofs << "  <DataSet timestep=\"" << *collector << "\" file=\"" << collection.str() << "\"/>\n";
   }
   allCollectors_.clear();

   pvdEnd_ = ofs.tellp();
   WALBERLA_ASSERT_GREATER( pvdEnd_, 0 );
   ofs << " </Collection>\n"
       << "</VTKFile>\n";

   ofs.close();
}


void VTKOutput::writeVTHBSeries()
{
   if ( !useMPIIO_  && collectorsToWrite_.empty() )
      return;

   std::string file( baseFolder_ + "/" + identifier_ + ".vthb.series" );
   std::fstream ofs( file.c_str() );

   if( !ofs ) // failed because file does not yet exist
   {
      ofs.open( file.c_str(), std::ios::out );

      ofs << "{\n"
          << "   \"file-series-version\" : \"1.0\",\n"
          << "   \"files\" : [\n";
   }
   else if( pvdEnd_ == std::streampos(-2) )
   {
      for( std::string line; std::getline(ofs, line); )
      {
         if( line.find("]") != std::string::npos )
         {
            WALBERLA_ASSERT_GREATER( ofs.tellg(), 0 );
            pvdEnd_ = ofs.tellg();
            pvdEnd_ -= int_c(line.size()) + 1;
            break;
         }
      }
      WALBERLA_ASSERT_GREATER( pvdEnd_, 0 );
      ofs.seekp(pvdEnd_);
   }
   else
   {
      ofs.seekp(pvdEnd_);
   }
   WALBERLA_ASSERT_GREATER(ofs.tellp(), 0);

   std::string ending = ".vthb";

   for( auto collector = allCollectors_.begin(); collector != allCollectors_.end(); ++collector )
   {
      std::ostringstream collection;
      collection << identifier_ << "/" << executionFolder_ << "_" << *collector << ending;
      ofs << "      { \"name\" : \"" << collection.str() << "\", \"time\":" << *collector << " },\n";
   }
   allCollectors_.clear();

   pvdEnd_ = ofs.tellp();
   WALBERLA_ASSERT_GREATER( pvdEnd_, 0 );
   ofs << "   ]\n"
       << "}\n";

   ofs.close();
}



void VTKOutput::writePVTI( const uint_t collector ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ );

   std::ostringstream collection;
   collection << baseFolder_ << "/" << identifier_ << "/" << executionFolder_ << "_" << collector << ".pvti";
   std::ofstream ofs( collection.str().c_str() );

   const CellInterval& cellBB = blockStorage_->getDomainCellBB();
   const AABB&         domain = blockStorage_->getDomain();

   ofs << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
       << " <PImageData WholeExtent=\"" << cellBB.xMin() << " " << (cellBB.xMax()+1) << " "
                                        << cellBB.yMin() << " " << (cellBB.yMax()+1) << " "
                                        << cellBB.zMin() << " " << (cellBB.zMax()+1) << "\""
                        << " Origin=\"" << domain.xMin() << " " << domain.yMin() << " " << domain.zMin() << "\""
                       << " Spacing=\"" << blockStorage_->dx() << " " << blockStorage_->dy() << " " << blockStorage_->dz() << "\">\n"
       << "  <PCellData>\n";

   writePCellData( ofs );

   ofs << "  </PCellData>\n";

   std::vector< filesystem::path > files;
   getFilenames( files, collector );

   for( auto file = files.begin(); file != files.end(); ++file )
   {
      std::ifstream ifs( file->string().c_str() );

      std::string piece;
      for( uint_t i = 0; i != 4; ++i )
         std::getline( ifs, piece );
      ifs.close();

      piece.erase( piece.length()-1, 1 );

      ofs << piece << " Source=\"" << executionFolder_ << "_" << collector << "/" << file->filename().string() << "\"/>\n";
   }

   ofs << " </PImageData>\n"
       << "</VTKFile>\n";

   ofs.close();
}


void VTKOutput::writeVTHB( const uint_t collector ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ );

   std::ostringstream collection;
   collection << baseFolder_ << "/" << identifier_ << "/" << executionFolder_ << "_" << collector << ".vthb";
   std::ofstream ofs( collection.str().c_str() );

   ofs << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"vtkNonOverlappingAMR\" version=\"1.1\" byte_order=\"" << endianness_ << "\"" << " header_type=\"UInt32\" compressor=\"vtkZLibDataCompressor\">\n"
       << " <vtkNonOverlappingAMR>" << "\n";

   std::vector< std::vector< filesystem::path >> files;
   uint_t levels = blockStorage_->getNumberOfLevels();
   files.resize(levels);
   getFilenamesSortedByLevel( files, collector );

   for( uint_t level = 0; level < files.size(); level++){
      ofs << "  <Block level=\"" << level << "\">\n";
      uint index = 0;
      for( auto file = files[level].begin(); file != files[level].end(); ++file ){
         ofs << "   <DataSet index=\"" << index << "\" file=\"" << executionFolder_ << "_" << collector << "/" << file->filename().string() << "\"/>\n";
         index++;
      }
      ofs << "  </Block>\n";
   }

   ofs << " </vtkNonOverlappingAMR>\n"
       << "</VTKFile>\n";

   ofs.close();
}

void VTKOutput::writePVTI_sampled( const uint_t collector ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ );
   WALBERLA_ASSERT_GREATER( samplingDx_, real_t(0) );
   WALBERLA_ASSERT_GREATER( samplingDy_, real_t(0) );
   WALBERLA_ASSERT_GREATER( samplingDz_, real_t(0) );

   std::ostringstream collection;
   collection << baseFolder_ << "/" << identifier_ << "/" << executionFolder_ << "_" << collector << ".pvti";
   std::ofstream ofs( collection.str().c_str() );

   const AABB&         domain = blockStorage_->getDomain();
   const CellInterval  cellBB = getSampledCellInterval( domain );

   ofs << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
       << " <PImageData WholeExtent=\"" << cellBB.xMin() << " " << (cellBB.xMax()+1) << " "
                                        << cellBB.yMin() << " " << (cellBB.yMax()+1) << " "
                                        << cellBB.zMin() << " " << (cellBB.zMax()+1) << "\""
       << " Origin=\"" << domain.xMin() << " " << domain.yMin() << " " << domain.zMin() << "\""
       << " Spacing=\"" << samplingDx_ << " " << samplingDy_ << " " << samplingDz_ << "\">\n"
       << "  <PCellData>\n";

   writePCellData( ofs );

   ofs << "  </PCellData>\n";

   std::vector< filesystem::path > files;
   getFilenames( files, collector );

   for( auto file = files.begin(); file != files.end(); ++file )
   {
      std::ifstream ifs( file->string().c_str() );

      std::string piece;
      for( uint_t i = 0; i != 4; ++i )
         std::getline( ifs, piece );
      ifs.close();

      piece.erase( piece.length()-1, 1 );

      ofs << piece << " Source=\"" << executionFolder_ << "_" << collector << "/" << file->filename().string() << "\"/>\n";
   }

   ofs << " </PImageData>\n"
      << "</VTKFile>\n";

   ofs.close();
}

bool VTKOutput::writeCombinedVTI( std::string localPart, const uint_t collector ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ );

   const MPI_Comm comm    = MPIManager::instance()->comm();
   const int rank         = MPIManager::instance()->rank();
   const int numProcesses = MPIManager::instance()->numProcesses();

   const bool noData = mpi::allReduce( localPart.empty(), mpi::LOGICAL_AND, comm );
   if( noData )
      return false;

   std::ostringstream collection;
   collection << baseFolder_ << "/" << identifier_ << "/" << executionFolder_ << "_" << collector << ".vti";

   const CellInterval& cellBB = blockStorage_->getDomainCellBB();
   const AABB&         domain = blockStorage_->getDomain();

   if( rank == 0 )
   {
      std::ostringstream header;

      header << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
         << " <ImageData WholeExtent=\"" << cellBB.xMin() << " " << ( cellBB.xMax() + 1 ) << " "
         << cellBB.yMin() << " " << ( cellBB.yMax() + 1 ) << " "
         << cellBB.zMin() << " " << ( cellBB.zMax() + 1 ) << "\""
         << " Origin=\"" << domain.xMin() << " " << domain.yMin() << " " << domain.zMin() << "\""
         << " Spacing=\"" << blockStorage_->dx() << " " << blockStorage_->dy() << " " << blockStorage_->dz()
         << "\">\n\n";

      localPart.insert( 0, header.str() );
   }

   localPart.append( "\n\n" );

   if( rank == numProcesses - 1 )
   {
      localPart.append( " </ImageData>\n</VTKFile>\n" );
   }

   mpi::writeMPITextFile( collection.str(), localPart, comm );

   return true;
}

bool VTKOutput::writeCombinedVTI_sampled( std::string localPart, const uint_t collector ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ );
   WALBERLA_ASSERT_GREATER( samplingDx_, real_t( 0 ) );
   WALBERLA_ASSERT_GREATER( samplingDy_, real_t( 0 ) );
   WALBERLA_ASSERT_GREATER( samplingDz_, real_t( 0 ) );

   const MPI_Comm comm    = MPIManager::instance()->comm();
   const int rank         = MPIManager::instance()->rank();
   const int numProcesses = MPIManager::instance()->numProcesses();

   const bool noData = mpi::allReduce( localPart.empty(), mpi::LOGICAL_AND, comm );
   if( noData )
      return false;

   std::ostringstream collection;
   collection << baseFolder_ << "/" << identifier_ << "/" << executionFolder_ << "_" << collector << ".vti";

   if( rank == 0 )
   {
      std::ostringstream header;

      const AABB&         domain = blockStorage_->getDomain();
      const CellInterval  cellBB = getSampledCellInterval( domain );

      header << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
         << " <ImageData WholeExtent=\"" << cellBB.xMin() << " " << ( cellBB.xMax() + 1 ) << " "
         << cellBB.yMin() << " " << ( cellBB.yMax() + 1 ) << " "
         << cellBB.zMin() << " " << ( cellBB.zMax() + 1 ) << "\""
         << " Origin=\"" << domain.xMin() << " " << domain.yMin() << " " << domain.zMin() << "\""
         << " Spacing=\"" << samplingDx_ << " " << samplingDy_ << " " << samplingDz_ << "\">\n\n";

      localPart.insert( 0, header.str() );
   }

   localPart.append( "\n\n" );

   if( rank == numProcesses - 1 )
   {
      localPart.append( " </ImageData>\n</VTKFile>\n" );
   }

   mpi::writeMPITextFile( collection.str(), localPart, comm );

   return true;
}






void VTKOutput::writePVTU( const uint_t collector ) const
{
   std::ostringstream collection;
   collection << baseFolder_ << "/" << identifier_ << "/" << executionFolder_ << "_" << collector << ".pvtu";
   std::ofstream ofs( collection.str().c_str() );

   ofs << "<?xml version=\"1.0\"?>\n"
       << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
       << " <PUnstructuredGrid GhostLevel=\"" << ghostLayers_ << "\">\n"
       << "  <PPoints>\n"
       << "   <PDataArray type=\"" << vtk::typeToString< float >() << "\" NumberOfComponents=\"3\" format=\"" << format_ << "\"/>\n"
       << "  </PPoints>\n"
       << "  <PPointData>\n";

   writePPointData( ofs );

   ofs << "  </PPointData>\n"
       << "  <PCellData>\n";

   if( ghostLayers_ > 0 )
      ofs << "   <PDataArray type=\"" << vtk::typeToString< uint8_t >()
                                      << "\" Name=\"vtkGhostLevels\" NumberOfComponents=\"1\" format=\"" << format_ << "\"/>\n";
   writePCellData( ofs );

   ofs << "  </PCellData>\n";

   std::vector< filesystem::path > files;
   getFilenames( files, collector );

   for( auto file = files.begin(); file != files.end(); ++file )
      ofs << "  <Piece Source=\"" << executionFolder_ << "_" << collector << "/" << file->filename().string() << "\"/>\n";

   ofs << " </PUnstructuredGrid>\n"
       << "</VTKFile>\n";

   ofs.close();
}




bool VTKOutput::writeCombinedVTU( std::string localPart, const uint_t collector ) const
{
   const MPI_Comm comm    = MPIManager::instance()->comm();
   const int rank         = MPIManager::instance()->rank();
   const int numProcesses = MPIManager::instance()->numProcesses();

   const bool noData = mpi::allReduce( localPart.empty(), mpi::LOGICAL_AND, comm );
   if( noData )
      return false;

   std::ostringstream collection;
   collection << baseFolder_ << "/" << identifier_ << "/" << executionFolder_ << "_" << collector << ".vtu";

   if( rank == 0 )
   {
      std::ostringstream header;
      header << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << endianness_ << "\">\n"
         << " <UnstructuredGrid GhostLevel=\"" << ghostLayers_ << "\">\n\n";

      localPart.insert( 0, header.str() );
   }

   localPart.append( "\n\n" );

   if( rank == numProcesses - 1 )
   {
      localPart.append( " </UnstructuredGrid>\n</VTKFile>\n" );
   }

   mpi::writeMPITextFile( collection.str(), localPart, comm );

   return true;
}



void VTKOutput::getFilenames( std::vector< filesystem::path >& files, const uint_t collector ) const
{
   std::ostringstream path;
   path << baseFolder_ << "/" << identifier_ << "/" << executionFolder_ << "_" << collector;
   filesystem::path directory( path.str() );

   WALBERLA_ASSERT( filesystem::exists( directory ) );

   for( filesystem::directory_iterator file( directory ); file != filesystem::directory_iterator(); ++file )
   {
      WALBERLA_ASSERT( filesystem::is_regular_file( file->path() ) && !filesystem::is_directory( file->path() ) );
      files.push_back( file->path() );
   }
}


void VTKOutput::getFilenamesSortedByLevel( std::vector< std::vector< filesystem::path >>& blocks, const uint_t collector ) const
{
   std::ostringstream path;
   path << baseFolder_ << "/" << identifier_ << "/" << executionFolder_ << "_" << collector;
   filesystem::path directory( path.str() );

   WALBERLA_ASSERT( filesystem::exists( directory ) );

   for( filesystem::directory_iterator file( directory ); file != filesystem::directory_iterator(); ++file )
   {
      std::string fileName = file->path().string();
      WALBERLA_ASSERT( filesystem::is_regular_file( file->path() ) && !filesystem::is_directory( file->path() ) );

      std::size_t pos1 = fileName.find("level[");
      WALBERLA_ASSERT_UNEQUAL(pos1, std::string::npos, "file names of the block data must contain the block level for AMR data files")
      std::size_t pos2 = fileName.find("].vti");
      WALBERLA_ASSERT_UNEQUAL(pos2, std::string::npos, "files must be in vti format")
      std::size_t len = pos2 - (pos1 + 6);
      uint_t level = uint_c(std::stoi(fileName.substr(pos1 + 6, len)));
      WALBERLA_ASSERT_LESS(level, blocks.size())
      blocks[level].push_back(file->path());
   }
}



void VTKOutput::writePPointData( std::ofstream& ofs ) const
{
   WALBERLA_ASSERT( !( pointDataSource_ && polylineDataSource_ ) );

   if( pointDataSource_ )
   {
      auto attributes = pointDataSource_->getAttributes();
      for( auto dataArray = attributes.begin(); dataArray != attributes.end(); ++dataArray )
         ofs << "   <PDataArray type=\"" << dataArray->type << "\" Name=\"" << dataArray->name
                                         << "\" NumberOfComponents=\"" << dataArray->components << "\" format=\"" << format_ << "\"/>\n";
   }
   else if( polylineDataSource_ )
   {
      auto attributes = polylineDataSource_->getAttributes();
      for( auto dataArray = attributes.begin(); dataArray != attributes.end(); ++dataArray )
         ofs << "   <PDataArray type=\"" << dataArray->type << "\" Name=\"" << dataArray->name
                                         << "\" NumberOfComponents=\"" << dataArray->components << "\" format=\"" << format_ << "\"/>\n";
   }
}



void VTKOutput::writePCellData( std::ofstream& ofs ) const
{
   if( outputDomainDecomposition_ )
   {
      ofs << "   <PDataArray type=\"" << vtk::typeToString< uint8_t >()
                                      << "\" Name=\"Level\" NumberOfComponents=\"1\" format=\"" << format_ << "\"/>\n"
          << "   <PDataArray type=\"" << vtk::typeToString< int >()
                                      << "\" Name=\"Process\" NumberOfComponents=\"1\" format=\"" << format_ << "\"/>\n";
   }

   for( auto writer = cellDataWriter_.begin(); writer != cellDataWriter_.end(); ++writer )
   {
      ofs << "   <PDataArray type=\"" << (*writer)->typeString() << "\" Name=\"" << (*writer)->identifier()
                                      << "\" NumberOfComponents=\"" << (*writer)->fSize() << "\" format=\"" << format_ << "\"/>\n";
   }
}

CellInterval VTKOutput::getSampledCellInterval( const AABB & aabb ) const
{
   WALBERLA_ASSERT_GREATER( samplingDx_, real_t(0) );
   WALBERLA_ASSERT_GREATER( samplingDy_, real_t(0) );
   WALBERLA_ASSERT_GREATER( samplingDz_, real_t(0) );

   const AABB& domain = blockStorage_->getDomain();

   const real_t dx2 = real_c(0.5) * samplingDx_;
   const real_t dy2 = real_c(0.5) * samplingDy_;
   const real_t dz2 = real_c(0.5) * samplingDz_;

   CellInterval result;
   result.xMin() = cell_idx_c( std::floor( ( aabb.xMin() + dx2 - domain.xMin() ) / samplingDx_ ) );
   result.yMin() = cell_idx_c( std::floor( ( aabb.yMin() + dy2 - domain.yMin() ) / samplingDy_ ) );
   result.zMin() = cell_idx_c( std::floor( ( aabb.zMin() + dz2 - domain.zMin() ) / samplingDz_ ) );

   result.xMax() = cell_idx_c( std::floor( ( aabb.xMax() - dx2 - domain.xMin() ) / samplingDx_ ) );
   result.yMax() = cell_idx_c( std::floor( ( aabb.yMax() - dy2 - domain.yMin() ) / samplingDy_ ) );
   result.zMax() = cell_idx_c( std::floor( ( aabb.zMax() - dz2 - domain.zMin() ) / samplingDz_ ) );

   return result;
}



} // namespace vtk
} // namespace walberla
