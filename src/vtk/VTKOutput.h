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
//! \ingroup vtk
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "AABBCellFilter.h"
#include "Base64Writer.h"
#include "BlockCellDataWriter.h"
#include "CellBBCellFilter.h"
#include "PointDataSource.h"
#include "PolylineDataSource.h"

#include "core/Abort.h"
#include "core/DataTypes.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "core/Filesystem.h"
#include <functional>
#include <tuple>

#include <fstream>
#include <string>
#include <vector>


namespace walberla {
namespace vtk {



class VTKOutput : public NonCopyable {

private:

   class VTKGEN : public uid::IndexGenerator< VTKGEN, size_t >{};
   typedef UID< VTKGEN > VTKUID;

   // types used during vertex-index mapping procedure when writing (P)VTU files
   typedef std::tuple< cell_idx_t, cell_idx_t, cell_idx_t > Vertex;
   typedef std::tuple< real_t,     real_t,     real_t >     VertexCoord;
   typedef int32_t Index;

   struct VertexCompare {
      bool operator()( const Vertex& lhs, const Vertex& rhs ) const
      {
         if( std::get<0>(lhs) < std::get<0>(rhs) ||
             ( std::get<0>(lhs) == std::get<0>(rhs) && std::get<1>(lhs) < std::get<1>(rhs) ) ||
             ( std::get<0>(lhs) == std::get<0>(rhs) && std::get<1>(lhs) == std::get<1>(rhs) && std::get<2>(lhs) < std::get<2>(rhs) ) )
            return true;
         return false;
      }
   };

   struct SamplingCell
   {
      Cell coordinates_; // global cell coordinates of this "sampling cell" in "sampling space"

      Cell localCell_; // block local cell that contains the center of this "sampling cell"
      real_t localCellX_;
      real_t localCellY_;
      real_t localCellZ_;

      AABB aabb_;  // AABB of this "sampling cell" in global coordinates
      real_t globalX_; // center of ...
      real_t globalY_; // ... this "sampling cell" ...
      real_t globalZ_; // ... in global coordinates
   };

public:

   class Write {
   public:
      Write( const shared_ptr<VTKOutput>& vtk, const bool immediatelyWriteCollectors = true,
             const int simultaneousIOOperations = 0,
             const Set<SUID>& requiredStates     = Set<SUID>::emptySet(),
             const Set<SUID>& incompatibleStates = Set<SUID>::emptySet() ) :
         vtk_( vtk ), immediatelyWriteCollectors_( immediatelyWriteCollectors ), simultaneousIOOperations_( simultaneousIOOperations ),
         requiredStates_( requiredStates ), incompatibleStates_( incompatibleStates ) {}
      void operator()() const { vtk_->write( immediatelyWriteCollectors_, simultaneousIOOperations_, requiredStates_, incompatibleStates_ ); }
      void execute() const { this->operator()(); }
      void forceWriteNextStep() { vtk_->forceWriteNextStep(); }
   private:
      const shared_ptr<VTKOutput> vtk_;
      const bool immediatelyWriteCollectors_;
      const int simultaneousIOOperations_;
      const Set<SUID> requiredStates_;
      const Set<SUID> incompatibleStates_;
   };

   class WriteCollectors {
   public:
      WriteCollectors(  const shared_ptr<VTKOutput>& vtk, const bool barrier ) : vtk_( vtk ), barrier_( barrier ) {}
      void operator()() const { vtk_->writeCollectors( barrier_ ); }
   private:
      const shared_ptr<VTKOutput> vtk_;
      const bool barrier_;
   };

   /// creates a VTKOutput object that is supposed to output the domain decomposition
   friend inline shared_ptr<VTKOutput> createVTKOutput_DomainDecomposition( const BlockStorage & bs, const std::string & identifier,
                                                                            const uint_t writeFrequency,
                                                                            const std::string & baseFolder, const std::string & executionFolder,
                                                                            const bool continuousNumbering, const bool binary, const bool littleEndian,
                                                                            const bool useMPIIO, const uint_t initialExecutionCount );

   /// creates a VTKOutput object that is supposed to output block data/cell data
   friend inline shared_ptr<VTKOutput> createVTKOutput_BlockData( const StructuredBlockStorage & sbs, const std::string & identifier,
                                                                  const uint_t writeFrequency, const uint_t ghostLayers, const bool forcePVTU,
                                                                  const std::string & baseFolder, const std::string & executionFolder,
                                                                  const bool continuousNumbering, const bool binary, const bool littleEndian,
                                                                  const bool useMPIIO, const uint_t initialExecutionCount );

   /// creates a VTKOutput object that is supposed to output arbitrary point data
   friend inline shared_ptr<VTKOutput> createVTKOutput_PointData( const shared_ptr< PointDataSource > pds, const std::string & identifier,
                                                                  const uint_t writeFrequency,
                                                                  const std::string & baseFolder, const std::string & executionFolder,
                                                                  const bool continuousNumbering, const bool binary, const bool littleEndian,
                                                                  const bool useMPIIO, const uint_t initialExecutionCount );

   /// creates a VTKOutput object that is supposed to output arbitrary polyline data
   friend inline shared_ptr<VTKOutput> createVTKOutput_PolylineData( const shared_ptr< PolylineDataSource > pds, const std::string & identifier,
                                                                     const uint_t writeFrequency,
                                                                     const std::string & baseFolder, const std::string & executionFolder,
                                                                     const bool continuousNumbering, const bool binary, const bool littleEndian,
                                                                     const bool useMPIIO, const uint_t initialExecutionCount );

   typedef std::function< void () > BeforeFunction;
   typedef std::function< void ( CellSet& filteredCells, const IBlock& block,
                                   const StructuredBlockStorage& storage, const uint_t ghostLayers ) >  CellFilter;

   ~VTKOutput();

   inline void setInitialWriteCallsToSkip( const uint_t initialWriteCallsToSkip );

   void addBeforeFunction( BeforeFunction f ) { beforeFunctions_.push_back( f ); }

   inline void addAABBInclusionFilter( const AABB & aabb );
   inline void addAABBExclusionFilter( const AABB & aabb );

   // 'filteredCells' contains block local cell coordinates of all cells that are written to file (if the cell is no excluded!)
   inline void addCellInclusionFilter( CellFilter f );
   // 'filteredCells' contains block local cell coordinates of all cells that must not be written to file
   inline void addCellExclusionFilter( CellFilter f );

   inline void addCellDataWriter( const shared_ptr< BlockCellDataWriterInterface > & writer );

   inline void setSamplingResolution( const real_t spacing );
   inline void setSamplingResolution( const real_t dx, const real_t dy, const real_t dz );

   void write( const bool immediatelyWriteCollectors = true,
               const int simultaneousIOOperations = 0,
               const Set<SUID>& requiredStates     = Set<SUID>::emptySet(),
               const Set<SUID>& incompatibleStates = Set<SUID>::emptySet() );

   void forceWrite( uint_t number,
                    const bool immediatelyWriteCollectors = true,
                    const int simultaneousIOOperations = 0,
                    const Set<SUID>& requiredStates     = Set<SUID>::emptySet(),
                    const Set<SUID>& incompatibleStates = Set<SUID>::emptySet() );

   void writeCollectors( const bool barrier );

   inline void forceWriteNextStep() { writeNextStep_ = true; }

private:

   VTKOutput();

   /// creates a VTKOutput object that is supposed to output the domain decomposition
   VTKOutput( const BlockStorage & sbs, const std::string & identifier, const uint_t writeFrequency,
              const std::string & baseFolder, const std::string & executionFolder,
              const bool continuousNumbering, const bool binary, const bool littleEndian, const bool useMPIIO,
              const uint_t initialExecutionCount = 0 );

   /// creates a VTKOutput object that is supposed to output block data/cell data
   VTKOutput( const StructuredBlockStorage & sbs, const std::string & identifier, const uint_t writeFrequency,
              const std::string & baseFolder, const std::string & executionFolder,
              const bool continuousNumbering, const bool binary, const bool littleEndian, const bool useMPIIO,
              const uint_t ghostLayers, const bool forcePVTU, const uint_t initialExecutionCount = 0 );

   /// creates a VTKOutput object that is supposed to output arbitrary point data
   VTKOutput( const shared_ptr< PointDataSource >& pds, const std::string & identifier, const uint_t writeFrequency,
              const std::string & baseFolder, const std::string & executionFolder,
              const bool continuousNumbering, const bool binary, const bool littleEndian, const bool useMPIIO,
              const uint_t initialExecutionCount = 0 );

   /// creates a VTKOutput object that is supposed to output arbitrary polyline data
   VTKOutput( const shared_ptr< PolylineDataSource >& pds, const std::string & identifier, const uint_t writeFrequency,
              const std::string & baseFolder, const std::string & executionFolder,
              const bool continuousNumbering, const bool binary, const bool littleEndian, const bool useMPIIO,
              const uint_t initialExecutionCount = 0 );

   void init( const std::string & identifier );

   void writeDomainDecomposition( const std::string& path, const Set<SUID>& requiredStates, const Set<SUID>& incompatibleStates ) const;
   void writeDomainDecompositionPieces( std::ostream& ofs, const Set<SUID>& requiredStates, const Set<SUID>& incompatibleStates ) const;

   void computeOutputPoints( std::vector<Vector3<real_t> > & points, std::vector<bool> & outputPoint,
                             uint_t & numberOfPoints ) const;
   void writePointDataPieceHelper( const std::vector<Vector3<real_t> > & points, const std::vector<bool> & outputPoint,
                                   const uint_t numberOfPoints, std::ostream & ofs ) const;
   void writePointData( const std::string& path );
   void writePointDataPieces(std::ostream& ofs);


   void computeOutputPolylines( std::vector< std::vector< Vector3< real_t > > > & lines,
      std::vector< std::vector< bool > > & outputPolylinePoint, std::vector< size_t > & polylineSize,
      uint_t & numberOfPolylines, uint_t & numberOfPolylinePoints ) const;

   void writePolylineDataPieceHelper( const std::vector< std::vector< Vector3< real_t > > > & lines,
      const std::vector< std::vector< bool > > & outputPolylinePoint, const std::vector< size_t > & polylineSize,
      const uint_t numberOfPolylines, const uint_t numberOfPolylinePoints, std::ostream & ofs ) const;

   void writePolylineData( const std::string& path );
   void writePolylineDataPieces( std::ostream& ofs );

   void computeVTUCells( const IBlock& block, CellVector & cellsOut ) const;

   void writeBlocks( const std::string& path, const Set<SUID>& requiredStates, const Set<SUID>& incompatibleStates );
   void writeBlockPieces( std::ostream & oss, const Set<SUID>& requiredStates, const Set<SUID>& incompatibleStates );

   void writeVTI( std::ostream& ofs, const IBlock& block ) const;
   void writeVTI_sampling( std::ostream& ofs, const IBlock& block ) const;

   void writeVTIPiece( std::ostream& ofs, const IBlock& block ) const;
   void writeVTIPiece_sampling( std::ostream& ofs, const IBlock& block ) const;

   void writeVTU( std::ostream& ofs, const IBlock& block, const CellVector& cells ) const;
   void writeVTU_sampling( std::ostream& ofs, const IBlock& block, const CellVector& cells ) const;

   void writeVTUPiece(std::ostream& ofs, const IBlock& block, const CellVector& cells) const;
   void writeVTUPiece_sampling(std::ostream& ofs, const IBlock& block, const CellVector& cells) const;

   void writeVTUHeader( std::ofstream& ofs, const uint_t numberOfCells, const std::vector< VertexCoord > & vc, const std::vector< Index > & ci ) const;
   void writeVTUHeaderPiece (std::ostream& ofs, const uint_t numberOfCells, const std::vector< VertexCoord > & vc, const std::vector< Index > & ci) const;

   uint8_t ghostLayerNr( const IBlock& block, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;

   std::vector< SamplingCell > getSamplingCells( const IBlock& block, const CellVector& cells ) const;

   void writeCellData( std::ostream& ofs, const IBlock& block, const CellVector& cells ) const;
   void writeCellData( std::ostream& ofs, const IBlock& block, const std::vector< SamplingCell >& cells ) const;

   void writePVD();

   void writePVTI( const uint_t collector ) const;
   void writePVTI_sampled( const uint_t collector ) const;
   void writePVTU( const uint_t collector ) const;

   bool writeCombinedVTI( std::string localPart, const uint_t collector ) const;
   bool writeCombinedVTI_sampled( std::string localPart, const uint_t collector ) const;
   bool writeCombinedVTU(std::string localPart, const uint_t collector) const;

   void getFilenames( std::vector< filesystem::path >& blocks, const uint_t collector ) const;
   void writePPointData( std::ofstream& ofs ) const;
   void writePCellData( std::ofstream& ofs ) const;

   CellInterval getSampledCellInterval( const AABB & aabb ) const;


   std::string identifier_;

   std::vector< BeforeFunction > beforeFunctions_;

   const BlockStorage * const unstructuredBlockStorage_;
   const StructuredBlockStorage * const blockStorage_;

   const shared_ptr< PointDataSource >    pointDataSource_;
   const shared_ptr< PolylineDataSource > polylineDataSource_;

   const std::string baseFolder_;
   const std::string executionFolder_;

   uint_t executionCounter_;
   uint_t initialWriteCallsToSkip_;
   const uint_t writeFrequency_;
   const bool continuousNumbering_;

   std::streampos pvdEnd_;
   std::vector< uint_t > allCollectors_;
   std::vector< uint_t > collectorsToWrite_;

   const bool binary_;
   const std::string format_;     // "binary" or "ascii"
   const std::string endianness_; // "LittleEndian" or "BigEndian"

   const bool useMPIIO_;

   const bool outputDomainDecomposition_; // if true, only the block structure (= the domain decomposition) is written to file

   real_t samplingDx_;
   real_t samplingDy_;
   real_t samplingDz_;

   const bool forcePVTU_;
         bool configured_;
         bool uniformGrid_;

   const uint_t ghostLayers_;

   std::vector< AABB >  aabbInclusionFilters_;
   std::vector< AABB >  aabbExclusionFilters_;

   std::vector< CellFilter >  cellInclusionFunctions_;
   std::vector< CellFilter >  cellExclusionFunctions_;

   std::vector< shared_ptr< BlockCellDataWriterInterface > >  cellDataWriter_;

   bool writeNextStep_;

}; // class VTKOutput



inline void VTKOutput::setInitialWriteCallsToSkip( const uint_t initialWriteCallsToSkip )
{
   if( executionCounter_ > 0 )
      WALBERLA_ABORT( "Setting the number of initial write calls to skip is only possible as long as \"write\" has not yet been called!" );

   initialWriteCallsToSkip_ = initialWriteCallsToSkip;
}



inline void VTKOutput::addAABBInclusionFilter( const AABB & aabb )
{
   if( outputDomainDecomposition_ )
      WALBERLA_ABORT( "You are trying to add an AABB inclusion filter to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is intended for outputting the domain decomposition, AABB filters won't have any effect." );
   if( configured_ )
      WALBERLA_ABORT( "You are trying to add an AABB inclusion filter to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is already configured.\nYou can only add filters as long as no output has been written." );

   if( pointDataSource_ || polylineDataSource_ )
      aabbInclusionFilters_.push_back( aabb );
   else
      cellInclusionFunctions_.push_back( AABBCellFilter(aabb) );
}



inline void VTKOutput::addAABBExclusionFilter( const AABB & aabb )
{
   if( outputDomainDecomposition_ )
      WALBERLA_ABORT( "You are trying to add an AABB exclusion filter to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is intended for outputting the domain decomposition, AABB filters won't have any effect." );
   if( configured_ )
      WALBERLA_ABORT( "You are trying to add an AABB exclusion filter to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is already configured.\nYou can only add filters as long as no output has been written." );

   if( pointDataSource_ || polylineDataSource_ )
      aabbExclusionFilters_.push_back( aabb );
   else
      cellExclusionFunctions_.push_back( AABBCellFilter(aabb) );
}



inline void VTKOutput::addCellInclusionFilter( CellFilter f )
{
   if( outputDomainDecomposition_ )
      WALBERLA_ABORT( "You are trying to add a cell inclusion filter to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is intended for outputting the domain decomposition, cell filters won't have any effect." );
   if( pointDataSource_ )
      WALBERLA_ABORT( "You are trying to add a cell inclusion filter to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is intended for outputting arbitrary point data, cell filters won't have any effect." );

   if( polylineDataSource_ )
      WALBERLA_ABORT( "You are trying to add a cell inclusion filter to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is intended for outputting arbitrary polyline data, cell filters won't have any effect." );

   if( configured_ )
      WALBERLA_ABORT( "You are trying to add a cell inclusion filter to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is already configured.\nYou can only add filters as long as no output has been written." );

   cellInclusionFunctions_.push_back(f);
}



inline void VTKOutput::addCellExclusionFilter( CellFilter f )
{
   if( outputDomainDecomposition_ )
      WALBERLA_ABORT( "You are trying to add a cell exclusion filter to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is intended for outputting the domain decomposition, cell filters won't have any effect." );
   if( pointDataSource_ )
      WALBERLA_ABORT( "You are trying to add a cell exclusion filter to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is intended for outputting arbitrary point data, cell filters won't have any effect." );

   if( polylineDataSource_ )
      WALBERLA_ABORT( "You are trying to add a cell exclusion filter to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is intended for outputting arbitrary polyline data, cell filters won't have any effect." );

   if( configured_ )
      WALBERLA_ABORT( "You are trying to add a cell exclusion filter to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is already configured.\nYou can only add filters as long as no output has been written." );

   cellExclusionFunctions_.push_back(f);
}



inline void VTKOutput::addCellDataWriter( const shared_ptr< BlockCellDataWriterInterface > & writer )
{
   if( outputDomainDecomposition_ )
      WALBERLA_ABORT( "You are trying to add a block cell data writer to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is intended for outputting the domain decomposition, cell data writers won't have any effect." );
   if( pointDataSource_ )
      WALBERLA_ABORT( "You are trying to add a block cell data writer to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is intended for outputting arbitrary point data, cell data writers won't have any effect." );

   if( polylineDataSource_ )
      WALBERLA_ABORT( "You are trying to add a block cell data writer to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is intended for outputting arbitrary polyline data, cell data writers won't have any effect." );

   if( configured_ )
      WALBERLA_ABORT( "You are trying to add a block cell data writer to VTKOutput \"" << identifier_ << "\", "
                      "but this VTKOutput is already configured.\nYou can only add cell data writers as long as no output has been written." );

   cellDataWriter_.push_back( writer );
}



inline void VTKOutput::setSamplingResolution( const real_t spacing )
{
   if( spacing > real_c(0) )
   {
      if( outputDomainDecomposition_ )
         WALBERLA_ABORT( "You are trying to set a sampling resolution of " << spacing << " to VTKOutput \"" << identifier_ << "\", "
                         "but this VTKOutput is intended for outputting the domain decomposition, setting a sampling resolution won't have any effect." );
      if( pointDataSource_ )
         WALBERLA_ABORT( "You are trying to set a sampling resolution of " << spacing << " to VTKOutput \"" << identifier_ << "\", "
                         "but this VTKOutput is intended for outputting arbitrary point data, setting a sampling resolution won't have any effect." );

      if( polylineDataSource_ )
         WALBERLA_ABORT( "You are trying to set a sampling resolution of " << spacing << " to VTKOutput \"" << identifier_ << "\", "
                         "but this VTKOutput is intended for outputting arbitrary polyline data, setting a sampling resolution won't have any effect." );

      if( configured_ )
         WALBERLA_ABORT( "You are trying to set a sampling resolution of " << spacing << " to VTKOutput \"" << identifier_ << "\", "
                         "but this VTKOutput is already configured.\nYou can only change the sampling resolution as long as no output has been written." );
      if( ghostLayers_ > 0 )
         WALBERLA_ABORT( "You are trying to set a sampling resolution of " << spacing << " to VTKOutput \"" << identifier_ << "\", "
                         "but this VTKOutput is configured to output " << ghostLayers_ << " ghost layer(s).\n"
                         "Writing ghost layers and using sampling at the same time is not supported!" );
   }

   samplingDx_ = spacing;
   samplingDy_ = spacing;
   samplingDz_ = spacing;
}



inline void VTKOutput::setSamplingResolution( const real_t dx, const real_t dy, const real_t dz )
{
   if( dx > real_c(0) && dy > real_c(0) && dz > real_c(0) )
   {
      if( outputDomainDecomposition_ )
         WALBERLA_ABORT( "You are trying to set a sampling resolution of ( " << dx << ", " << dy << ", " << dz << " ) "
                         "to VTKOutput \"" << identifier_ << "\", "
                         "but this VTKOutput is intended for outputting the domain decomposition, setting a sampling resolution won't have any effect." );
      if( pointDataSource_ )
         WALBERLA_ABORT( "You are trying to set a sampling resolution of ( " << dx << ", " << dy << ", " << dz << " ) "
                         "to VTKOutput \"" << identifier_ << "\", "
                         "but this VTKOutput is intended for outputting arbitrary point data, setting a sampling resolution won't have any effect." );

      if( polylineDataSource_ )
         WALBERLA_ABORT( "You are trying to set a sampling resolution of ( " << dx << ", " << dy << ", " << dz << " ) "
                         "to VTKOutput \"" << identifier_ << "\", "
                         "but this VTKOutput is intended for outputting arbitrary polyline data, setting a sampling resolution won't have any effect." );

      if( configured_ )
         WALBERLA_ABORT( "You are trying to set a sampling resolution of ( " << dx << ", " << dy << ", " << dz << " ) "
                         "to VTKOutput \"" << identifier_ << "\", "
                         "but this VTKOutput is already configured.\nYou can only change the sampling resolution as long as no output has been written." );
      if( ghostLayers_ > 0 )
         WALBERLA_ABORT( "You are trying to set a sampling resolution of ( " << dx << ", " << dy << ", " << dz << " ) "
                         "to VTKOutput \"" << identifier_ << "\", but this VTKOutput is configured to output " << ghostLayers_ << " ghost layer(s).\n"
                         "Writing ghost layers and using sampling at the same time is not supported!" );
   }

   samplingDx_ = dx;
   samplingDy_ = dy;
   samplingDz_ = dz;
}






////////////////////
// FREE FUNCTIONS //
////////////////////



inline uint_t determineWriteFrequency( const real_t dt_SI, const uint_t fps )
{
   return uint_c( real_c(1.0) / ( dt_SI * real_c(fps) ) + real_t(0.5) );
}



inline shared_ptr<VTKOutput> createVTKOutput_DomainDecomposition( const BlockStorage & bs,
                                                                  const std::string & identifier = std::string("domain_decomposition"),
                                                                  const uint_t writeFrequency = 1,
                                                                  const std::string & baseFolder = std::string("vtk_out"),
                                                                  const std::string & executionFolder = std::string("simulation_step"),
                                                                  const bool continuousNumbering = false, const bool binary = true,
                                                                  const bool littleEndian = true, const bool useMPIIO = true,
                                                                  const uint_t initialExecutionCount = 0 )
{
   return shared_ptr<VTKOutput>( new VTKOutput( bs, identifier, writeFrequency, baseFolder, executionFolder,
                                                continuousNumbering, binary, littleEndian, useMPIIO, initialExecutionCount ) );
}
inline shared_ptr<VTKOutput> createVTKOutput_DomainDecomposition( const StructuredBlockStorage & sbs,
                                                                  const std::string & identifier = std::string("domain_decomposition"),
                                                                  const uint_t writeFrequency = 1,
                                                                  const std::string & baseFolder = std::string("vtk_out"),
                                                                  const std::string & executionFolder = std::string("simulation_step"),
                                                                  const bool continuousNumbering = false, const bool binary = true,
                                                                  const bool littleEndian = true, const bool useMPIIO = true,
                                                                  const uint_t initialExecutionCount = 0 )
{
   return createVTKOutput_DomainDecomposition( sbs.getBlockStorage(), identifier, writeFrequency, baseFolder, executionFolder,
                                               continuousNumbering, binary, littleEndian, useMPIIO, initialExecutionCount );
}

inline shared_ptr<VTKOutput> createVTKOutput_DomainDecomposition( const shared_ptr< const BlockStorage > & sbs,
                                                                  const std::string & identifier = std::string("domain_decomposition"),
                                                                  const uint_t writeFrequency = 1,
                                                                  const std::string & baseFolder = std::string("vtk_out"),
                                                                  const std::string & executionFolder = std::string("simulation_step"),
                                                                  const bool continuousNumbering = false, const bool binary = true,
                                                                  const bool littleEndian = true, const bool useMPIIO = true,
                                                                  const uint_t initialExecutionCount = 0 )
{
   if( !sbs )
      WALBERLA_ABORT( "creating VTK output for domain decomposition failed (StructuredBlockStorage shared pointer is NULL)" );

   return createVTKOutput_DomainDecomposition( *sbs, identifier, writeFrequency, baseFolder, executionFolder,
                                               continuousNumbering, binary, littleEndian, useMPIIO, initialExecutionCount );
}
inline shared_ptr<VTKOutput> createVTKOutput_DomainDecomposition( const shared_ptr< const StructuredBlockStorage > & sbs,
                                                                  const std::string & identifier = std::string("domain_decomposition"),
                                                                  const uint_t writeFrequency = 1,
                                                                  const std::string & baseFolder = std::string("vtk_out"),
                                                                  const std::string & executionFolder = std::string("simulation_step"),
                                                                  const bool continuousNumbering = false, const bool binary = true,
                                                                  const bool littleEndian = true, const bool useMPIIO = true,
                                                                  const uint_t initialExecutionCount = 0 )
{
   if( !sbs )
      WALBERLA_ABORT( "creating VTK output for domain decomposition failed (StructuredBlockStorage shared pointer is NULL)" );

   return createVTKOutput_DomainDecomposition( sbs->getBlockStorage(), identifier, writeFrequency, baseFolder, executionFolder,
                                               continuousNumbering, binary, littleEndian, useMPIIO, initialExecutionCount );
}



inline shared_ptr<VTKOutput> createVTKOutput_BlockData( const StructuredBlockStorage & sbs,
                                                        const std::string & identifier = std::string("block_data"),
                                                        const uint_t writeFrequency = 1, const uint_t ghostLayers = 0, const bool forcePVTU = false,
                                                        const std::string & baseFolder = std::string("vtk_out"),
                                                        const std::string & executionFolder = std::string("simulation_step"),
                                                        const bool continuousNumbering = false, const bool binary = true,
                                                        const bool littleEndian = true, const bool useMPIIO = true,
                                                        const uint_t initialExecutionCount = 0 )
{
   return shared_ptr<VTKOutput>( new VTKOutput( sbs, identifier, writeFrequency, baseFolder, executionFolder,
                                                continuousNumbering, binary, littleEndian, useMPIIO, ghostLayers, forcePVTU, initialExecutionCount ) );
}

inline shared_ptr<VTKOutput> createVTKOutput_BlockData( const shared_ptr< const StructuredBlockStorage > & sbs,
                                                        const std::string & identifier = std::string("block_data"),
                                                        const uint_t writeFrequency = 1, const uint_t ghostLayers = 0, const bool forcePVTU = false,
                                                        const std::string & baseFolder = std::string("vtk_out"),
                                                        const std::string & executionFolder = std::string("simulation_step"),
                                                        const bool continuousNumbering = false, const bool binary = true,
                                                        const bool littleEndian = true, const bool useMPIIO = true,
                                                        const uint_t initialExecutionCount = 0 )
{
   if( !sbs )
      WALBERLA_ABORT( "creating VTK output for block data failed (StructuredBlockStorage shared pointer is NULL)" );

   return createVTKOutput_BlockData( *sbs, identifier, writeFrequency, ghostLayers, forcePVTU, baseFolder, executionFolder,
                                               continuousNumbering, binary, littleEndian, useMPIIO, initialExecutionCount );
}



inline shared_ptr<VTKOutput> createVTKOutput_PointData( const shared_ptr< PointDataSource > pds,
                                                        const std::string & identifier = std::string("point_data"),
                                                        const uint_t writeFrequency = 1,
                                                        const std::string & baseFolder = std::string("vtk_out"),
                                                        const std::string & executionFolder = std::string("simulation_step"),
                                                        const bool continuousNumbering = false, const bool binary = true,
                                                        const bool littleEndian = true, const bool useMPIIO = true,
                                                        const uint_t initialExecutionCount = 0 )
{
   return shared_ptr<VTKOutput>( new VTKOutput( pds, identifier, writeFrequency, baseFolder, executionFolder,
                                                continuousNumbering, binary, littleEndian, useMPIIO, initialExecutionCount ) );
}


inline shared_ptr<VTKOutput> createVTKOutput_PolylineData( const shared_ptr< PolylineDataSource > pds,
                                                           const std::string & identifier = std::string("point_data"),
                                                           const uint_t writeFrequency = 1,
                                                           const std::string & baseFolder = std::string("vtk_out"),
                                                           const std::string & executionFolder = std::string("simulation_step"),
                                                           const bool continuousNumbering = false, const bool binary = true,
                                                           const bool littleEndian = true, const bool useMPIIO = true,
                                                           const uint_t initialExecutionCount = 0 )
{
   return shared_ptr<VTKOutput>( new VTKOutput( pds, identifier, writeFrequency, baseFolder, executionFolder,
                                                continuousNumbering, binary, littleEndian, useMPIIO, initialExecutionCount ) );
}



inline void writeDomainDecomposition( const StructuredBlockStorage & sbs,
                                      const std::string & identifier = std::string("domain_decomposition"),
                                      const std::string & baseFolder = std::string("vtk_out"),
                                      const std::string & executionFolder = std::string("write_call"),
                                      const bool binary = true, const bool littleEndian = true,
                                      const int simultaneousIOOperations = 0,
                                      const Set<SUID>& requiredStates     = Set<SUID>::emptySet(),
                                      const Set<SUID>& incompatibleStates = Set<SUID>::emptySet(),
                                      bool useMPIIO = true )
{
   auto vtkOutput = createVTKOutput_DomainDecomposition( sbs, identifier, 1, baseFolder, executionFolder, true, binary, littleEndian, useMPIIO );
   vtkOutput->write( true, simultaneousIOOperations, requiredStates, incompatibleStates );
}
inline void writeDomainDecomposition( const BlockStorage & bs,
                                      const std::string & identifier = std::string("domain_decomposition"),
                                      const std::string & baseFolder = std::string("vtk_out"),
                                      const std::string & executionFolder = std::string("write_call"),
                                      const bool binary = true, const bool littleEndian = true,
                                      const int simultaneousIOOperations = 0,
                                      const Set<SUID>& requiredStates     = Set<SUID>::emptySet(),
                                      const Set<SUID>& incompatibleStates = Set<SUID>::emptySet(),
                                      bool useMPIIO = true )
{
   auto vtkOutput = createVTKOutput_DomainDecomposition( bs, identifier, 1, baseFolder, executionFolder, true, binary, littleEndian, useMPIIO );
   vtkOutput->write( true, simultaneousIOOperations, requiredStates, incompatibleStates );
}

inline void writeDomainDecomposition( const shared_ptr< const BlockStorage > & bs,
                                      const std::string & identifier = std::string("domain_decomposition"),
                                      const std::string & baseFolder = std::string("vtk_out"),
                                      const std::string & executionFolder = std::string("write_call"),
                                      const bool binary = true, const bool littleEndian = true,
                                      const int simultaneousIOOperations = 0,
                                      const Set<SUID>& requiredStates     = Set<SUID>::emptySet(),
                                      const Set<SUID>& incompatibleStates = Set<SUID>::emptySet(),
                                      bool useMPIIO = true )
{
   if( !bs )
      WALBERLA_ABORT( "Writing domain decomposition failed (StructuredBlockStorage shared pointer is NULL)" );

   writeDomainDecomposition( *bs, identifier, baseFolder, executionFolder, binary, littleEndian,
                             simultaneousIOOperations, requiredStates, incompatibleStates, useMPIIO );
}
inline void writeDomainDecomposition( const shared_ptr< const StructuredBlockStorage > & sbs,
                                      const std::string & identifier = std::string("domain_decomposition"),
                                      const std::string & baseFolder = std::string("vtk_out"),
                                      const std::string & executionFolder = std::string("write_call"),
                                      const bool binary = true, const bool littleEndian = true,
                                      const int simultaneousIOOperations = 0,
                                      const Set<SUID>& requiredStates     = Set<SUID>::emptySet(),
                                      const Set<SUID>& incompatibleStates = Set<SUID>::emptySet(),
                                      bool useMPIIO = true )
{
   if( !sbs )
      WALBERLA_ABORT( "Writing domain decomposition failed (StructuredBlockStorage shared pointer is NULL)" );

   writeDomainDecomposition( *sbs, identifier, baseFolder, executionFolder, binary, littleEndian,
                             simultaneousIOOperations, requiredStates, incompatibleStates, useMPIIO );
}



inline VTKOutput::Write writeFiles( const shared_ptr<VTKOutput> & vtk,
                                    const bool immediatelyWriteCollectors = true,
                                    const int simultaneousIOOperations = 0,
                                    const Set<SUID>& requiredStates     = Set<SUID>::emptySet(),
                                    const Set<SUID>& incompatibleStates = Set<SUID>::emptySet() )
{
   return VTKOutput::Write( vtk, immediatelyWriteCollectors, simultaneousIOOperations, requiredStates, incompatibleStates );
}



inline VTKOutput::WriteCollectors writeCollectorFiles( const shared_ptr<VTKOutput> & vtk, const bool barrier )
{
   return VTKOutput::WriteCollectors( vtk, barrier );
}



} // namespace vtk
} // namespace walberla


