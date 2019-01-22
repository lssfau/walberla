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
//! \file StructuredBlockStorage.h
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "BlockStorage.h"

#include "core/DataTypes.h"
#include "core/NonCopyable.h"
#include "core/cell/Cell.h"
#include "core/cell/CellInterval.h"
#include "core/debug/Debug.h"

#include <functional>


namespace walberla {
namespace domain_decomposition {



class StructuredBlockStorage;

template< typename T >
struct StructuredBlockDataCreator {

   StructuredBlockDataCreator(
         const std::function< T * ( IBlock * const block, StructuredBlockStorage * const storage ) > & function,
         const std::string & identifier            = std::string(),
         const Set<SUID> &   requiredSelectors     = Set<SUID>::emptySet(),
         const Set<SUID> &   incompatibleSelectors = Set<SUID>::emptySet() ) :
      function_( function ), identifier_( identifier ),
      requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors ) {}

   std::function< T * ( IBlock * const block, StructuredBlockStorage * const storage ) > function_;

   std::string identifier_;
   Set<SUID>   requiredSelectors_;
   Set<SUID>   incompatibleSelectors_;
};



//**********************************************************************************************************************
/*!
*   \brief Base class for structured block storages (the entire simulation space is partitioned into blocks [see class
*          'IBlock'], the blocks are partitioned into cells, and structured block storages manage those blocks)
*
*   The class 'StructuredBlockStorage' already provides some basic functionality that every structured block storage
*   will possess - regardless of the actual implementation of any derived class. Every structured block storage ...
*
*      - ... encapsulates an object of type BlockStorage and thereby provides the same functionality than the class
*            BlockStorage. For documentation of this functionality see class BlockStorage.
*      - ... may have multiple different grid refinement levels.
*      - ... possesses an axis-aligned cell bounding box for every grid refinement level. Each of these cell bounding
*            boxes covers the entire simulation space/domain and starts at (0,0,0). A cell bounding box is always
*            represented as a cell interval (see class 'CellInterval').
*      - ... knows the cell size dx/dy/dz (these three values may be identical) of each grid refinement level.
*      - ... provides a mechanism for adding a cell bounding box to every block as a block data "item" (see member
*            function "createCellBoundingBoxes()"). This cell bounding box represents the region within the global
*            cell coordinate space that is covered by this block. The cell bounding box is always given in relation to
*            the grid level the corresponding block is located on (-> every block is required to be assigned to one of
*            the available grid refinement levels, see member function "getLevel( const IBlock& block )").
*      - ... provides mapping mechanisms between the cell coordinate space and the underlying 3D simulation space
*      - ... provides mapping mechanisms between block local and simulation global cell coordinates (every block has its
*            own local cell coordinate space starting at (0,0,0)).
*
*   Structured block storages build upon block storages (see class 'BlockStorage'), every structured block storage
*   encapsulates an object of type BlockStorage (see member function getBlockStorage()).
*   'StructuredBlockStorage' also acts as an interface class which requires every derived class to implement a set of
*   functions that, for example, can be used to retrieve locally allocated blocks that are located at a certain position
*   within the global cell coordinate space. For more information on which functions must be implemented by a derived
*   class refer to the documentation of the public interface of 'StructuredBlockStorage'.
*
*   Attention: Immediately after constructing an instance of a class derived from StructuredBlockStorage the member
*              method "createCellBoundingBoxes()" should be called. Otherwise no block is assigned a cell bounding box
*              and all the transformation methods (-> transforming block local cells to a global cells and vice versa)
*              are not available.
*/
//**********************************************************************************************************************

class StructuredBlockStorage : private NonCopyable {

public:

   /// helper class for adding multiple block data initialization functions
   class StructuredBlockDataAdder {
   public:
      StructuredBlockDataAdder( StructuredBlockStorage & storage, const std::string & identifier = std::string() )
         : storage_( storage ), identifier_( identifier ) {}
      template< typename T >
      StructuredBlockDataAdder & operator<<( const BlockDataCreator<T> & bdc ) {
         dataHandling_.add( walberla::make_shared< internal::BlockDataHandlingHelper<T> >( bdc.dataHandling_ ), bdc.requiredSelectors_,
                            bdc.incompatibleSelectors_, bdc.identifier_ ); return *this; }
      template< typename T >
      StructuredBlockDataAdder & operator<<( const StructuredBlockDataCreator<T> & sbdc ) {
         dataHandling_.add( walberla::make_shared< internal::BlockDataHandlingHelper<T> >( walberla::make_shared< internal::BlockDataHandlingFunctionAdaptor<T> >( std::bind( sbdc.function_,  std::placeholders::_1, &storage_ ) ) ),
                            sbdc.requiredSelectors_, sbdc.incompatibleSelectors_, sbdc.identifier_ ); return *this; }
      operator BlockDataID() { return storage_.getBlockStorage().addBlockData( dataHandling_, identifier_ ); }
   private:
      StructuredBlockStorage & storage_;
      std::string identifier_;
      internal::SelectableBlockDataHandlingWrapper dataHandling_;
   };

   typedef BlockStorage::const_iterator const_iterator;
   typedef BlockStorage::iterator             iterator;

   const BlockStorage& getBlockStorage() const { return *blockStorage_; }
         BlockStorage& getBlockStorage()       { return *blockStorage_; }

    const shared_ptr<       BlockStorage > & getBlockStoragePointer()       { return blockStorage_; }
          shared_ptr< const BlockStorage >   getBlockStoragePointer() const { return blockStorage_; }

    uint_t getNumberOfLevels() const { return levels_; }

   inline bool operator==( const StructuredBlockStorage& rhs ) const;
   inline bool operator!=( const StructuredBlockStorage& rhs ) const { return !operator==( rhs ); }



   // direct access to all member functions of encapsulated class BlockStorage (-> for documentation of these functions see class BlockStorage) >>>>>>

   const AABB& getDomain() const { return blockStorage_->getDomain(); }

   bool isXPeriodic() const { return blockStorage_->isXPeriodic(); }
   bool isYPeriodic() const { return blockStorage_->isYPeriodic(); }
   bool isZPeriodic() const { return blockStorage_->isZPeriodic(); }
   bool isPeriodic( const uint_t index ) const { return blockStorage_->isPeriodic( index ); }

   void mapToPeriodicDomain( real_t & x, real_t & y, real_t & z ) const { blockStorage_->mapToPeriodicDomain(x,y,z); }
   void mapToPeriodicDomain( Vector3< real_t > & p ) const { blockStorage_->mapToPeriodicDomain(p); }

   bool periodicIntersect( const math::AABB & box1, const math::AABB & box2 ) const { return blockStorage_->periodicIntersect(box1, box2); }
   bool periodicIntersect( const math::AABB & box1, const math::AABB & box2, const real_t _dx ) const { return blockStorage_->periodicIntersect(box1, box2, _dx); }

   iterator begin( const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                   const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
      { return blockStorage_->begin( requiredSelectors, incompatibleSelectors ); }
   iterator   end() { return blockStorage_->end();   }

   const_iterator begin( const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                         const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) const
      { return blockStorage_->begin( requiredSelectors, incompatibleSelectors); }
   const_iterator   end() const { return blockStorage_->end();   }

   uint_t getNumberOfBlocks() const { return blockStorage_->getNumberOfBlocks(); }
   uint_t size()              const { return blockStorage_->size(); }
   bool   empty()             const { return blockStorage_->empty(); }

   void getBlocks( std::vector< const IBlock* >& blocks ) const { blockStorage_->getBlocks( blocks ); }
   void getBlocks( std::vector<       IBlock* >& blocks )       { blockStorage_->getBlocks( blocks ); }

   void getBlocksContainedWithinAABB( std::vector< const IBlock* >& blocks, const AABB& aabb ) const
      { blockStorage_->getBlocksContainedWithinAABB( blocks, aabb ); }
   void getBlocksContainedWithinAABB( std::vector<       IBlock* >& blocks, const AABB& aabb )
      { blockStorage_->getBlocksContainedWithinAABB( blocks, aabb ); }

   void getBlocksOverlappedByAABB( std::vector< const IBlock* >& blocks, const AABB& aabb ) const
      { blockStorage_->getBlocksOverlappedByAABB( blocks, aabb ); }
   void getBlocksOverlappedByAABB( std::vector<       IBlock* >& blocks, const AABB& aabb )
      { blockStorage_->getBlocksOverlappedByAABB( blocks, aabb ); }

   const IBlock* getBlock( const IBlockID& id ) const { return blockStorage_->getBlock( id ); }
         IBlock* getBlock( const IBlockID& id )       { return blockStorage_->getBlock( id ); }

   const IBlock* getBlock( const IBlockID::IDType& id ) const { return blockStorage_->getBlock( id ); }
         IBlock* getBlock( const IBlockID::IDType& id )       { return blockStorage_->getBlock( id ); }

   const IBlock* getBlock( const real_t x, const real_t y, const real_t z ) const { return blockStorage_->getBlock(x,y,z); }
         IBlock* getBlock( const real_t x, const real_t y, const real_t z )       { return blockStorage_->getBlock(x,y,z); }

   const IBlock* getBlock( const Vector3< real_t > & p ) const { return blockStorage_->getBlock(p); }
         IBlock* getBlock( const Vector3< real_t > & p )       { return blockStorage_->getBlock(p); }

   bool containsGlobalBlockInformation() const { return blockStorage_->containsGlobalBlockInformation(); }

   void getAllBlocks( std::vector< shared_ptr< IBlockID > >& blocks ) const { blockStorage_->getAllBlocks( blocks ); }

   bool blockExists        ( const real_t x, const real_t y, const real_t z ) const { return blockStorage_->blockExists        (x,y,z); }
   bool blockExistsLocally ( const real_t x, const real_t y, const real_t z ) const { return blockStorage_->blockExistsLocally (x,y,z); }
   bool blockExistsRemotely( const real_t x, const real_t y, const real_t z ) const { return blockStorage_->blockExistsRemotely(x,y,z); }

   bool blockExists        ( const Vector3< real_t > & p ) const { return blockStorage_->blockExists        (p); }
   bool blockExistsLocally ( const Vector3< real_t > & p ) const { return blockStorage_->blockExistsLocally (p); }
   bool blockExistsRemotely( const Vector3< real_t > & p ) const { return blockStorage_->blockExistsRemotely(p); }

   bool blockExists        ( const IBlockID & id ) const { return blockStorage_->blockExists        ( id ); }
   bool blockExistsLocally ( const IBlockID & id ) const { return blockStorage_->blockExistsLocally ( id ); }
   bool blockExistsRemotely( const IBlockID & id ) const { return blockStorage_->blockExistsRemotely( id ); }

   void getBlockID( IBlockID & id, const real_t x, const real_t y, const real_t z ) const { blockStorage_->getBlockID( id, x, y, z ); }
   void getBlockID( IBlockID & id, const Vector3< real_t > & p ) const { blockStorage_->getBlockID( id, p ); }

   AABB      getAABB       ( const IBlockID & id ) const { return blockStorage_->getAABB       ( id ); }
   Set<SUID> getState      ( const IBlockID & id ) const { return blockStorage_->getState      ( id ); }
   uint_t    getProcessRank( const IBlockID & id ) const { return blockStorage_->getProcessRank( id ); }

   void getAABB       ( AABB &      aabb,  const IBlockID & id ) const { blockStorage_->getAABB       ( aabb,  id ); }
   void getState      ( Set<SUID> & state, const IBlockID & id ) const { blockStorage_->getState      ( state, id ); }
   void getProcessRank( uint_t &    rank,  const IBlockID & id ) const { blockStorage_->getProcessRank( rank,  id ); }

   bool atDomainXMinBorder( const IBlock & block ) const { return blockStorage_->atDomainXMinBorder( block ); }
   bool atDomainXMaxBorder( const IBlock & block ) const { return blockStorage_->atDomainXMaxBorder( block ); }
   bool atDomainYMinBorder( const IBlock & block ) const { return blockStorage_->atDomainYMinBorder( block ); }
   bool atDomainYMaxBorder( const IBlock & block ) const { return blockStorage_->atDomainYMaxBorder( block ); }
   bool atDomainZMinBorder( const IBlock & block ) const { return blockStorage_->atDomainZMinBorder( block ); }
   bool atDomainZMaxBorder( const IBlock & block ) const { return blockStorage_->atDomainZMaxBorder( block ); }

   bool atDomainMinBorder( const uint_t index, const IBlock & block ) const { return blockStorage_->atDomainMinBorder( index, block ); }
   bool atDomainMaxBorder( const uint_t index, const IBlock & block ) const { return blockStorage_->atDomainMaxBorder( index, block ); }

   const std::vector< uint_t > & getNeighboringProcesses() const { return blockStorage_->getNeighboringProcesses(); }

   std::map< uint_t, std::vector< Vector3<real_t> > > getNeighboringProcessOffsets() const
      { return blockStorage_->getNeighboringProcessOffsets(); }

   internal::BlockDataHandlingAdder addBlockData( const std::string& identifier = std::string() ) { return blockStorage_->addBlockData( identifier ); }

   template< typename T >
   inline BlockDataID addBlockData( const shared_ptr< T > & dataHandling,
                                    const std::string & identifier          = std::string(),
                                    const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
   { return blockStorage_->addBlockData( dataHandling, identifier, requiredSelectors, incompatibleSelectors ); }

   template< typename T >
   inline BlockDataID addBlockData( std::function< T* ( IBlock* const block ) > function,
                                    const std::string& identifier = std::string(),
                                    const Set<SUID>& requiredSelectors     = Set<SUID>::emptySet(),
                                    const Set<SUID>& incompatibleSelectors = Set<SUID>::emptySet() )
   { return blockStorage_->addBlockData( function, identifier, requiredSelectors, incompatibleSelectors ); }

   template< typename T >
   inline BlockDataID loadBlockData( const std::string & file, const shared_ptr< T > & dataHandling,
                                     const std::string & identifier          = std::string(),
                                     const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                     const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
   { return blockStorage_->loadBlockData( file, dataHandling, identifier, requiredSelectors, incompatibleSelectors ); }
   
   void saveBlockData( const std::string & file, const BlockDataID & id ) { blockStorage_->saveBlockData( file, id ); }
   void serializeBlockData( const BlockDataID & id, mpi::SendBuffer & buffer ) { blockStorage_->serializeBlockData(id, buffer); }
   void deserializeBlockData( const BlockDataID & id, mpi::RecvBuffer & buffer )  { blockStorage_->deserializeBlockData(id, buffer); }

   void clearBlockData( const BlockDataID & id ) { blockStorage_->clearBlockData(id); }

   uint_t numberOfBlockDataItems() const { return blockStorage_->numberOfBlockDataItems(); }

   std::vector< std::string > getBlockDataIdentifiers() const { return blockStorage_->getBlockDataIdentifiers(); }
   const std::string&         getBlockDataIdentifier( const ConstBlockDataID & id ) const { return blockStorage_->getBlockDataIdentifier( id ); }

   // <<<<<< direct access to all member functions of encapsulated class BlockStorage (-> for documentation of these functions see class BlockStorage)



   // additional member functions provided by StructuredBlockStorage (some are pure virtual and must be implemented by any derived class)

   inline const CellInterval& getDomainCellBB( const uint_t level = 0 ) const; ///< the global cell coordinate space of the simulation on level "level"

   inline uint_t getNumberOfXCells( const uint_t level = 0 ) const; ///< number of global cells in x direction (within the simulation domain on level "level")
   inline uint_t getNumberOfYCells( const uint_t level = 0 ) const; ///< number of global cells in y direction (within the simulation domain on level "level")
   inline uint_t getNumberOfZCells( const uint_t level = 0 ) const; ///< number of global cells in z direction (within the simulation domain on level "level")

   /// number of global cells in x/y/z direction (within the simulation domain on level "level")
   inline uint_t getNumberOfCells( const uint_t index, const uint_t level = 0 ) const;



   real_t dx( const uint_t level = 0 ) const { WALBERLA_ASSERT_LESS( level, dx_.size() ); return dx_[ level ]; } ///< cell size on level "level" in x direction
   real_t dy( const uint_t level = 0 ) const { WALBERLA_ASSERT_LESS( level, dy_.size() ); return dy_[ level ]; } ///< cell size on level "level" in y direction
   real_t dz( const uint_t level = 0 ) const { WALBERLA_ASSERT_LESS( level, dz_.size() ); return dz_[ level ]; } ///< cell size on level "level" in z direction

   void mapToPeriodicDomain( Cell& cell, const uint_t level = 0 ) const; // -> for documentation of this function see StructuredBlockStorage.cpp



   inline Cell getCell( const real_t x, const real_t y, const real_t z, const uint_t level = 0 ) const;
   inline Cell getCell( const Vector3< real_t > & p,  const uint_t level = 0 ) const { return getCell( p[0], p[1], p[2], level ); }

   inline void getCell( Cell & cell, const real_t x, const real_t y, const real_t z, const uint_t level = 0 ) const;
   inline void getCell( Cell & cell, const Vector3< real_t > & p,  const uint_t level = 0 ) const { getCell( cell, p[0], p[1], p[2], level ); }

   inline Vector3<real_t> getCellCenter( const Cell & cell, const uint_t level = 0 ) const;
   inline void getCellCenter( real_t & x, real_t & y, real_t & z, const Cell & cell, const uint_t level = 0 ) const;
   inline void getCellCenter( Vector3< real_t > & p, const Cell & cell, const uint_t level = 0 ) const { getCellCenter( p[0], p[1], p[2], cell, level ); }

   inline AABB getCellAABB( const Cell & cell, const uint_t level = 0 ) const;
   inline void getCellAABB( AABB & aabb, const Cell & cell, const uint_t level = 0 ) const;

   inline CellInterval getCellBBFromAABB           ( const AABB& aabb, const uint_t level = 0 ) const;
   inline CellInterval getCellBBFromCellAlignedAABB( const AABB& aabb, const uint_t level = 0 ) const;

   void getCellBBFromAABB           ( CellInterval& cellBB, const AABB& aabb, const uint_t level = 0 ) const;
   void getCellBBFromCellAlignedAABB( CellInterval& cellBB, const AABB& aabb, const uint_t level = 0 ) const;

   bool isCellAlignedAABB( const AABB& aabb, const uint_t level = 0 ) const;

   inline AABB getAABBFromCellBB( const CellInterval& cellBB, const uint_t level = 0 ) const;
          void getAABBFromCellBB( AABB& aabb, const CellInterval& cellBB, const uint_t level = 0 ) const;



   inline const IBlock* getBlock( const Cell& cell, const uint_t level = 0 ) const;
   inline       IBlock* getBlock( const Cell& cell, const uint_t level = 0 );

   /// Returns true if there exists a block on level "level" at cell coordinate "cell" [Periodicity is not
   /// considered! For mapping cell coordinates to the periodic simulation space see 'mapToPeriodicDomain'.].
   /// This member function is guaranteed to work properly only if 'containsGlobalBlockInformation() == true'.
   virtual bool blockExists        ( const Cell& cell, const uint_t level = 0 ) const = 0;

   /// Returns true if locally there exists a block on level "level" at cell coordinate "cell" [Periodicity is not
   /// considered! For mapping cell coordinates to the periodic simulation space see 'mapToPeriodicDomain'.].
   /// This member function is always guaranteed to work properly, even if 'containsGlobalBlockInformation() == false'.
   virtual bool blockExistsLocally ( const Cell& cell, const uint_t level = 0 ) const = 0;

   /// Returns true if remotely there exists a block on level "level" at cell coordinate "cell" [Periodicity is not
   /// considered! For mapping cell coordinates to the periodic simulation space see 'mapToPeriodicDomain'.].
   /// This member function is guaranteed to work properly only if 'containsGlobalBlockInformation() == true'.
   virtual bool blockExistsRemotely( const Cell& cell, const uint_t level = 0 ) const = 0;

   /// Returns the block ID that corresponds to the block located on level "level" at cell coordinate "cell" [Periodicity
   /// is not considered! For mapping cell coordinates to the periodic simulation space see 'mapToPeriodicDomain'.].
   /// For local blocks, this function is always guaranteed to work. For remote blocks, this function is guaranteed to work
   /// properly only if 'containsGlobalBlockInformation() == true'.
   /// If the request cannot be satisfied (for example if no block exists on level "level" at cell coordinate "cell"), the
   /// simulation must be aborted and the call to this function must not return!
   virtual void getBlockID( IBlockID& id, const Cell& cell, const uint_t level = 0 ) const = 0;



   /// must return the level the block "block" is assigned to (must be an unsigned integer in the range [0,'number-of-levels') )
   virtual uint_t getLevel( const IBlock& block ) const = 0;

   void createCellBoundingBoxes();

   /// Returns the block data ID required for accessing the cell bounding box of blocks - fails in debug mode if no block cell bounding boxes
   /// have been created via "createCellBoundingBoxes()". (remember: every block resides on exactly one grid level, and all blocks managed by a
   //  structured block storage are assigned a corresponding cell bounding box as block data once "createCellBoundingBoxes()" is called.)
   inline ConstBlockDataID getBlockCellBBId() const { WALBERLA_ASSERT( blockCellBBCreated_ ); return blockCellBBId_; }

   inline const CellInterval& getBlockCellBB( const IBlock& block ) const;

   virtual uint_t getNumberOfXCells( const IBlock& block ) const = 0; ///< number of local cells of block "block" in x direction
   virtual uint_t getNumberOfYCells( const IBlock& block ) const = 0; ///< number of local cells of block "block" in y direction
   virtual uint_t getNumberOfZCells( const IBlock& block ) const = 0; ///< number of local cells of block "block" in z direction

   virtual uint_t getNumberOfCells( const IBlock& block, const uint_t index ) const = 0; ///< number of local cells of block "block" in x/y/z direction



   inline Cell getBlockLocalCell( const IBlock& block, const real_t x, const real_t y, const real_t z ) const;
   inline Cell getBlockLocalCell( const IBlock& block, const Vector3< real_t > & p  ) const { return getBlockLocalCell( block, p[0], p[1], p[2] ); }

   inline void getBlockLocalCell( Cell& localCell, const IBlock& block, const real_t x, const real_t y, const real_t z ) const;
   inline void getBlockLocalCell( Cell& localCell, const IBlock& block, const Vector3< real_t > & p  ) const { getBlockLocalCell( localCell, block, p[0], p[1], p[2] ); }

   inline Vector3< real_t > getBlockLocalCellCenter( const IBlock & block, const Cell & localCell ) const;
   inline void getBlockLocalCellCenter( const IBlock & block, const Cell & localCell, real_t & x, real_t & y, real_t & z ) const;
   inline void getBlockLocalCellCenter( const IBlock & block, const Cell & localCell, Vector3< real_t > & p ) const { getBlockLocalCellCenter( block, localCell, p[0], p[1], p[2] ); }
   
   inline AABB getBlockLocalCellAABB( const IBlock & block, const Cell & localCell ) const;
   inline void getBlockLocalCellAABB( const IBlock & block, const Cell & localCell, AABB & aabb ) const;

   // global <-> local

   inline void transformGlobalToBlockLocal( Vector3<real_t> & local, const IBlock& block, const Vector3<real_t> & global ) const;
   inline void transformGlobalToBlockLocal( Vector3<real_t> & point, const IBlock& block ) const;

   inline void transformBlockLocalToGlobal( Vector3<real_t> & global, const IBlock& block, const Vector3<real_t> & local ) const;
   inline void transformBlockLocalToGlobal( Vector3<real_t> &  point, const IBlock& block ) const;
   
   // global <-> local (cell)

   inline void transformGlobalToBlockLocalCell( Cell& local, const IBlock& block, const Cell& global ) const;
   inline void transformGlobalToBlockLocalCell( Cell& cell,  const IBlock& block ) const;

   inline void transformBlockLocalToGlobalCell( Cell& global, const IBlock& block, const Cell& local ) const;
   inline void transformBlockLocalToGlobalCell( Cell& cell,   const IBlock& block ) const;

   // global <-> local (cell interval/bounding box)

   inline void transformGlobalToBlockLocalCellInterval( CellInterval& local,    const IBlock& block, const CellInterval& global ) const;
   inline void transformGlobalToBlockLocalCellInterval( CellInterval& interval, const IBlock& block ) const;

   inline void transformBlockLocalToGlobalCellInterval( CellInterval& global,   const IBlock& block, const CellInterval& local ) const;
   inline void transformBlockLocalToGlobalCellInterval( CellInterval& interval, const IBlock& block ) const;



   /// Must be used if multiple initialization functions with different selection attributes are registered for initializing
   /// the same block data "item". A 'BlockDataID' that corresponds to this block data "item" is returned.
   /// "BlockDataCreator" and "StructuredBlockDataCreator" structures may be mixed.
   ///
   /// Usage: BlockDataID id = blockStorage.addBlockData( "[optional block data identifier]" ) << StructuredBlockDataCreator( ... )
   ///                                                                                         << BlockDataCreator( ... ) << ... ;
   StructuredBlockDataAdder addStructuredBlockData( const std::string& identifier = std::string() )
                                                                               { return StructuredBlockDataAdder( *this, identifier ); }

   template< typename T >
   inline BlockDataID addStructuredBlockData( std::function< T* ( IBlock* const block, StructuredBlockStorage* const storage ) > function,
                                              const std::string& identifier          = std::string(),
                                              const Set<SUID>& requiredSelectors     = Set<SUID>::emptySet(),
                                              const Set<SUID>& incompatibleSelectors = Set<SUID>::emptySet() );

protected:

   /// Every derived class must call this constructor!
   StructuredBlockStorage( const shared_ptr<BlockStorage> & blockStorage,
                           const std::vector< uint_t > & xCells, const std::vector< uint_t > & yCells, const std::vector< uint_t > & zCells );

   virtual ~StructuredBlockStorage() {} ///< Must not be made public! No one should be allowed to delete a variable of type 'StructuredBlockStorage*'

   void resetCellDecomposition( const std::vector< uint_t > & xCells, const std::vector< uint_t > & yCells, const std::vector< uint_t > & zCells );

   virtual bool equal( const StructuredBlockStorage* rhs ) const = 0;

   // helper class for 'StructuredBlockStorage::addCellBoundingBoxesAsBlockData'
   class CellBoundingBoxHandling : public AlwaysInitializeBlockDataHandling< CellInterval >
   {
   public:
      CellBoundingBoxHandling( const StructuredBlockStorage & storage ) : storage_( storage ) {}
      CellInterval * initialize( IBlock * const block ) { return storage_.initializeCellBoundingBox( block ); }
   private:
      const StructuredBlockStorage & storage_;
   };
   friend class CellBoundingBoxHandling;

   CellInterval * initializeCellBoundingBox( IBlock * const block ) const;
   virtual inline BlockDataID addCellBoundingBoxesAsBlockData( const std::string & identifier );

private:

   StructuredBlockStorage(); ///< Must not be made public or protected! Derived classes must call one of the available public/protected constructors.



   shared_ptr<BlockStorage> blockStorage_; ///< reference to an encapsulated object of type class BlockStorage (the object itself must be stored as a member in the derived class)

   uint_t levels_; ///< number of different grid levels managed by this block storage (every grid level has its own cell size dx/dy/dz)

   std::vector< CellInterval > domainCellBB_; ///< cell bounding box for the entire simulation domain for every level (always starts at (0,0,0)) [min & max included!]

   std::vector< real_t > dx_; ///< for every level: ( domain x width ) / ( number of cells in x direction )
   std::vector< real_t > dy_; ///< for every level: ( domain y width ) / ( number of cells in y direction )
   std::vector< real_t > dz_; ///< for every level: ( domain z width ) / ( number of cells in z direction )

   bool        blockCellBBCreated_; ///< flag for checking whether or not a cell bounding box has already been added as a block data "item" to all blocks
   BlockDataID blockCellBBId_;      ///< block data ID for accessing a block's cell bounding box

}; // class StructuredBlockStorage



/// The following members are not used for checking if two StructuredBlockStorage objects are equal: blockCellBBId_
inline bool StructuredBlockStorage::operator==( const StructuredBlockStorage& rhs ) const {

   if( blockStorage_ != rhs.blockStorage_ || levels_ != rhs.levels_ )
      return false;

   for( uint_t i = 0; i != levels_; ++i ) {
      if( domainCellBB_[i] != rhs.domainCellBB_[i] ||
          !realIsEqual( dx_[i], rhs.dx_[i] ) || !realIsEqual( dy_[i], rhs.dy_[i] ) || !realIsEqual( dz_[i], rhs.dz_[i] ) )
         return false;
   }

   return blockCellBBCreated_ == rhs.blockCellBBCreated_ && equal( &rhs );
}



inline const CellInterval& StructuredBlockStorage::getDomainCellBB( const uint_t level ) const {

   WALBERLA_ASSERT_LESS( level, domainCellBB_.size() );

   return domainCellBB_[ level ];
}



inline uint_t StructuredBlockStorage::getNumberOfXCells( const uint_t level ) const {

   WALBERLA_ASSERT_LESS( level, domainCellBB_.size() );

   return uint_c( domainCellBB_[ level ].xMax() + 1 );
}



inline uint_t StructuredBlockStorage::getNumberOfYCells( const uint_t level ) const {

   WALBERLA_ASSERT_LESS( level, domainCellBB_.size() );

   return uint_c( domainCellBB_[ level ].yMax() + 1 );
}



inline uint_t StructuredBlockStorage::getNumberOfZCells( const uint_t level ) const {

   WALBERLA_ASSERT_LESS( level, domainCellBB_.size() );

   return uint_c( domainCellBB_[ level ].zMax() + 1 );
}



inline uint_t StructuredBlockStorage::getNumberOfCells( const uint_t index, const uint_t level ) const {

   WALBERLA_ASSERT_LESS( index, uint_t(3) );
   WALBERLA_ASSERT_LESS( level, domainCellBB_.size() );

   return uint_c( domainCellBB_[ level ].max()[ index ] + 1 );
}



//**********************************************************************************************************************
/*!
*   For documentation see member function
*   "void getCell( Cell & cell, const real_t x, const real_t y, const real_t z, const uint_t level ) const"
*/
//**********************************************************************************************************************
inline Cell StructuredBlockStorage::getCell( const real_t x, const real_t y, const real_t z, const uint_t level ) const
{
   Cell cell;
   getCell( cell, x, y, z, level );
   return cell;
}



//**********************************************************************************************************************
/*!
*   \brief Maps any point in 3D space to a cell coordinate within the global cell coordinate space on level "level".
*
*   Points located on any of the three "min" faces of a cell (= the three faces that include the lower left/min corner
*   of the cell) are considered to be part of the cell - points located on any of the three "max" faces are NOT!
*
*   [Attention: Periodicity is not considered! For mapping cell coordinates to the periodic simulation space see
*    member function 'mapToPeriodicDomain'.]
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::getCell( Cell & cell, const real_t x, const real_t y, const real_t z, const uint_t level ) const
{
   const AABB& domain = getDomain();

   cell.x() = cell_idx_c( std::floor( ( x - domain.xMin() ) / dx( level ) ) );
   cell.y() = cell_idx_c( std::floor( ( y - domain.yMin() ) / dy( level ) ) );
   cell.z() = cell_idx_c( std::floor( ( z - domain.zMin() ) / dz( level ) ) );
}



//**********************************************************************************************************************
/*!
*   For documentation see member function
*   "void getCellCenter( real_t & x, real_t & y, real_t & z, const Cell & cell, const uint_t level ) const"
*/
//**********************************************************************************************************************
inline Vector3<real_t> StructuredBlockStorage::getCellCenter( const Cell & cell, const uint_t level ) const
{
   Vector3<real_t> center;
   getCellCenter( center[0], center[1], center[2], cell, level );
   return center;
}



//**********************************************************************************************************************
/*!
*   Returns the location (within the global 3D simulation space) of the center of cell "cell" on level "level".
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::getCellCenter( real_t & x, real_t & y, real_t & z, const Cell & cell, const uint_t level ) const
{
   const AABB& domain = getDomain();

   x = domain.xMin() + ( real_c( cell.x() ) + real_c(0.5) ) * dx( level );
   y = domain.yMin() + ( real_c( cell.y() ) + real_c(0.5) ) * dy( level );
   z = domain.zMin() + ( real_c( cell.z() ) + real_c(0.5) ) * dz( level );
}



//**********************************************************************************************************************
/*!
*   For documentation see member function "void getCellAABB( AABB& aabb, const Cell& cell, const uint_t level ) const"
*/
//**********************************************************************************************************************
inline AABB StructuredBlockStorage::getCellAABB( const Cell & cell, const uint_t level ) const
{
   AABB aabb;
   getCellAABB( aabb, cell, level );
   return aabb;
}



//**********************************************************************************************************************
/*!
*   Returns the axis-aligned bounding box (with respect to the global 3D simulation space) that covers the cell "cell"
*   on level "level".
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::getCellAABB( AABB& aabb, const Cell& cell, const uint_t level ) const {

   const AABB& domain = getDomain();

   const real_t x = domain.xMin() + real_c( cell.x() ) * dx( level );
   const real_t y = domain.yMin() + real_c( cell.y() ) * dy( level );
   const real_t z = domain.zMin() + real_c( cell.z() ) * dz( level );

   aabb.init( x, y, z, x + dx( level ), y + dy( level ), z + dz( level ) );
}



//**********************************************************************************************************************
/*!
*   For documentation see member function
*   "void getCellBBFromAABB( CellInterval& cellBB, const AABB& aabb, const uint_t level ) const"
*/
//**********************************************************************************************************************
inline CellInterval StructuredBlockStorage::getCellBBFromAABB( const AABB& aabb, const uint_t level ) const
{
   CellInterval interval;
   getCellBBFromAABB( interval, aabb, level );
   return interval;
}



//**********************************************************************************************************************
/*!
*   For documentation see member function
*   "void getCellBBFromCellAlignedAABB( CellInterval& cellBB, const AABB& aabb, const uint_t level ) const"
*/
//**********************************************************************************************************************
inline CellInterval StructuredBlockStorage::getCellBBFromCellAlignedAABB( const AABB& aabb, const uint_t level ) const
{
   CellInterval interval;
   getCellBBFromCellAlignedAABB( interval, aabb, level );
   return interval;
}



//**********************************************************************************************************************
/*!
*   For documentation see member function
*   "void getAABBFromCellBB( AABB& aabb, const CellInterval& cellBB, const uint_t level ) const"
*/
//**********************************************************************************************************************
inline AABB StructuredBlockStorage::getAABBFromCellBB( const CellInterval& cellBB, const uint_t level ) const
{
   AABB aabb;
   getAABBFromCellBB( aabb, cellBB, level );
   return aabb;
}



//**********************************************************************************************************************
/*!
*   \brief Returns the block located at global cell coordinate "cell" on grid level "level" (returns 'NULL' if the
*          block doesn't exist locally!).
*
*   Attention: Periodicity is not considered! For mapping cell coordinates to the periodic simulation space see member
*              function 'mapToPeriodicDomain'.
*/
//**********************************************************************************************************************
inline const IBlock* StructuredBlockStorage::getBlock( const Cell& cell, const uint_t level ) const {

   real_t x, y, z;
   getCellCenter( x, y, z, cell, level );

   const IBlock* block = blockStorage_->getBlock(x,y,z);
   if( block == NULL )
      return NULL;

   return ( getLevel( *block ) == level ) ? block : NULL;
}



//**********************************************************************************************************************
/*!
*   \brief Returns the block located at global cell coordinate "cell" on grid level "level" (returns 'NULL' if the
*          block doesn't exist locally!).
*
*   Attention: Periodicity is not considered! For mapping cell coordinates to the periodic simulation space see member
*              function 'mapToPeriodicDomain'.
*/
//**********************************************************************************************************************
inline IBlock* StructuredBlockStorage::getBlock( const Cell& cell, const uint_t level ) {

   real_t x, y, z;
   getCellCenter( x, y, z, cell, level );

   IBlock* block = blockStorage_->getBlock(x,y,z);
   if( block == NULL )
      return NULL;

   return ( getLevel( *block ) == level ) ? block : NULL;
}



//**********************************************************************************************************************
/*!
*   \brief Returns the cell bounding box of block "block"
*
*   In debug mode: fails if no block cell bounding boxes have been created via "createCellBoundingBoxes()".
*
*   Remember: Every block resides on exactly one grid level, and all blocks managed by a structured block storage are
*             assigned a corresponding cell bounding box as block data once "createCellBoundingBoxes()" is called.
*/
//**********************************************************************************************************************
inline const CellInterval& StructuredBlockStorage::getBlockCellBB( const IBlock& block ) const
{
   WALBERLA_ASSERT_EQUAL( blockStorage_.get(), &(block.getBlockStorage()) );
   WALBERLA_ASSERT( blockCellBBCreated_ );

   return *(block.uncheckedFastGetData< CellInterval >( blockCellBBId_ ));
}



//**********************************************************************************************************************
/*!
*   For documentation see member function
*   "void getBlockLocalCell( Cell& localCell, const IBlock& block, const real_t x, const real_t y, const real_t z ) const"
*/
//**********************************************************************************************************************
inline Cell StructuredBlockStorage::getBlockLocalCell( const IBlock& block, const real_t x, const real_t y, const real_t z ) const
{
   Cell cell;
   getBlockLocalCell( cell, block, x, y, z );
   return cell;
}



//**********************************************************************************************************************
/*!
*   \brief Maps any point in 3D space to a cell coordinate within the local cell coordinate space of the block "block".
*
*   Points located on any of the three "min" faces of a cell (= the three faces that include the lower left/min corner
*   of the cell) are considered to be part of the cell - points located on any of the three "max" faces are NOT!
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::getBlockLocalCell( Cell& localCell, const IBlock& block, const real_t x, const real_t y, const real_t z ) const
{
   WALBERLA_ASSERT_EQUAL( blockStorage_.get(), &(block.getBlockStorage()) );

   const AABB & aabb  = block.getAABB();
   const uint_t level = getLevel( block );

   WALBERLA_ASSERT_LESS( level, levels_ );

   localCell.x() = cell_idx_c( std::floor( ( x - aabb.xMin() ) / dx( level ) ) );
   localCell.y() = cell_idx_c( std::floor( ( y - aabb.yMin() ) / dy( level ) ) );
   localCell.z() = cell_idx_c( std::floor( ( z - aabb.zMin() ) / dz( level ) ) );
}



//**********************************************************************************************************************
/*!
*   For documentation see member function
*   "void getBlockLocalCellCenter( const IBlock & block, const Cell & localCell, real_t & x, real_t & y, real_t & z ) const"
*/
//**********************************************************************************************************************
inline Vector3< real_t > StructuredBlockStorage::getBlockLocalCellCenter( const IBlock & block, const Cell & localCell ) const
{
   Vector3< real_t > center;
   getBlockLocalCellCenter( block, localCell, center[0], center[1], center[2] );
   return center;
}



//**********************************************************************************************************************
/*!
*   Returns the location (within the global 3D simulation space) of the center of block "block"s local cell "localCell".
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::getBlockLocalCellCenter( const IBlock & block, const Cell & localCell, real_t & x, real_t & y, real_t & z ) const
{
   WALBERLA_ASSERT_EQUAL( blockStorage_.get(), &(block.getBlockStorage()) );

   const AABB & aabb  = block.getAABB();
   const uint_t level = getLevel( block );

   WALBERLA_ASSERT_LESS( level, levels_ );

   x = aabb.xMin() + ( real_c( localCell.x() ) + real_c(0.5) ) * dx( level );
   y = aabb.yMin() + ( real_c( localCell.y() ) + real_c(0.5) ) * dy( level );
   z = aabb.zMin() + ( real_c( localCell.z() ) + real_c(0.5) ) * dz( level );
}



//**********************************************************************************************************************
/*!
*   For documentation see member function
*   "void getBlockLocalCellAABB( const IBlock & block, const Cell & localCell, AABB & aabb ) const"
*/
//**********************************************************************************************************************
inline AABB StructuredBlockStorage::getBlockLocalCellAABB( const IBlock & block, const Cell & localCell ) const
{
   AABB aabb;
   getBlockLocalCellAABB( block, localCell, aabb );
   return aabb;
}



//**********************************************************************************************************************
/*!
*   Returns the axis-aligned bounding box (with respect to the global 3D simulation space) that covers block "block"s
*   local cell "localCell"
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::getBlockLocalCellAABB( const IBlock & block, const Cell & localCell, AABB & aabb ) const
{
   WALBERLA_ASSERT_EQUAL( blockStorage_.get(), &(block.getBlockStorage()) );

   const AABB& blockAABB = block.getAABB();
   const uint_t level = getLevel( block );

   WALBERLA_ASSERT_LESS( level, levels_ );

   const real_t x = blockAABB.xMin() + real_c( localCell.x() ) * dx( level );
   const real_t y = blockAABB.yMin() + real_c( localCell.y() ) * dy( level );
   const real_t z = blockAABB.zMin() + real_c( localCell.z() ) * dz( level );

   aabb.init( x, y, z, x + dx( level ), y + dy( level ), z + dz( level ) );
}



//**********************************************************************************************************************
/*!
*   Transforms a global point "global" (assumed to be located on the same grid level than the block "block"
*   resides on) into the block local point "local".
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::transformGlobalToBlockLocal( Vector3<real_t> & local, const IBlock& block, const Vector3<real_t> & global ) const
{
   WALBERLA_ASSERT_EQUAL( blockStorage_.get(), &(block.getBlockStorage()) );

   const uint_t level = getLevel( block );

   local = global - block.getAABB().minCorner();
   local[0] /= dx( level );
   local[1] /= dy( level );
   local[2] /= dz( level );
}



//**********************************************************************************************************************
/*!
*   Transforms a global point (assumed to be located on the same grid level than the block "block" resides on)
*   into a block local point.
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::transformGlobalToBlockLocal( Vector3<real_t> & point, const IBlock& block ) const
{
   WALBERLA_ASSERT_EQUAL( blockStorage_.get(), &(block.getBlockStorage()) );

   const uint_t level = getLevel( block );

   point -= block.getAABB().minCorner();
   point[0] /= dx( level );
   point[1] /= dy( level );
   point[2] /= dz( level );
}



//**********************************************************************************************************************
/*!
*   Transforms the block local point "local" into the global point "global" (the global point is
*   given with respect to the grid level the block "block" resides on).
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::transformBlockLocalToGlobal( Vector3<real_t> & global, const IBlock& block, const Vector3<real_t> & local ) const
{
   WALBERLA_ASSERT_EQUAL( blockStorage_.get(), &(block.getBlockStorage()) );

   const uint_t level = getLevel( block );

   global[0] = local[0] * dx( level );
   global[1] = local[1] * dy( level );
   global[2] = local[2] * dz( level );

   global += block.getAABB().minCorner();
}



//**********************************************************************************************************************
/*!
*   Transforms a block local point into a global point (the global point is given with respect to the
*   grid level the block "block" resides on).
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::transformBlockLocalToGlobal( Vector3<real_t> &  point, const IBlock& block ) const
{
   WALBERLA_ASSERT_EQUAL( blockStorage_.get(), &(block.getBlockStorage()) );

   const uint_t level = getLevel( block );
      
   point[0] *= dx( level );
   point[1] *= dy( level );
   point[2] *= dz( level );

   point += block.getAABB().minCorner();
}




//**********************************************************************************************************************
/*!
*   Transforms the global cell coordinate "global" (assumed to be located on the same grid level than the block "block"
*   resides on) into the block local cell coordinate "local".
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::transformGlobalToBlockLocalCell( Cell& local, const IBlock& block, const Cell& global ) const {

   WALBERLA_ASSERT_EQUAL( blockStorage_.get(), &(block.getBlockStorage()) );

   const CellInterval& cellBB = getBlockCellBB( block );

   local.x() = global.x() - cellBB.xMin();
   local.y() = global.y() - cellBB.yMin();
   local.z() = global.z() - cellBB.zMin();
}



//**********************************************************************************************************************
/*!
*   Transforms a global cell coordinate (assumed to be located on the same grid level than the block "block" resides on)
*   into a block local cell coordinate.
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::transformGlobalToBlockLocalCell( Cell& cell, const IBlock& block ) const {

   transformGlobalToBlockLocalCell( cell, block, cell );
}



//**********************************************************************************************************************
/*!
*   Transforms the block local cell coordinate "local" into the global cell coordinate "global" (the global cell is
*   given with respect to the grid level the block "block" resides on).
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::transformBlockLocalToGlobalCell( Cell& global, const IBlock& block, const Cell& local ) const {

   WALBERLA_ASSERT_EQUAL( blockStorage_.get(), &(block.getBlockStorage()) );

   const CellInterval& cellBB = getBlockCellBB( block );

   global.x() = local.x() + cellBB.xMin();
   global.y() = local.y() + cellBB.yMin();
   global.z() = local.z() + cellBB.zMin();
}



//**********************************************************************************************************************
/*!
*   Transforms a block local cell coordinate into a global cell coordinate (the global cell is given with respect to the
*   grid level the block "block" resides on).
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::transformBlockLocalToGlobalCell( Cell& cell, const IBlock& block ) const {

   transformBlockLocalToGlobalCell( cell, block, cell );
}



//**********************************************************************************************************************
/*!
*   Transforms the global cell interval "global" (assumed to be located on the same grid level than the block "block"
*   resides on) into the block local cell interval "local".
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::transformGlobalToBlockLocalCellInterval( CellInterval& local, const IBlock& block,
                                                                             const CellInterval& global ) const {

   transformGlobalToBlockLocalCell( local.min(), block, global.min() );
   transformGlobalToBlockLocalCell( local.max(), block, global.max() );
}



//**********************************************************************************************************************
/*!
*   Transforms a global cell interval (assumed to be located on the same grid level than the block "block" resides on)
*   into a block local cell interval.
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::transformGlobalToBlockLocalCellInterval( CellInterval& interval, const IBlock& block ) const {

   transformGlobalToBlockLocalCell( interval.min(), block );
   transformGlobalToBlockLocalCell( interval.max(), block );
}



//**********************************************************************************************************************
/*!
*   Transforms the block local cell interval "local" into the global cell interval "global" (the global cell interval is
*   given with respect to the grid level the block "block" resides on).
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::transformBlockLocalToGlobalCellInterval( CellInterval& global, const IBlock& block,
                                                                             const CellInterval& local ) const {

   transformBlockLocalToGlobalCell( global.min(), block, local.min() );
   transformBlockLocalToGlobalCell( global.max(), block, local.max() );
}



//**********************************************************************************************************************
/*!
*   Transforms a block local cell interval into a global cell interval (the global cell interval is given with respect
*   to the grid level the block "block" resides on).
*/
//**********************************************************************************************************************
inline void StructuredBlockStorage::transformBlockLocalToGlobalCellInterval( CellInterval& interval, const IBlock& block ) const {

   transformBlockLocalToGlobalCell( interval.min(), block );
   transformBlockLocalToGlobalCell( interval.max(), block );
}



//**********************************************************************************************************************
/*!
*   This function can be used for initializing a new block data "item". A 'BlockDataID' that corresponds to this block
*   data "item" is returned. This block data ID must be used to access/retrieve block data (see function 'getData' of
*   class 'IBlock').
*   If multiple initialization functions with different selection attributes shall be used for initializing the same
*   block data "item" see member function "addStructuredBlockData( const std::string& identifier )"
*/
//**********************************************************************************************************************
template< typename T >
inline BlockDataID StructuredBlockStorage::addStructuredBlockData(
      std::function< T* ( IBlock* const block, StructuredBlockStorage* const storage ) > function,
      const std::string& identifier, const Set<SUID>& requiredSelectors, const Set<SUID>& incompatibleSelectors )
{
   internal::SelectableBlockDataHandlingWrapper dataHandling(
            walberla::make_shared< internal::BlockDataHandlingHelper<T> >( walberla::make_shared< internal::BlockDataHandlingFunctionAdaptor<T> >( std::bind( function,  std::placeholders::_1, this ) ) ),
            requiredSelectors, incompatibleSelectors, identifier );


   return blockStorage_->addBlockData( dataHandling, identifier );
}



inline BlockDataID StructuredBlockStorage::addCellBoundingBoxesAsBlockData( const std::string & identifier )
{
   return addBlockData( walberla::make_shared< CellBoundingBoxHandling >( *this ), identifier );
}



} // namespace domain_decomposition

using domain_decomposition::StructuredBlockDataCreator;
using domain_decomposition::StructuredBlockStorage;

} // namespace walberla
