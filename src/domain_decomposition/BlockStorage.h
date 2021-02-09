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
//! \file BlockStorage.h
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "BlockDataHandling.h"
#include "IBlock.h"

#include "core/DataTypes.h"
#include "core/NonCopyable.h"
#include "core/Set.h"
#include "core/debug/Debug.h"
#include "core/math/Vector3.h"
#include "core/mpi/MPIWrapper.h"
#include "core/selectable/IsSetSelected.h"
#include "core/uid/GlobalState.h"
#include "core/uid/SUID.h"

#include <functional>
#include <string>
#include <vector>


namespace walberla {
namespace domain_decomposition {



//**********************************************************************************************************************
/*!
*   \brief Base class for block storage data structures (the entire simulation space is partitioned into blocks [see
*          class 'IBlock'], and block storage data structures manage these blocks)
*
*   The class 'BlockStorage' already provides some basic functionality that every block storage data structure will
*   possess - regardless of the actual implementation of any derived class. Every block storage data structure ...
*
*      - ... possesses an axis-aligned bounding box that covers the entire simulation space/domain.
*      - ... holds information about whether or not the simulation is periodic in x-/y-/z-direction.
*      - ... provides iterators for traversing all locally allocated blocks (that is all blocks that are assigned to
*            this process).
*      - ... provides a mechanism for registering block data handling objects. These objects are used to add
*            any kind of data to blocks.
*
*   'BlockStorage' also acts as an interface class which requires every derived class to implement a set of functions
*   that can be used to retrieve locally allocated blocks that match certain criteria (for example all blocks that are
*   contained within a certain region of the simulation space). For more information on which functions must be
*   implemented by a derived class refer to the documentation of the public interface of 'BlockStorage'.
*/
//**********************************************************************************************************************

class BlockStorage : private NonCopyable {

   friend class IBlock;
   friend class StructuredBlockStorage;

public:

   typedef std::map< IBlockID::IDType, IBlock* > BlockContainerType ;

   class const_iterator;

   class iterator {
      friend class const_iterator;
      friend class BlockStorage;
   public:
      using iterator_category = std::forward_iterator_tag;
      using value_type = IBlock;
      using difference_type = std::ptrdiff_t;
      using pointer = IBlock*;
      using reference = IBlock&;

      iterator( const iterator & it )  = default;

      iterator & operator++()    { ++it_; checkStateAndAdapt(); return *this; }      // prefix ++X
      iterator   operator++(int) { iterator it( *this ); operator++(); return it; }; // postfix X++

      bool operator==( const iterator & rhs ) const { return it_ == rhs.it_; }
      bool operator!=( const iterator & rhs ) const { return it_ != rhs.it_; }

      IBlock * get() { return it_->second; }

      IBlock & operator*()  { return *(it_->second); }
      IBlock * operator->() { return   it_->second; }

   private:

      iterator( BlockContainerType::iterator it, BlockContainerType::iterator end,
                const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
         it_( it ), end_( end ), requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors ) { checkStateAndAdapt(); }

      void checkStateAndAdapt()
      {
         if( it_ != end_ && !( requiredSelectors_.empty() && incompatibleSelectors_.empty() ) )
         {
            while( it_ != end_ && !selectable::isSetSelected( it_->second->getState() + uid::globalState(), requiredSelectors_, incompatibleSelectors_ ) )
            {
               ++it_;
            }
         }
      }

      BlockContainerType::iterator it_;
      BlockContainerType::iterator end_;

      Set<SUID> requiredSelectors_;
      Set<SUID> incompatibleSelectors_;
   };


   class const_iterator {
      friend class BlockStorage;
   public:

      const_iterator( const iterator & it ) :
         it_( it.it_ ), end_( it.end_ ), requiredSelectors_( it.requiredSelectors_ ), incompatibleSelectors_( it.incompatibleSelectors_ ) {}
      const_iterator( const const_iterator & it )  = default;

      const_iterator & operator++()    { ++it_; checkStateAndAdapt(); return *this; }            // prefix ++X
      const_iterator   operator++(int) { const_iterator it( *this ); operator++(); return it; }; // postfix X++

      bool operator==( const const_iterator & rhs ) const { return it_ == rhs.it_; }
      bool operator!=( const const_iterator & rhs ) const { return it_ != rhs.it_; }

      const IBlock * get() const { return it_->second; }

      const IBlock & operator*()  const { return *(it_->second); }
      const IBlock * operator->() const { return  it_->second; }

   private:

      const_iterator( BlockContainerType::const_iterator it, BlockContainerType::const_iterator end,
                      const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                      const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet()) :
         it_( it ), end_( end ), requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors ) { checkStateAndAdapt(); }

      void checkStateAndAdapt()
      {
         if( it_ != end_ && !( requiredSelectors_.empty() && incompatibleSelectors_.empty() ) )
         {
            while( it_ != end_ && !selectable::isSetSelected( it_->second->getState() + uid::globalState(), requiredSelectors_, incompatibleSelectors_ ) )
            {
               ++it_;
            }
         }
      }

      BlockContainerType::const_iterator it_;
      BlockContainerType::const_iterator end_;

      Set<SUID> requiredSelectors_;
      Set<SUID> incompatibleSelectors_;
   };



   const AABB& getDomain() const { return domain_; } ///< returns an axis-aligned bounding box that covers the entire simulation space/domain

   bool isXPeriodic() const { return periodic_[0]; }
   bool isYPeriodic() const { return periodic_[1]; }
   bool isZPeriodic() const { return periodic_[2]; }
   bool isPeriodic( const uint_t index ) const { WALBERLA_ASSERT_LESS( index, uint_t(3) ); return periodic_[index]; }

   void mapToPeriodicDomain( real_t & x, real_t & y, real_t & z ) const; // -> for documentation of this function see BlockStorage.cpp
   void mapToPeriodicDomain( Vector3< real_t > & p ) const { mapToPeriodicDomain( p[0], p[1], p[2] ); }

   bool periodicIntersect( const math::AABB & box1, const math::AABB & box2 ) const; // -> for documentation of this function see BlockStorage.cpp
   bool periodicIntersect( const math::AABB & box1, const math::AABB & box2, const real_t dx ) const; // -> for documentation of this function see BlockStorage.cpp


   inline bool operator==( const BlockStorage& rhs ) const;
   inline bool operator!=( const BlockStorage& rhs ) const { return !operator==( rhs ); }

   /// iterator for traversing all locally allocated blocks
   iterator begin( const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                   const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
      { return iterator( iBlocks_.begin(), iBlocks_.end(), requiredSelectors, incompatibleSelectors ); }
   iterator   end() { return iterator( iBlocks_.end(), iBlocks_.end() ); } ///< iterator for traversing all locally allocated blocks

   /// iterator for traversing all locally allocated blocks
   const_iterator begin( const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                         const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) const
      { return const_iterator( iBlocks_.begin(), iBlocks_.end(), requiredSelectors, incompatibleSelectors ); }
   const_iterator   end() const { return const_iterator( iBlocks_.end(), iBlocks_.end() ); } ///< iterator for traversing all locally allocated blocks

   uint_t getNumberOfBlocks() const { return iBlocks_.size(); } ///< number of locally allocated blocks
   uint_t size()              const { return iBlocks_.size(); } ///< number of locally allocated blocks
   bool   empty()             const { return iBlocks_.empty(); }

   /// inserts all locally allocated blocks into vector 'blocks'
   void getBlocks( std::vector< const IBlock* >& blocks ) const
   { for (auto blkIt = iBlocks_.begin(); blkIt != iBlocks_.end(); ++blkIt ) { blocks.push_back(blkIt->second); } }

   /// inserts all locally allocated blocks into vector 'blocks'
   void getBlocks( std::vector<       IBlock* >& blocks )
   { for (auto blkIt = iBlocks_.begin(); blkIt != iBlocks_.end(); ++blkIt ) { blocks.push_back(blkIt->second); } }

   //*******************************************************************************************************************
   /*!
   *   Calling this function causes all locally allocated blocks that are completely contained within the axis-aligned
   *   bounding box 'aabb' to be inserted at the end of vector 'blocks'. The behavior of this function is identical to:
   *
   *   \code
   *     for( const_iterator block = begin(); block != end(); ++block )
   *        if( aabb.contains( block->getAABB() ) ) blocks.push_back( *block );
   *   \endcode
   */
   //*******************************************************************************************************************
   virtual void getBlocksContainedWithinAABB( std::vector< const IBlock* >& blocks, const AABB& aabb ) const = 0;

   //*******************************************************************************************************************
   /*!
   *   Calling this function causes all locally allocated blocks that are completely contained within the axis-aligned
   *   bounding box 'aabb' to be inserted at the end of vector 'blocks'. The behavior of this function is identical to:
   *
   *   \code
   *     for( iterator block = begin(); block != end(); ++block )
   *        if( aabb.contains( block->getAABB() ) ) blocks.push_back( *block );
   *   \endcode
   */
   //*******************************************************************************************************************
   virtual void getBlocksContainedWithinAABB( std::vector<       IBlock* >& blocks, const AABB& aabb )       = 0;

   //*******************************************************************************************************************
   /*!
   *   Calling this function causes all locally allocated blocks that intersect with the axis-aligned bounding box
   *   'aabb' to be inserted at the end of vector 'blocks'. The behavior of this function is identical to:
   *
   *   \code
   *     for( const_iterator block = begin(); block != end(); ++block )
   *        if( aabb.intersects( block->getAABB() ) ) blocks.push_back( *block );
   *   \endcode
   */
   //*******************************************************************************************************************
   virtual void getBlocksOverlappedByAABB( std::vector< const IBlock* >& blocks, const AABB& aabb ) const = 0;

   //*******************************************************************************************************************
   /*!
   *   Calling this function causes all locally allocated blocks that intersect with the axis-aligned bounding box
   *   'aabb' to be inserted at the end of vector 'blocks'. The behavior of this function is identical to:
   *
   *   \code
   *     for( iterator block = begin(); block != end(); ++block )
   *        if( aabb.intersects( block->getAABB() ) ) blocks.push_back( *block );
   *   \endcode
   */
   //*******************************************************************************************************************
   virtual void getBlocksOverlappedByAABB( std::vector<       IBlock* >& blocks, const AABB& aabb )       = 0;

   //*******************************************************************************************************************
   /*!
   *   \brief Returns the block associated with block ID 'id' (returns 'NULL' if the block doesn't exist locally!). The
   *   behavior of this function is identical to:
   *
   *   \code
   *     for( const_iterator block = begin(); block != end(); ++block )
   *        if( block->getId() == id ) return *block;
   *     return NULL;
   *   \endcode
   */
   //*******************************************************************************************************************
   virtual const IBlock* getBlock( const IBlockID& id ) const = 0;
   inline  const IBlock* getBlock( const IBlockID::IDType& id) const { 
     auto it = iBlocks_.find(id);
    
     if( it != iBlocks_.end()) return it->second;

     return nullptr;
   }

   //*******************************************************************************************************************
   /*!
   *   \brief Returns the block associated with block ID 'id' (returns 'NULL' if the block doesn't exist locally!). The
   *   behavior of this function is identical to:
   *
   *   \code
   *     for( iterator block = begin(); block != end(); ++block )
   *        if( block->getId() == id ) return *block;
   *     return NULL;
   *   \endcode
   */
   //*******************************************************************************************************************
   virtual       IBlock* getBlock( const IBlockID& id )       = 0;
   inline        IBlock* getBlock( const IBlockID::IDType& id) { 
     auto it = iBlocks_.find(id);
    
     if( it != iBlocks_.end()) return it->second;

     return nullptr;
   }

   //*******************************************************************************************************************
   /*!
   *   \brief Returns the block located at position (x,y,z) (returns 'NULL' if the block doesn't exist locally!). The
   *   behavior of this function is identical to:
   *
   *   \code
   *     for( const_iterator block = begin(); block != end(); ++block )
   *        if( block->getAABB().contains(x,y,z) ) return *block;
   *     return NULL;
   *   \endcode
   *
   *   Periodicity is not considered! For mapping points to the periodic simulation space see 'mapToPeriodicDomain'.
   */
   //*******************************************************************************************************************
   virtual const IBlock* getBlock( const real_t x, const real_t y, const real_t z ) const = 0;
           const IBlock* getBlock( const Vector3< real_t > & p ) const { return getBlock( p[0], p[1], p[2] ); }

   //*******************************************************************************************************************
   /*!
   *   \brief Returns the block located at position (x,y,z) (returns 'NULL' if the block doesn't exist locally!). The
   *   behavior of this function is identical to:
   *
   *   \code
   *     for( iterator block = begin(); block != end(); ++block )
   *        if( block->getAABB().contains(x,y,z) ) return *block;
   *     return NULL;
   *   \endcode
   *
   *   Periodicity is not considered! For mapping points to the periodic simulation space see 'mapToPeriodicDomain'.
   */
   //*******************************************************************************************************************
   virtual       IBlock* getBlock( const real_t x, const real_t y, const real_t z )       = 0;
                 IBlock* getBlock( const Vector3< real_t > & p ) { return getBlock( p[0], p[1], p[2] ); }



   /// Indicates whether or not information about remote blocks (blocks that reside on other processes) is available.
   /// This information includes the process rank, the state, and the axis-aligned bounding box of any block (local or remote).
   virtual bool containsGlobalBlockInformation() const = 0;

   /// Returns the block ID of every block in the simulation (global and remote).
   /// This member function is guaranteed to work properly only if 'containsGlobalBlockInformation() == true'.
   virtual void getAllBlocks( std::vector< shared_ptr< IBlockID > >& blocks ) const = 0;

   /// Returns true if there exists a block at position (x,y,z) [Periodicity is not considered! For mapping points to the
   /// periodic simulation space see 'mapToPeriodicDomain'.].
   /// This member function is guaranteed to work properly only if 'containsGlobalBlockInformation() == true'.
   virtual bool blockExists        ( const real_t x, const real_t y, const real_t z ) const = 0;
           bool blockExists        ( const Vector3< real_t > & p ) const { return blockExists( p[0], p[1], p[2] ); }

   /// Returns true if locally there exists a block at position (x,y,z) [Periodicity is not considered! For mapping points
   /// to the periodic simulation space see 'mapToPeriodicDomain'.].
   /// This member function is always guaranteed to work properly, even if 'containsGlobalBlockInformation() == false'.
   virtual bool blockExistsLocally ( const real_t x, const real_t y, const real_t z ) const = 0;
           bool blockExistsLocally ( const Vector3< real_t > & p ) const { return blockExistsLocally( p[0], p[1], p[2] ); }

   /// Returns true if remotely there exists a block at position (x,y,z) [Periodicity is not considered! For mapping points
   /// to the periodic simulation space see 'mapToPeriodicDomain'.].
   /// This member function is guaranteed to work properly only if 'containsGlobalBlockInformation() == true'.
   virtual bool blockExistsRemotely( const real_t x, const real_t y, const real_t z ) const = 0;
           bool blockExistsRemotely( const Vector3< real_t > & p ) const { return blockExistsRemotely( p[0], p[1], p[2] ); }

   /// Returns true if there exists a block with block ID 'id'.
   /// This member function is guaranteed to work properly only if 'containsGlobalBlockInformation() == true'.
   virtual bool blockExists        ( const IBlockID & id ) const = 0;

   /// Returns true if locally there exists a block with block ID 'id'.
   /// This member function is always guaranteed to work properly, even if 'containsGlobalBlockInformation() == false'.
   virtual bool blockExistsLocally ( const IBlockID & id ) const = 0;

   /// Returns true if remotely there exists a block with block ID 'id'.
   /// This member function is guaranteed to work properly only if 'containsGlobalBlockInformation() == true'.
   virtual bool blockExistsRemotely( const IBlockID & id ) const = 0;

   /// Returns the block ID that corresponds to the block located at position (x,y,z) [Periodicity is not considered!
   /// For mapping points to the periodic simulation space see 'mapToPeriodicDomain'.]. For local blocks, this function
   /// is always guaranteed to work. For remote blocks, this function is guaranteed to work properly only if
   /// 'containsGlobalBlockInformation() == true'.
   /// If the request cannot be satisfied (for example if no block exists at location (x,y,z)), the simulation must be
   /// aborted and the call to this function must not return!
   virtual void getBlockID( IBlockID & id, const real_t x, const real_t y, const real_t z ) const = 0;
           void getBlockID( IBlockID & id, const Vector3< real_t > & p ) const { getBlockID( id, p[0], p[1], p[2] ); }

   /// must return the level the block "block" is assigned to (must be an unsigned integer in the range [0,'number-of-levels') )
   virtual uint_t getLevel( const IBlock& /*block*/ ) const { return 0; }

   /// Returns the axis-aligned bounding box that corresponds to the block with ID 'id'. For local blocks, this function
   /// is always guaranteed to work. For remote blocks, this function is guaranteed to work properly only if
   /// 'containsGlobalBlockInformation() == true'.
   /// If the request cannot be satisfied, the simulation must be aborted and the call to this function must not return!
   virtual void getAABB( AABB & aabb, const IBlockID & id ) const = 0;

   inline AABB getAABB( const IBlockID & id ) const { AABB aabb; getAABB( aabb, id ); return aabb; }

   /// Returns the block state that corresponds to the block with ID 'id'. For local blocks, this function is always
   /// guaranteed to work. For remote blocks, this function is guaranteed to work properly only if
   /// 'containsGlobalBlockInformation() == true'.
   /// If the request cannot be satisfied, the simulation must be aborted and the call to this function must not return!
   virtual void getState( Set<SUID> & state, const IBlockID & id ) const = 0;

   inline Set<SUID> getState( const IBlockID & id ) const { Set<SUID> state; getState( state, id ); return state; }

   /// Returns the rank of the process the block with ID 'id' resides on. For local blocks, this function is always
   /// guaranteed to work. For remote blocks, this function is guaranteed to work properly only if
   /// 'containsGlobalBlockInformation() == true'.
   /// If the request cannot be satisfied, the simulation must be aborted and the call to this function must not return!
   virtual void getProcessRank( uint_t & rank, const IBlockID & id ) const = 0;

   inline uint_t getProcessRank( const IBlockID & id ) const { uint_t rank; getProcessRank( rank, id ); return rank; }



   /// Returns true if the block 'block' is located at the lower x-axis border of the domain
   inline bool atDomainXMinBorder( const IBlock & block ) const;
   /// Returns true if the block 'block' is located at the upper x-axis border of the domain
   inline bool atDomainXMaxBorder( const IBlock & block ) const;
   /// Returns true if the block 'block' is located at the lower y-axis border of the domain
   inline bool atDomainYMinBorder( const IBlock & block ) const;
   /// Returns true if the block 'block' is located at the upper y-axis border of the domain
   inline bool atDomainYMaxBorder( const IBlock & block ) const;
   /// Returns true if the block 'block' is located at the lower z-axis border of the domain
   inline bool atDomainZMinBorder( const IBlock & block ) const;
   /// Returns true if the block 'block' is located at the upper z-axis border of the domain
   inline bool atDomainZMaxBorder( const IBlock & block ) const;

   inline bool atDomainMinBorder( const uint_t index, const IBlock & block ) const;
   inline bool atDomainMaxBorder( const uint_t index, const IBlock & block ) const;



   /// Returns all neighboring process IDs
   virtual const std::vector< uint_t > & getNeighboringProcesses() const = 0;

   //*******************************************************************************************************************
   /*!
   *   For every neighboring process one or more offsets must be returned. Normally, for every neighboring process
   *   exactly one offset is returned: (0,0,0). Only in case of periodicity this offset can be different from (0,0,0)!
   *   If the offset is different from (0,0,0), the x-component must be identical to either 0 or +/- the x-width of the
   *   domain. The same applies for y and z.
   *   What this offset is supposed to do:
   *   Each block in 3D has 26 neighborhood regions. Suppose we leave block A in direction (-1,0,0) and enter the
   *   neighboring block B, which is stored on another process. If A is located at the x-min border of the domain and
   *   the simulation is only periodic in x-direction, block B must be located at the x-max border of the domain and the
   *   corresponding offset to block B - or more precisely to the process storing block B - is (+ x-width of the domain, 0, 0).
   *   Since each process can store multiple blocks and a simulation can be periodic in more than just one direction,
   *   there may be multiple offsets to another process - hence the std::vector.
   *   For an actual implementation of "getNeighboringProcessOffsets()" see BlockForest.cpp.
   */
   //*******************************************************************************************************************
   virtual std::map< uint_t, std::vector< Vector3<real_t> > > getNeighboringProcessOffsets() const = 0;



   /// Must be used if multiple block data handling objects with different selection attributes are registered for initializing
   /// the same block data "item". A 'BlockDataID' that corresponds to this block data "item" is returned.
   ///
   /// Usage: BlockDataID id = blockStorage.addBlockData( "[optional block data identifier]" ) << BlockDataCreator( ... )
   ///                                                                                         << BlockDataCreator( ... ) << ... ;
   internal::BlockDataHandlingAdder addBlockData( const std::string & identifier = std::string() ) { return internal::BlockDataHandlingAdder( *this, identifier ); }

   template< typename T >
   inline BlockDataID addBlockData( const shared_ptr< T > & dataHandling,
                                    const std::string & identifier          = std::string(),
                                    const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() );

   template< typename T >
   inline BlockDataID addBlockData( std::function< T* ( IBlock* const block ) > function,
                                    const std::string & identifier          = std::string(),
                                    const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                    const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() );
                                    
   BlockDataID addBlockData( const internal::SelectableBlockDataHandlingWrapper & dataHandling, const std::string & identifier = std::string() );

   template< typename T >
   inline BlockDataID loadBlockData( const std::string & file, const shared_ptr< T > & dataHandling,
                                     const std::string & identifier          = std::string(),
                                     const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                                     const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() );
                                    
   BlockDataID loadBlockData( const std::string & file,
                              const internal::SelectableBlockDataHandlingWrapper & dataHandling, const std::string & identifier = std::string() );
                              
   void saveBlockData( const std::string & file, const BlockDataID & id );

   void serializeBlockData( const BlockDataID & id, mpi::SendBuffer & buffer );
   void deserializeBlockData( const BlockDataID & id, mpi::RecvBuffer & buffer );

   inline void clearBlockData( const BlockDataID & id );

   uint_t numberOfBlockDataItems() const { return blockDataItem_.size(); }

   inline std::vector< std::string > getBlockDataIdentifiers() const;
   inline const std::string &        getBlockDataIdentifier( const ConstBlockDataID & id ) const;

   void rebuildProcessesWithBlocksCommunicator() { rebuildProcessesWithBlocksCommunicator_ = true; }

   MPI_Comm processesWithBlocksCommunicator();

protected:

   /// Every derived class must call this constructor!
   inline BlockStorage( const AABB& domain, const bool xPeriodic, const bool yPeriodic, const bool zPeriodic );

   virtual inline ~BlockStorage(); ///< Must not be made public! No one should be allowed to delete a variable of type 'BlockStorage*'

   virtual bool equal( const BlockStorage* rhs ) const = 0;

   inline void addBlockData( IBlock * const block, const BlockDataID & index, internal::BlockData * const data );



   AABB domain_;      ///< axis-aligned bounding box for the entire simulation space/domain
   bool periodic_[3]; ///< periodicity flags
   
   std::vector< internal::BlockDataItem > blockDataItem_;

private:

   BlockStorage(); ///< Must not be made public or protected! Derived classes must call one of the available public/protected constructors.

   inline void registerBlock( const std::pair<IBlockID::IDType, IBlock*>& block ); // All three functions must not be made public!
   inline void removeBlock  ( const IBlock* block );                               // All three functions are intended for internal use only.
   inline void removeBlock  ( const IBlockID::IDType& blockID );                   // All three functions are intended for internal use only.


   BlockContainerType iBlocks_; ///< holds pointers to all locally allocated blocks

   bool     rebuildProcessesWithBlocksCommunicator_;
   MPI_Comm processesWithBlocksCommunicator_; ///<  MPI communicator that only contains processes that possess blocks

}; // class BlockStorage



inline BlockStorage::BlockStorage( const AABB& domain, const bool xPeriodic, const bool yPeriodic, const bool zPeriodic ) :

   domain_( domain ), rebuildProcessesWithBlocksCommunicator_( true ), processesWithBlocksCommunicator_( MPI_COMM_NULL ) {

   periodic_[0] = xPeriodic;
   periodic_[1] = yPeriodic;
   periodic_[2] = zPeriodic;
}



inline BlockStorage::~BlockStorage()
{
   WALBERLA_MPI_SECTION()
   {
      if( processesWithBlocksCommunicator_ != MPI_COMM_NULL )
         MPI_Comm_free( &processesWithBlocksCommunicator_ );
   }
}



/// The following members are not used for checking if two BlockStorage objects are equal: iBlocks_
inline bool BlockStorage::operator==( const BlockStorage& rhs ) const {

   if( domain_ != rhs.domain_ || periodic_[0] != rhs.periodic_[0] || periodic_[1] != rhs.periodic_[1] || periodic_[2] != rhs.periodic_[2] )
      return false;

   return blockDataItem_ == rhs.blockDataItem_ && equal( &rhs );
}



inline bool BlockStorage::atDomainXMinBorder( const IBlock & block ) const
{
   WALBERLA_ASSERT_EQUAL( this, &block.getBlockStorage() );
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.xMin(), domain_.xMin(), real_c( 1.0E-6 ) * ( blockAABB.xMax() - blockAABB.xMin() ) );
}



inline bool BlockStorage::atDomainXMaxBorder( const IBlock & block ) const
{
   WALBERLA_ASSERT_EQUAL( this, &block.getBlockStorage() );
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.xMax(), domain_.xMax(), real_c( 1.0E-6 ) * ( blockAABB.xMax() - blockAABB.xMin() ) );
}



inline bool BlockStorage::atDomainYMinBorder( const IBlock & block ) const
{
   WALBERLA_ASSERT_EQUAL( this, &block.getBlockStorage() );
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.yMin(), domain_.yMin(), real_c( 1.0E-6 ) * ( blockAABB.yMax() - blockAABB.yMin() ) );
}



inline bool BlockStorage::atDomainYMaxBorder( const IBlock & block ) const
{
   WALBERLA_ASSERT_EQUAL( this, &block.getBlockStorage() );
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.yMax(), domain_.yMax(), real_c( 1.0E-6 ) * ( blockAABB.yMax() - blockAABB.yMin() ) );
}



inline bool BlockStorage::atDomainZMinBorder( const IBlock & block ) const
{
   WALBERLA_ASSERT_EQUAL( this, &block.getBlockStorage() );
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.zMin(), domain_.zMin(), real_c( 1.0E-6 ) * ( blockAABB.zMax() - blockAABB.zMin() ) );
}



inline bool BlockStorage::atDomainZMaxBorder( const IBlock & block ) const
{
   WALBERLA_ASSERT_EQUAL( this, &block.getBlockStorage() );
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.zMax(), domain_.zMax(), real_c( 1.0E-6 ) * ( blockAABB.zMax() - blockAABB.zMin() ) );
}



inline bool BlockStorage::atDomainMinBorder( const uint_t index, const IBlock & block ) const
{
   WALBERLA_ASSERT_LESS( index, uint_t(3) );
   WALBERLA_ASSERT_EQUAL( this, &block.getBlockStorage() );
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.min( index ), domain_.min( index ), real_c( 1.0E-6 ) * ( blockAABB.max( index ) - blockAABB.min( index) ) );
}



inline bool BlockStorage::atDomainMaxBorder( const uint_t index, const IBlock & block ) const
{
   WALBERLA_ASSERT_LESS( index, uint_t(3) );
   WALBERLA_ASSERT_EQUAL( this, &block.getBlockStorage() );
   const AABB & blockAABB = block.getAABB();
   return realIsEqual( blockAABB.max( index ), domain_.max( index ), real_c( 1.0E-6 ) * ( blockAABB.max( index ) - blockAABB.min( index ) ) );
}



//**********************************************************************************************************************
/*!
*   This function can be used for initializing a new block data "item". A 'BlockDataID' that corresponds to this block
*   data "item" is returned. This block data ID must be used to access/retrieve block data (see function 'getData' of
*   class 'IBlock').
*   If multiple block data handling objects with different selection attributes shall be used for initializing the same
*   block data "item" see member function "addBlockData( const std::string & identifier )"
*/
//**********************************************************************************************************************
template< typename T >
inline BlockDataID BlockStorage::addBlockData( const shared_ptr< T > & dataHandling, const std::string & identifier,
                                               const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors )
{
   //static_assert( std::is_base_of< BlockDataHandling<typename T::value_type>, T >::value );
   internal::SelectableBlockDataHandlingWrapper sbdhw( walberla::make_shared< internal::BlockDataHandlingHelper<typename T::value_type> >( dataHandling ),
                                                       requiredSelectors, incompatibleSelectors, identifier );

   return addBlockData( sbdhw, identifier );
}



//**********************************************************************************************************************
/*!
*   This function can be used for initializing a new block data "item". A 'BlockDataID' that corresponds to this block
*   data "item" is returned. This block data ID must be used to access/retrieve block data (see function 'getData' of
*   class 'IBlock').
*   If multiple initialization functions with different selection attributes shall be used for initializing the same
*   block data "item" see member function "addBlockData( const std::string & identifier )"
*/
//**********************************************************************************************************************
template< typename T >
inline BlockDataID BlockStorage::addBlockData( std::function< T* ( IBlock* const block ) > function, const std::string & identifier,
                                               const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors )
{
   internal::SelectableBlockDataHandlingWrapper dataHandling(
            walberla::make_shared< internal::BlockDataHandlingHelper<T> >( walberla::make_shared< internal::BlockDataHandlingFunctionAdaptor<T> >( function ) ),
            requiredSelectors, incompatibleSelectors, identifier );

   return addBlockData( dataHandling, identifier );
}



//**********************************************************************************************************************
/*!
*   This function can be used for initializing a new block data "item" by loading the data from file.
*   A 'BlockDataID' that corresponds to this block data "item" is returned. This block data ID must be used to
*   access/retrieve block data (see function 'getData' of class 'IBlock').
*/
//**********************************************************************************************************************
template< typename T >
inline BlockDataID BlockStorage::loadBlockData( const std::string & file, const shared_ptr< T > & dataHandling, const std::string & identifier,
                                                const Set<SUID> & requiredSelectors, const Set<SUID> & incompatibleSelectors )
{
   internal::SelectableBlockDataHandlingWrapper sbdhw( walberla::make_shared< internal::BlockDataHandlingHelper<typename T::value_type> >( dataHandling ),
                                                       requiredSelectors, incompatibleSelectors, identifier );

   return loadBlockData( file, sbdhw, identifier );
}



//**********************************************************************************************************************
/*!
*   This function can be used for removing all data that corresponds to block data ID 'id'.
*   Please note: The block data ID 'id' will still be valid, but blocks won't return anything anymore,
*                they will only return NULL for 'id'.
*/
//**********************************************************************************************************************
inline void BlockStorage::clearBlockData( const BlockDataID & id )
{
   for( auto block = begin(); block != end(); ++block )
      block->deleteData( id );
}



inline std::vector< std::string > BlockStorage::getBlockDataIdentifiers() const
{
   std::vector< std::string > identifiers;

   for( auto it = blockDataItem_.begin(); it != blockDataItem_.end(); ++it )
      identifiers.push_back( it->getIdentifier() );

   return identifiers;
}



inline const std::string & BlockStorage::getBlockDataIdentifier( const ConstBlockDataID & id ) const
{
   static std::string noData( "[no block data]" );

   if( !(id < blockDataItem_.size()) )
      return noData;

   return blockDataItem_[id].getIdentifier();
}



inline void BlockStorage::addBlockData( IBlock * const block, const BlockDataID & index, internal::BlockData * const data )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   block->addData( index, data );
}



//**********************************************************************************************************************
/*!
*   This function is called every time a new local block that belongs to this block storage data structure is created
*   (see the constructor of 'IBlock'). By doing so, the base class 'BlockStorage' is able to provide iterators for
*   traversing all locally allocated blocks.
*/
//**********************************************************************************************************************
inline void BlockStorage::registerBlock( const std::pair<IBlockID::IDType, IBlock*>& block ) {

   WALBERLA_ASSERT_NOT_NULLPTR( block.second );

   iBlocks_.insert(block);
}



//**********************************************************************************************************************
/*!
*   This function is called every time a block that belongs to this block storage data structure is destroyed (see the
*   destructor of 'IBlock'). By doing so, the base class 'BlockStorage' is able to provide iterators for traversing all
*   locally allocated blocks.
*/
//**********************************************************************************************************************
inline void BlockStorage::removeBlock( const IBlock* block  )
{
   auto it = iBlocks_.end();
   for (auto blkIt = iBlocks_.begin(); blkIt != iBlocks_.end(); ++blkIt)
   {
      if (blkIt->second == block)
      {
         it = blkIt;
         break;
      }
   }
   WALBERLA_ASSERT_UNEQUAL( it, iBlocks_.end() );
   removeBlock( it->first );
}



//**********************************************************************************************************************
/*!
*   This function is called every time a block that belongs to this block storage data structure is destroyed (see the
*   destructor of 'IBlock'). By doing so, the base class 'BlockStorage' is able to provide iterators for traversing all
*   locally allocated blocks.
*/
//**********************************************************************************************************************
inline void BlockStorage::removeBlock( const IBlockID::IDType& blockID  )
{
   iBlocks_.erase( blockID );
}



} // namespace domain_decomposition

using domain_decomposition::BlockDataCreator;
using domain_decomposition::BlockStorage;

} // namespace walberla
