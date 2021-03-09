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
//! \file BoundaryHandling.h
//! \ingroup boundary
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "Boundary.h"

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/cell/CellSet.h"
#include "core/cell/CellVector.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/uid/UID.h"
#include "core/uid/UIDGenerators.h"

#include "domain_decomposition/IBlock.h"

#include "field/FlagField.h"

#include "stencil/Directions.h"

#include <tuple>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>


namespace walberla {
namespace boundary {

namespace internal {
#if defined(__GLIBCXX__) && (!defined(_GLIBCXX_RELEASE) || _GLIBCXX_RELEASE < 7)
template< typename T >
struct tuple_size : std::tuple_size<std::tuple<>>
{};

template< typename... Members >
struct tuple_size< std::tuple<Members...> > : std::tuple_size<std::tuple<Members...>>
{};
#else
using std::tuple_size;
#endif
}



class BHUIDGenerator : public uid::IndexGenerator< BHUIDGenerator, uint_t >{};
using BoundaryHandlingUID = UID<BHUIDGenerator>;



template< typename FlagField_T, typename Stencil, typename... Boundaries > // Boundaries: all the boundary classes that are considered by this boundary handler
class BoundaryHandling
{
public:

   template< typename F, typename... T >
   friend class BoundaryHandlingCollection;

   using FlagField = FlagField_T;
   using flag_t = typename FlagField_T::flag_t;
   using ConstFlagFieldBaseIterator = typename FlagField_T::const_base_iterator;

   enum Mode { OPTIMIZED_SPARSE_TRAVERSAL, ENTIRE_FIELD_TRAVERSAL };



   class BlockSweep {
   public:
      BlockSweep( const BlockDataID & handling, const uint_t numberOfGhostLayersToInclude = 0 ) :
         handling_( handling ), numberOfGhostLayersToInclude_( numberOfGhostLayersToInclude ) {}
      void operator()( IBlock * block ) const {
         BoundaryHandling * handling = block->getData< BoundaryHandling >( handling_ );
         (*handling)( numberOfGhostLayersToInclude_ );
      }
   protected:
      const BlockDataID handling_;
      const uint_t numberOfGhostLayersToInclude_;
   };



   BoundaryHandling( const std::string & identifier, FlagField_T * const flagField, const flag_t domain, const Boundaries & ... boundaryConditions,
                     const Mode mode );
   BoundaryHandling( const std::string & identifier, FlagField_T * const flagField, const flag_t domain, const Boundaries & ... boundaryConditions ) :
      BoundaryHandling( identifier, flagField, domain, boundaryConditions..., OPTIMIZED_SPARSE_TRAVERSAL )
   {}

   bool operator==( const BoundaryHandling & rhs ) const { WALBERLA_CHECK( false, "You are trying to compare boundary handling " << uid_ <<                        // For testing purposes, block data items must be comparable with operator "==".
                                                                                  " with boundary handling " << rhs.getUID() <<                                    // Since instances of type "BoundaryHandling" are registered as block data items,
                                                                                  ".\nHowever, boundary handling instances are not comparable!" ); return false; } // "BoundaryHandling" must implement operator "==". As of right now, comparing
                                                                                                                                                                   // two boundary handling instances will always fail... :-) TODO: fixit?
   bool operator!=( const BoundaryHandling & rhs ) const { return !operator==( rhs ); }

   const BoundaryHandlingUID & getUID() const { return uid_; }

   const FlagField_T * getFlagField() const { return flagField_; }
         FlagField_T * getFlagField()       { return flagField_; }

   flag_t getNearBoundaryFlag() const { return nearBoundary_; } ///< Never (ever) set near boundary flags manually outside the boundary handler!
   flag_t getBoundaryMask() const { return boundary_; }
   flag_t getDomainMask() const { return domain_; }

   inline bool isEmpty       ( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline bool isNearBoundary( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline bool isBoundary    ( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline bool isDomain      ( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;

   inline bool isEmpty       ( const Cell & cell ) const;
   inline bool isNearBoundary( const Cell & cell ) const;
   inline bool isBoundary    ( const Cell & cell ) const;
   inline bool isDomain      ( const Cell & cell ) const;

   inline bool isEmpty       ( const ConstFlagFieldBaseIterator & it ) const;
   inline bool isNearBoundary( const ConstFlagFieldBaseIterator & it ) const;
   inline bool isBoundary    ( const ConstFlagFieldBaseIterator & it ) const;
   inline bool isDomain      ( const ConstFlagFieldBaseIterator & it ) const;

   inline bool containsBoundaryCondition( const BoundaryUID & uid ) const;
   inline bool containsBoundaryCondition( const FlagUID & flag ) const;
   inline bool containsBoundaryCondition( const flag_t    flag ) const { return ( boundary_ & flag ) == flag; }

   inline flag_t getBoundaryMask( const BoundaryUID & uid ) const { return getBoundaryMask( boundaryConditions_, uid ); }

   template< typename Boundary_T >
   inline const Boundary_T & getBoundaryCondition( const BoundaryUID & uid ) const; ///< You most likely have to call this function via "handling.template getBoundaryCondition< Boundary_T >( uid )".
   template< typename Boundary_T >
   inline       Boundary_T & getBoundaryCondition( const BoundaryUID & uid );       ///< You most likely have to call this function via "handling.template getBoundaryCondition< Boundary_T >( uid )".

   inline BoundaryUID getBoundaryUID( const FlagUID & flag ) const;
   inline BoundaryUID getBoundaryUID( const flag_t    flag ) const;

   inline uint_t numberOfMatchingBoundaryConditions( const flag_t mask ) const;

   inline bool checkConsistency( const uint_t numberOfGhostLayersToInclude = 0 ) const;
          bool checkConsistency( const CellInterval & cells ) const;

   inline void refresh( const uint_t numberOfGhostLayersToInclude = 0 ); // reset near ...
          void refresh( const CellInterval & cells );                    // ... boundary flags

   void refreshOutermostLayer( cell_idx_t thickness = 1 ); // reset near boundary flags in the outermost "inner" layers

   //** Set Domain Cells ***********************************************************************************************
   /*! \name Set Domain Cells */
   //@{
   inline void setDomain(                             const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   inline void setDomain( const flag_t domainSubFlag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   inline void setDomain(                             const CellInterval & cells );
   inline void setDomain( const flag_t domainSubFlag, const CellInterval & cells );

   template< typename CellIterator >
   inline void setDomain(                             const CellIterator & begin, const CellIterator & end );
   template< typename CellIterator >
   inline void setDomain( const flag_t domainSubFlag, const CellIterator & begin, const CellIterator & end );

   inline void forceDomain(                             const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   inline void forceDomain( const flag_t domainSubFlag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   inline void forceDomain(                             const CellInterval & cells );
   inline void forceDomain( const flag_t domainSubFlag, const CellInterval & cells );

   template< typename CellIterator >
   inline void forceDomain(                             const CellIterator & begin, const CellIterator & end );
   template< typename CellIterator >
   inline void forceDomain( const flag_t domainSubFlag, const CellIterator & begin, const CellIterator & end );

   inline void fillWithDomain(                             const uint_t numberOfGhostLayersToInclude = 0 );
   inline void fillWithDomain( const flag_t domainSubFlag, const uint_t numberOfGhostLayersToInclude     );

   inline void fillWithDomain(                             const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   inline void fillWithDomain( const flag_t domainSubFlag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   inline void fillWithDomain(                             const CellInterval & cells );
   inline void fillWithDomain( const flag_t domainSubFlag, const CellInterval & cells );

   template< typename CellIterator >
   inline void fillWithDomain(                             const CellIterator & begin, const CellIterator & end );
   template< typename CellIterator >
   inline void fillWithDomain( const flag_t domainSubFlag, const CellIterator & begin, const CellIterator & end );
   //@}
   //*******************************************************************************************************************

   //** Set Boundary Cells *********************************************************************************************
   /*! \name Set Boundary Cells */
   //@{
   inline shared_ptr<BoundaryConfiguration> createBoundaryConfiguration( const BoundaryUID & uid, const Config::BlockHandle & config ) const;

   inline void setBoundary( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                            const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   inline void setBoundary( const flag_t   flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                            const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   inline void setBoundary( const FlagUID & flag, const CellInterval & cells, const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   inline void setBoundary( const flag_t    flag, const CellInterval & cells, const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   template< typename CellIterator >
   inline void setBoundary( const FlagUID & flag, const CellIterator & begin, const CellIterator & end,
                            const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   template< typename CellIterator >
   inline void setBoundary( const flag_t    flag, const CellIterator & begin, const CellIterator & end,
                            const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   inline void forceBoundary( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                              const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   inline void forceBoundary( const flag_t    flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                              const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   inline void forceBoundary( const FlagUID & flag, const CellInterval & cells, const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   inline void forceBoundary( const flag_t    flag, const CellInterval & cells, const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   template< typename CellIterator >
   inline void forceBoundary( const FlagUID & flag, const CellIterator & begin, const CellIterator & end,
                              const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   template< typename CellIterator >
   inline void forceBoundary( const flag_t    flag, const CellIterator & begin, const CellIterator & end,
                              const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   //@}
   //*******************************************************************************************************************

   //** Remove Domain Cells ********************************************************************************************
   /*! \name Remove Domain Cells */
   //@{
   inline void removeDomain(                    const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   inline void removeDomain( const flag_t mask, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   inline void removeDomain( const uint_t numberOfGhostLayersToInclude = 0 );
   inline void removeDomain(                    const CellInterval & cells );
   inline void removeDomain( const flag_t mask, const CellInterval & cells );

   template< typename CellIterator >
   inline void removeDomain(                    const CellIterator & begin, const CellIterator & end );
   template< typename CellIterator >
   inline void removeDomain( const flag_t mask, const CellIterator & begin, const CellIterator & end );
   //@}
   //*******************************************************************************************************************

   //** Remove Boundary Cells ******************************************************************************************
   /*! \name Remove Boundary Cells */
   //@{
   inline void removeBoundary(                       const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   inline void removeBoundary( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   inline void removeBoundary( const flag_t    mask, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   inline void removeBoundary( const uint_t numberOfGhostLayersToInclude = 0 );
   inline void removeBoundary(                       const CellInterval & cells );
   inline void removeBoundary( const FlagUID & flag, const CellInterval & cells );
   inline void removeBoundary( const flag_t    mask, const CellInterval & cells );

   template< typename CellIterator >
   inline void removeBoundary(                       const CellIterator & begin, const CellIterator & end );
   template< typename CellIterator >
   inline void removeBoundary( const FlagUID & flag, const CellIterator & begin, const CellIterator & end );
   template< typename CellIterator >
   inline void removeBoundary( const flag_t    mask, const CellIterator & begin, const CellIterator & end );
   //@}
   //*******************************************************************************************************************

   //** General Flag Handling ******************************************************************************************
   /*! \name General Flag Handling */
   //@{
   inline void setFlag( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                        const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   inline void setFlag( const flag_t    flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                        const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   inline void setFlag( const FlagUID & flag, const CellInterval & cells,
                        const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   inline void setFlag( const flag_t    flag, const CellInterval & cells,
                        const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   template< typename CellIterator >
   inline void setFlag( const FlagUID & flag, const CellIterator & begin, const CellIterator & end,
                        const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   template< typename CellIterator >
   inline void setFlag( const flag_t    flag, const CellIterator & begin, const CellIterator & end,
                        const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   inline void forceFlag( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                          const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   inline void forceFlag( const flag_t    flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                          const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   inline void forceFlag( const FlagUID & flag, const CellInterval & cells,
                          const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   inline void forceFlag( const flag_t    flag, const CellInterval & cells,
                          const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   template< typename CellIterator >
   inline void forceFlag( const FlagUID & flag, const CellIterator & begin, const CellIterator & end,
                          const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   template< typename CellIterator >
   inline void forceFlag( const flag_t    flag, const CellIterator & begin, const CellIterator & end,
                          const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   inline void removeFlag( const FlagUID & flag, const uint_t numberOfGhostLayersToInclude = 0 );
   inline void removeFlag( const flag_t    flag, const uint_t numberOfGhostLayersToInclude = 0 );

   inline void removeFlag( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   inline void removeFlag( const flag_t    flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   inline void removeFlag( const FlagUID & flag, const CellInterval & cells );
   inline void removeFlag( const flag_t    flag, const CellInterval & cells );

   template< typename CellIterator >
   inline void removeFlag( const FlagUID & flag, const CellIterator & begin, const CellIterator & end );
   template< typename CellIterator >
   inline void removeFlag( const flag_t    flag, const CellIterator & begin, const CellIterator & end );
   //@}
   //*******************************************************************************************************************

   //** Clear Cells ****************************************************************************************************
   /*! \name Clear Cells */
   //@{
   inline void clear( const uint_t numberOfGhostLayersToInclude = 0 );
   inline void clear( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   inline void clear( const CellInterval & cells );
   template< typename CellIterator >
   inline void clear( const CellIterator & begin, const CellIterator & end );
   //@}
   //*******************************************************************************************************************

   //** Boundary Treatment *********************************************************************************************
   /*! \name Boundary Treatment */
   //@{
   static BlockSweep getBlockSweep( const BlockDataID handling, const uint_t numberOfGhostLayersToInclude = 0 )
      { return BlockSweep( handling, numberOfGhostLayersToInclude ); }

   inline void operator()( const uint_t numberOfGhostLayersToInclude = 0 );
   inline void operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   inline void operator()( const CellInterval & cells );
   template< typename CellIterator >
   inline void operator()( const CellIterator & begin, const CellIterator & end );

   inline void beforeBoundaryTreatment() { beforeBoundaryTreatment( boundaryConditions_ ); }
   inline void  afterBoundaryTreatment() {  afterBoundaryTreatment( boundaryConditions_ ); }
   //@}
   //*******************************************************************************************************************

   //** Pack / Unpack boundary handling ********************************************************************************
   /*! \name Pack / Unpack boundary handling */
   //@{
   template< typename Buffer_T >
   void   pack( Buffer_T & buffer, const CellInterval & interval, const bool assumeIdenticalFlagMapping = true ) const;
   template< typename Buffer_T >
   void unpack( Buffer_T & buffer, const CellInterval & interval, const bool assumeIdenticalFlagMapping = true );

   template< typename Buffer_T >
   void   pack( Buffer_T & buffer, stencil::Direction direction, const uint_t numberOfLayers = 1, const bool assumeIdenticalFlagMapping = true ) const;
   template< typename Buffer_T >
   void unpack( Buffer_T & buffer, stencil::Direction direction, const uint_t numberOfLayers = 1, const bool assumeIdenticalFlagMapping = true );
   //@}
   //*******************************************************************************************************************

   inline void        toStream( std::ostream & os ) const;
   inline std::string toString() const;

private:

   CellInterval getGhostLayerCellInterval( const uint_t numberOfGhostLayersToInclude ) const;

   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   typename std::enable_if<(N!=-1), void>::type setupBoundaryConditions(       BoundariesTuple & boundaryConditions );
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1  >
   typename std::enable_if<(N==-1), void>::type setupBoundaryConditions( const BoundariesTuple & ) const {}

   inline std::vector< BoundaryUID > getBoundaryUIDs() const;
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type getBoundaryUIDs( const BoundariesTuple & boundaryConditions, std::vector< BoundaryUID > & uids ) const;
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type getBoundaryUIDs( const BoundariesTuple &, std::vector< BoundaryUID > & ) const {}

   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), BoundaryUID>::type getBoundaryUID( const BoundariesTuple & boundaryConditions, const flag_t flag ) const;
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), BoundaryUID>::type getBoundaryUID( const BoundariesTuple &, const flag_t flagUID ) const;

   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), bool>::type containsBoundaryCondition( const BoundariesTuple & boundaryConditions, const BoundaryUID & uid ) const;
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), bool>::type containsBoundaryCondition( const BoundariesTuple &, const BoundaryUID & ) const { return false; }

   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), flag_t>::type getBoundaryMask( const BoundariesTuple & boundaryConditions, const BoundaryUID & uid ) const;
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), flag_t>::type getBoundaryMask( const BoundariesTuple &, const BoundaryUID & ) const { return numeric_cast<flag_t>(0); }

   //** Get Boundary Class (private helper functions) ******************************************************************
   /*! \name Get Boundary Class (private helper functions) */
   //@{

   // matching type (-> Boundary_T) not yet found ...

   template< typename Boundary_T, typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline const typename std::enable_if<(N!=0), Boundary_T>::type & getBoundaryCondition( const BoundaryUID & uid, const BoundariesTuple & boundaryConditions,
                                                   typename std::enable_if< std::is_same< Boundary_T, typename std::tuple_element<N, BoundariesTuple>::type >::value >::type* /*dummy*/ = nullptr ) const
   {
      if( uid == std::get<N>( boundaryConditions ).getUID() )
         return std::get<N>( boundaryConditions );
      else
         return getBoundaryCondition_TypeExists< Boundary_T, BoundariesTuple, N-1 >( uid, boundaryConditions );
   }

   template< typename Boundary_T, typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline const typename std::enable_if<(N==0), Boundary_T>::type & getBoundaryCondition( const BoundaryUID & uid, const BoundariesTuple & boundaryConditions,
                                                   typename std::enable_if< std::is_same< Boundary_T, typename std::tuple_element<N, BoundariesTuple>::type >::value >::type* /*dummy*/ = nullptr ) const
   {
      if( uid == std::get<N>( boundaryConditions ).getUID() )
         return std::get<N>( boundaryConditions );
      else
         WALBERLA_ABORT( "The requested boundary condition " << uid.getIdentifier() << " is not part of this boundary handling." );

#ifdef __IBMCPP__
      return *(reinterpret_cast< Boundary_T * >( NULL )); // silencing incorrect IBM compiler warning
#endif
   }

   template< typename Boundary_T, typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline const typename std::enable_if<(N!=0), Boundary_T>::type & getBoundaryCondition( const BoundaryUID & uid, const BoundariesTuple & boundaryConditions,
                                                   typename std::enable_if< std::is_same< typename std::is_same< Boundary_T, typename std::tuple_element<N, BoundariesTuple>::type >::type,
                                                                                              std::false_type >::value >::type* /*dummy*/ = nullptr,
                                                   typename std::enable_if< (N>0) >::type* /*dummy*/ = nullptr ) const
   {
      return getBoundaryCondition< Boundary_T, BoundariesTuple, N-1 >( uid, boundaryConditions );
   }

   template< typename Boundary_T, typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline const typename std::enable_if<(N==0), Boundary_T>::type & getBoundaryCondition( const BoundaryUID & /*uid*/, const BoundariesTuple & /*boundaryConditions*/,
                                                   typename std::enable_if< std::is_same< typename std::is_same< Boundary_T, typename std::tuple_element<0, BoundariesTuple>::type >::type,
                                                                                              std::false_type >::value >::type* /*dummy*/ = 0 ) const
   {
      static_assert( sizeof(Boundary_T) == 0, "The requested boundary class is not part of this boundary handling." );
   }

   // matching type (-> Boundary_T) exists!

   template< typename Boundary_T, typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline const typename std::enable_if<(N!=0), Boundary_T>::type & getBoundaryCondition_TypeExists( const BoundaryUID & uid, const BoundariesTuple & boundaryConditions,
                                                              typename std::enable_if< std::is_same< Boundary_T, typename std::tuple_element<N, BoundariesTuple>::type >::value >::type* /*dummy*/ = 0 ) const
   {
      if( uid == std::get<N>( boundaryConditions ).getUID() )
         return std::get<N>( boundaryConditions );
      else
         return getBoundaryCondition_TypeExists< Boundary_T, BoundariesTuple, N-1 >( uid, boundaryConditions );
   }

   template< typename Boundary_T, typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline const typename std::enable_if<(N==0), Boundary_T>::type & getBoundaryCondition_TypeExists( const BoundaryUID & uid, const BoundariesTuple & boundaryConditions,
                                                              typename std::enable_if< std::is_same< Boundary_T, typename std::tuple_element<0, BoundariesTuple>::type >::value >::type* /*dummy*/ = 0 ) const
   {
      if( uid == std::get<0>( boundaryConditions ).getUID() )
         return std::get<0>( boundaryConditions );
      else
         WALBERLA_ABORT( "The requested boundary condition " << uid.getIdentifier() << " is not part of this boundary handling." );

#ifdef __IBMCPP__
      return *(reinterpret_cast< Boundary_T * >( NULL )); // silencing incorrect IBM compiler warning
#endif
   }

   template< typename Boundary_T, typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline const typename std::enable_if<(N!=0), Boundary_T>::type & getBoundaryCondition_TypeExists( const BoundaryUID & uid, const BoundariesTuple & boundaryConditions,
                                                              typename std::enable_if< std::is_same< typename std::is_same< Boundary_T, typename std::tuple_element<N, BoundariesTuple>::type >::type,
                                                                                                         std::false_type >::value >::type* /*dummy*/ = 0 ) const
   {
      return getBoundaryCondition_TypeExists< Boundary_T, BoundariesTuple, N-1 >( uid, boundaryConditions );
   }

   template< typename Boundary_T, typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline const typename std::enable_if<(N==0), Boundary_T>::type & getBoundaryCondition_TypeExists( const BoundaryUID & uid, const BoundariesTuple & /*boundaryConditions*/,
                                                              typename std::enable_if< std::is_same< typename std::is_same< Boundary_T, typename std::tuple_element<0, BoundariesTuple>::type >::type,
                                                                                                         std::false_type >::value >::type* /*dummy*/ = nullptr ) const
   {
      WALBERLA_ABORT( "The requested boundary condition " << uid.getIdentifier() << " is not part of this boundary handling." );

#ifdef __IBMCPP__
      return *(reinterpret_cast< Boundary_T * >( NULL )); // silencing incorrect IBM compiler warning
#endif
   }

   //*******************************************************************************************************************

   inline uint_t numberOfMatchingBoundaryConditions( const BoundaryUID & uid ) const;

   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), uint_t>::type numberOfMatchingBoundaryConditions( const BoundariesTuple & boundaryConditions, const BoundaryUID & uid ) const;
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), uint_t>::type numberOfMatchingBoundaryConditions( const BoundariesTuple &, const BoundaryUID & ) const { return uint_c(0); }

   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), uint_t>::type numberOfMatchingBoundaryConditions( const BoundariesTuple & boundaryConditions, const flag_t mask ) const;
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), uint_t>::type numberOfMatchingBoundaryConditions( const BoundariesTuple &, const flag_t ) const { return uint_c(0); }

   inline bool checkFlagField( const uint_t numberOfGhostLayersToInclude = 0 ) const;

   inline void addDomain( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const flag_t domain );

   //** Set Boundary Cells (private helper functions) ******************************************************************
   /*! \name Set Boundary Cells (private helper functions) */
   //@{
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), shared_ptr<BoundaryConfiguration>>::type createBoundaryConfiguration( const BoundariesTuple & boundaryConditions,
                                                                         const BoundaryUID & uid, const Config::BlockHandle & config ) const;
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), shared_ptr<BoundaryConfiguration>>::type createBoundaryConfiguration( const BoundariesTuple &, const BoundaryUID & uid,
                                                                         const Config::BlockHandle & ) const;

   inline void addNearBoundary( const CellInterval & cells );
   inline void addBoundary( const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   template< typename BoundariesTuple, int N = internal::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type setBoundary(       BoundariesTuple & boundaryConditions, const flag_t flag,
                                                                                        const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                                        const BoundaryConfiguration & parameter );
   template< typename BoundariesTuple, int N = internal::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type setBoundary( const BoundariesTuple &, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t,
                                                              const BoundaryConfiguration & ) const;

   template< typename BoundariesTuple, int N = internal::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type setBoundary(       BoundariesTuple & boundaryConditions, const flag_t flag, const CellInterval & cells,
                                                                                        const BoundaryConfiguration & parameter );
   template< typename BoundariesTuple, int N = internal::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type setBoundary( const BoundariesTuple &, const flag_t, const CellInterval &, const BoundaryConfiguration & ) const;

   template< typename BoundariesTuple, int N = internal::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type setBoundary(       BoundariesTuple & boundaryConditions, const flag_t flag, const CellVector & cells,
                                                                                        const BoundaryConfiguration & parameter );
   template< typename BoundariesTuple, int N = internal::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type setBoundary( const BoundariesTuple &, const flag_t, const CellVector &, const BoundaryConfiguration & ) const;
   //@}
   //*******************************************************************************************************************

   //** Remove Boundary Cells (private helper functions) ***************************************************************
   /*! \name Remove Boundary Cells (private helper functions) */
   //@{
   template< typename BoundariesTuple, int N = internal::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type removeBoundary(       BoundariesTuple & boundaryConditions, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                                           const bool checkNearBoundaryFlags = true );
   template< typename BoundariesTuple, int N = internal::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type removeBoundary( const BoundariesTuple &, const cell_idx_t, const cell_idx_t, const cell_idx_t, const bool ) const { WALBERLA_CHECK( false ); }
   //@}
   //*******************************************************************************************************************

   //** Boundary Treatment (private helper functions) ******************************************************************
   /*! \name Boundary Treatment (private helper functions) */
   //@{
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type treatDirection( BoundariesTuple & boundaryConditions, const uint_t index,
                               const std::vector< std::vector< std::pair< Cell, stencil::Direction > > > & cellDirectionPairs );
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type treatDirection( const BoundariesTuple &, const uint_t,
                               const std::vector< std::vector< std::pair< Cell, stencil::Direction > > > & ) const {}

   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type treatDirection(       BoundariesTuple & boundaryConditions, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                                           const stencil::Direction dir,
                                                                                           const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz );
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type treatDirection( const BoundariesTuple & , const cell_idx_t, const cell_idx_t, const cell_idx_t, const stencil::Direction,
                                                                  const cell_idx_t, const cell_idx_t, const cell_idx_t ) const { WALBERLA_CHECK( false ); }

   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type beforeBoundaryTreatment(       BoundariesTuple & boundaryConditions );
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type beforeBoundaryTreatment( const BoundariesTuple & ) const {}

   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type afterBoundaryTreatment(       BoundariesTuple & boundaryConditions );
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type afterBoundaryTreatment( const BoundariesTuple & ) const {}
   //@}
   //*******************************************************************************************************************

   //** Pack / Unpack boundary handling (private helper functions) *****************************************************
   /*! \name Pack / Unpack boundary handling (private helper functions) */
   //@{
   inline std::map< std::string, flag_t > getFlagMapping() const;

   template< typename Buffer_T >
   std::vector< flag_t > getNeighborFlagMapping( Buffer_T & buffer, const bool assumeIdenticalFlagMapping, bool & identicalFlagMapping ) const;

   void translateMask( flag_t & mask, const std::vector< flag_t > & flagMapping ) const;

   CellInterval   getPackingInterval( stencil::Direction direction, const uint_t numberOfLayers ) const;
   CellInterval getUnpackingInterval( stencil::Direction direction, const uint_t numberOfLayers ) const;

   template< typename Buffer_T >
   inline void pack( Buffer_T & buffer, const flag_t mask, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;

   template< typename Buffer_T, typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type pack( const BoundariesTuple & boundaryConditions, Buffer_T & buffer,
                     const flag_t mask, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   template< typename Buffer_T, typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   typename std::enable_if<(N==-1), void>::type pack( const BoundariesTuple &, Buffer_T &, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t )
              const { WALBERLA_CHECK( false ); }

   template< typename Buffer_T >
   inline void unpack( Buffer_T & buffer, const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   template< typename Buffer_T >
   inline void unpackBoundary( Buffer_T & buffer, const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   template< typename Buffer_T, typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type unpackBoundary( BoundariesTuple & boundaryConditions, Buffer_T & buffer, const flag_t flag,
                               const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   template< typename Buffer_T, typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type unpackBoundary( const BoundariesTuple &, Buffer_T &, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t)
                               const { WALBERLA_CHECK( false ); }
   //@}
   //*******************************************************************************************************************

   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   typename std::enable_if<(N!=-1), void>::type getBoundaryConditions( const BoundariesTuple & boundaryConditions, std::vector< std::string > & bcs ) const;
   template< typename BoundariesTuple, int N = std::tuple_size<BoundariesTuple>::value - 1 >
   typename std::enable_if<(N==-1), void>::type getBoundaryConditions( const BoundariesTuple &, std::vector< std::string > & ) const {}

   template< typename T > static void valueToStream( std::ostream & os, const T       value ) { os << value; }
                          static void valueToStream( std::ostream & os, const  int8_t value ) { os <<  int_c( value ); }
                          static void valueToStream( std::ostream & os, const uint8_t value ) { os << uint_c( value ); }



   const BoundaryHandlingUID uid_;

   FlagField_T * const flagField_;

   const CellInterval innerBB_;
   const CellInterval outerBB_;

   const flag_t nearBoundary_; // "near boundary" flags are only set in cells that are marked as "domain" cells and that are ...
         flag_t boundary_;     // ... located next to "boundary" cells. "boundary_" and "domain_" must be disjoint!
   const flag_t domain_;       // Cells that are marked as neither "boundary" cells nor "domain" cells are ignored.

   const Mode mode_;

   bool dirty_; //!< Every time "near boundary" flags are set or removed, dirty_ is set to true.
                /*!< Triggers a rebuild of vector cellDirectionPairs_ in "OPTIMIZED_SPARSE_TRAVERSAL" mode. */

   std::vector< flag_t > bcMaskMapping_; // boundary condition index (position in tuple) to associated mask mapping

   std::vector< bool > rebuildCellDirectionPairs_;
   std::vector< std::vector< std::vector< std::pair< Cell, stencil::Direction > > > > cellDirectionPairs_; // 1st vector: numberOfGhostLayersToInclude
                                                                                                           // 2nd vector: boundary condition index
                                                                                                           // 3rd vector: vector of cell<->direction pairs
   using Tuple = std::tuple<Boundaries...>;
   Tuple boundaryConditions_;
   bool  threadSafeBCs_;

}; // class BoundaryHandling



template< typename FlagField_T, typename Stencil, typename... Boundaries >
BoundaryHandling< FlagField_T, Stencil, Boundaries... >::BoundaryHandling( const std::string & identifier, FlagField_T * const flagField,
                                                                   const flag_t domain, const Boundaries & ... boundaryConditions, const Mode mode ) :

   uid_( identifier ),
   flagField_( flagField ),
   innerBB_( 1 - cell_idx_c( flagField_->nrOfGhostLayers() ), 1 - cell_idx_c( flagField_->nrOfGhostLayers() ), 1 - cell_idx_c( flagField_->nrOfGhostLayers() ),
             cell_idx_c( flagField_->xSize() + flagField_->nrOfGhostLayers() ) - 2, cell_idx_c( flagField_->ySize() + flagField_->nrOfGhostLayers() ) - 2,
             cell_idx_c( flagField_->zSize() + flagField_->nrOfGhostLayers() ) - 2 ),
   outerBB_( -cell_idx_c( flagField_->nrOfGhostLayers() ), -cell_idx_c( flagField_->nrOfGhostLayers() ), -cell_idx_c( flagField_->nrOfGhostLayers() ),
             cell_idx_c( flagField_->xSize() + flagField_->nrOfGhostLayers() ) - 1, cell_idx_c( flagField_->ySize() + flagField_->nrOfGhostLayers() ) - 1,
             cell_idx_c( flagField_->zSize() + flagField_->nrOfGhostLayers() ) - 1 ),
   nearBoundary_( flagField_->registerFlag( std::string( "near boundary (" ) + identifier + std::string( ")" ) ) ),
   boundary_( 0 ),
   domain_( domain ),
   mode_( mode ),
   dirty_( false ),
   boundaryConditions_( std::make_tuple(boundaryConditions...) ),
   threadSafeBCs_( true )
{
   setupBoundaryConditions( boundaryConditions_ );

   if( flagField_->nrOfGhostLayers() < 1 )
      WALBERLA_ABORT( "The flag field passed to the boundary handling \"" << identifier << "\" must contain at least one ghost layer!" );

   rebuildCellDirectionPairs_.resize( flagField_->nrOfGhostLayers(), false );
   cellDirectionPairs_.resize( flagField_->nrOfGhostLayers() );

   WALBERLA_CHECK_EQUAL( nearBoundary_ & boundary_, flag_t(0), "The near boundary flag must not be identical to a flag used for marking a boundary cell.\n"
                                                               "This check failed for boundary handling " << uid_ << "!" );
   WALBERLA_CHECK_EQUAL( nearBoundary_ & domain_,   flag_t(0), "The near boundary flag must not be identical to a flag used for marking a domain cell.\n"
                                                               "This check failed for boundary handling " << uid_ << "!" );
   WALBERLA_CHECK_EQUAL( boundary_ & domain_,       flag_t(0), "Flags used for marking domain cells must be different to flags used for marking boundary cells.\n"
                                                               "This check failed for boundary handling " << uid_ << "!" );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::isEmpty( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   return !flagField_->isPartOfMaskSet( x, y, z, boundary_ ) && !flagField_->isPartOfMaskSet( x, y, z, domain_ );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::isNearBoundary( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   return flagField_->isFlagSet( x, y, z, nearBoundary_ );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::isBoundary( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   return flagField_->isPartOfMaskSet( x, y, z, boundary_ );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::isDomain( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   return flagField_->isPartOfMaskSet( x, y, z, domain_ );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::isEmpty( const Cell & cell ) const
{
   return !flagField_->isPartOfMaskSet( cell.x(), cell.y(), cell.z(), boundary_ ) && 
          !flagField_->isPartOfMaskSet( cell.x(), cell.y(), cell.z(), domain_ );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::isNearBoundary( const Cell & cell ) const
{
   return flagField_->isFlagSet( cell.x(), cell.y(), cell.z(), nearBoundary_ );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::isBoundary( const Cell & cell ) const
{
   return flagField_->isPartOfMaskSet( cell.x(), cell.y(), cell.z(), boundary_ );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::isDomain( const Cell & cell ) const
{
   return flagField_->isPartOfMaskSet( cell.x(), cell.y(), cell.z(), domain_ );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::isEmpty( const ConstFlagFieldBaseIterator & it ) const
{
   WALBERLA_ASSERT_EQUAL( it.getField(), flagField_ );

   return !isPartOfMaskSet( it, boundary_ ) && !isPartOfMaskSet( it, domain_ );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::isNearBoundary( const ConstFlagFieldBaseIterator & it ) const
{
   WALBERLA_ASSERT_EQUAL( it.getField(), flagField_ );

   return isFlagSet( it, nearBoundary_ );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::isBoundary( const ConstFlagFieldBaseIterator & it ) const
{
   WALBERLA_ASSERT_EQUAL( it.getField(), flagField_ );

   return isPartOfMaskSet( it, boundary_ );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::isDomain( const ConstFlagFieldBaseIterator & it ) const
{
   WALBERLA_ASSERT_EQUAL( it.getField(), flagField_ );

   return isPartOfMaskSet( it, domain_ );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::containsBoundaryCondition( const BoundaryUID & uid ) const
{
   return containsBoundaryCondition( boundaryConditions_, uid );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::containsBoundaryCondition( const FlagUID & flag ) const
{
   if( !flagField_->flagExists( flag ) )
      return false;

   return containsBoundaryCondition( flagField_->getFlag( flag ) );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename Boundary_T >
inline const Boundary_T & BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getBoundaryCondition( const BoundaryUID & uid ) const
{
   return getBoundaryCondition< Boundary_T, std::tuple< Boundaries... > >( uid, boundaryConditions_ );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename Boundary_T >
inline Boundary_T & BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getBoundaryCondition( const BoundaryUID & uid )
{
   return const_cast< Boundary_T & >( static_cast< const BoundaryHandling * >( this )->template getBoundaryCondition< Boundary_T >( uid ) );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline BoundaryUID BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getBoundaryUID( const FlagUID & flag ) const
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   return getBoundaryUID( flagField_->getFlag( flag ) );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline BoundaryUID BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getBoundaryUID( const flag_t flag ) const
{
   WALBERLA_ASSERT( field::isFlag( flag ) );
   WALBERLA_ASSERT( flagField_->isRegistered( flag ) );

   return getBoundaryUID( boundaryConditions_, flag );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline uint_t BoundaryHandling< FlagField_T, Stencil, Boundaries... >::numberOfMatchingBoundaryConditions( const flag_t mask ) const
{
   return numberOfMatchingBoundaryConditions( boundaryConditions_, mask );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::checkConsistency( const uint_t numberOfGhostLayersToInclude ) const
{
   CellInterval cells = getGhostLayerCellInterval( numberOfGhostLayersToInclude );
   return checkConsistency( cells );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::checkConsistency( const CellInterval & cells ) const
{
   CellInterval localCells( innerBB_ );
   localCells.intersect( cells );

   for( auto cell = flagField_->beginSliceXYZ( localCells ); cell != flagField_->end(); ++cell )
   {
      if( isPartOfMaskSet( cell, boundary_ ) ) // if the cell is marked as a boundary cell ...
      {
         // ... near boundary & domain must not be set
         WALBERLA_CHECK( !isPartOfMaskSet( cell, domain_ ) && !isPartOfMaskSet( cell, nearBoundary_ ) );
         if( isPartOfMaskSet( cell, domain_ ) || isPartOfMaskSet( cell, nearBoundary_ ) )
            return false;

         // ... exactly one boundary condition must match
         WALBERLA_CHECK_EQUAL( numberOfMatchingBoundaryConditions( *cell ), uint_t(1) );
         if( numberOfMatchingBoundaryConditions( *cell ) != 1 )
            return false;

         // ... only one bit of the mask may be set
         WALBERLA_CHECK( field::isFlag( *cell & boundary_ ) );
         if( !field::isFlag( *cell & boundary_ ) )
            return false;
      }
      else if( isPartOfMaskSet( cell, domain_ ) ) // if the cell is marked as a domain cell ...
      {
         // ... it must not be marked as a boundary cell
         WALBERLA_CHECK( !isPartOfMaskSet( cell, boundary_ ) );
         if( isPartOfMaskSet( cell, boundary_ ) )
            return false;

         bool boundaryNeighbor( false );
         for( auto d = Stencil::beginNoCenter(); d != Stencil::end(); ++d )
         {
            const flag_t & nMask = cell.neighbor(*d);
            if( isPartOfMaskSet( nMask, boundary_ ) )
            {
               // if a neighbor is marked as a boundary cell, exactly one boundary condition must match
               WALBERLA_CHECK_EQUAL( numberOfMatchingBoundaryConditions( nMask ), uint_t(1) );
               if( numberOfMatchingBoundaryConditions( nMask ) != 1 )
                  return false;
               boundaryNeighbor = true;
            }
         }

         // if at least one neighbor is marked as a boundary cell, near boundary must be set - otherwise it must not be set
         WALBERLA_CHECK( ( boundaryNeighbor && isFlagSet( cell, nearBoundary_ ) ) || ( !boundaryNeighbor && !isFlagSet( cell, nearBoundary_ ) ) );
         if( ( boundaryNeighbor && !isFlagSet( cell, nearBoundary_ ) ) || ( !boundaryNeighbor && isFlagSet( cell, nearBoundary_ ) ) )
            return false;

         // ... only one bit of the mask may be set
         WALBERLA_CHECK( field::isFlag( *cell & domain_ ) );
         if( !field::isFlag( *cell & domain_ ) )
            return false;
      }
   }

   return true;
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::refresh( const uint_t numberOfGhostLayersToInclude )
{
   CellInterval cells = getGhostLayerCellInterval( numberOfGhostLayersToInclude );
   refresh( cells );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::refresh( const CellInterval & cells )
{
   CellInterval localCells( innerBB_ );
   localCells.intersect( cells );

   for( auto cell = flagField_->beginSliceXYZ( localCells ); cell != flagField_->end(); ++cell )
   {
      if( isPartOfMaskSet( cell, domain_ ) )
      {
         field::removeFlag( cell, nearBoundary_ );

         for( auto d = Stencil::beginNoCenter(); d != Stencil::end(); ++d )
         {
            if( isPartOfMaskSet( cell.neighbor(*d), boundary_ ) )
            {
               addFlag( cell, nearBoundary_ );
               break;
            }
         }
      }
   }

   dirty_ = true;
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::refreshOutermostLayer( cell_idx_t thickness )
{
   uint_t extent = std::min( std::min( innerBB_.xSize(), innerBB_.ySize() ), innerBB_.zSize() );
   WALBERLA_ASSERT_GREATER( extent, uint_t(0) );

   if( extent == 1 )
   {
      refresh();
      return;
   }

   WALBERLA_ASSERT_GREATER_EQUAL( thickness, cell_idx_t(1) );
   WALBERLA_ASSERT_LESS_EQUAL( thickness, cell_idx_c( extent ) / cell_idx_t(2) );

   const cell_idx_t one = numeric_cast<cell_idx_t>(1);
   thickness -= one;

   CellInterval xlow ( innerBB_.xMin(), innerBB_.yMin(), innerBB_.zMin(), innerBB_.xMin() + thickness, innerBB_.yMax(), innerBB_.zMax() );
   CellInterval xhigh( innerBB_.xMax() - thickness, innerBB_.yMin(), innerBB_.zMin(), innerBB_.xMax(), innerBB_.yMax(), innerBB_.zMax() );

   CellInterval ylow ( innerBB_.xMin() + thickness + one, innerBB_.yMin(), innerBB_.zMin(),
                       innerBB_.xMax() - thickness - one, innerBB_.yMin() + thickness, innerBB_.zMax() );
   CellInterval yhigh( innerBB_.xMin() + thickness + one, innerBB_.yMax() - thickness, innerBB_.zMin(),
                       innerBB_.xMax() - thickness - one, innerBB_.yMax(), innerBB_.zMax() );

   CellInterval zlow ( innerBB_.xMin() + thickness + one, innerBB_.yMin() + thickness + one, innerBB_.zMin(),
                       innerBB_.xMax() - thickness - one, innerBB_.yMax() - thickness - one, innerBB_.zMin() + thickness );
   CellInterval zhigh( innerBB_.xMin() + thickness + one, innerBB_.yMin() + thickness + one, innerBB_.zMax() - thickness,
                       innerBB_.xMax() - thickness - one, innerBB_.yMax() - thickness - one, innerBB_.zMax() );

   refresh( xlow  );
   refresh( xhigh );
   refresh( ylow  );
   refresh( yhigh );
   refresh( zlow  );
   refresh( zhigh );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setDomain( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT( field::isFlag( domain_ ) );
   setDomain( domain_, x, y, z );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setDomain( const flag_t domainSubFlag,
                                                                        const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT_EQUAL( domain_ & domainSubFlag, domainSubFlag );
   WALBERLA_ASSERT( field::isFlag( domainSubFlag ) );

   if( outerBB_.contains(x,y,z) )
      addDomain( x, y, z, domainSubFlag );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setDomain( const CellInterval & cells )
{
   WALBERLA_ASSERT( field::isFlag( domain_ ) );
   setDomain( domain_, cells );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setDomain( const flag_t domainSubFlag, const CellInterval & cells )
{
   WALBERLA_ASSERT_EQUAL( domain_ & domainSubFlag, domainSubFlag );
   WALBERLA_ASSERT( field::isFlag( domainSubFlag ) );

   CellInterval localCells( outerBB_ );
   localCells.intersect( cells );

   if( localCells.empty() )
      return;

   CellInterval innerHull( localCells.xMin()+1, localCells.yMin()+1, localCells.zMin()+1,
                           localCells.xMax()-1, localCells.yMax()-1, localCells.zMax()-1 );

   for( auto cell = flagField_->beginSliceXYZ( localCells ); cell != flagField_->end(); ++cell )
   {
      const cell_idx_t x = cell.x();
      const cell_idx_t y = cell.y();
      const cell_idx_t z = cell.z();

      WALBERLA_ASSERT( isEmpty(x,y,z) );

      addFlag( cell, domainSubFlag );

      // near boundary?

      if( !innerHull.contains(x,y,z) && innerBB_.contains(x,y,z) )
      {
         for( auto d = Stencil::beginNoCenter(); d != Stencil::end(); ++d )
         {
            if( isPartOfMaskSet( cell.neighbor(*d), boundary_ ) )
            {
               addFlag( cell, nearBoundary_ );
               dirty_ = true;
               break;
            }
         }
      }
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setDomain( const CellIterator & begin, const CellIterator & end )
{
   WALBERLA_ASSERT( field::isFlag( domain_ ) );
   setDomain( domain_, begin, end );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setDomain( const flag_t domainSubFlag, const CellIterator & begin, const CellIterator & end )
{
   WALBERLA_ASSERT_EQUAL( domain_ & domainSubFlag, domainSubFlag );
   WALBERLA_ASSERT( field::isFlag( domainSubFlag ) );

   for( auto cell = begin; cell != end; ++cell )
   {
      const cell_idx_t x = cell->x();
      const cell_idx_t y = cell->y();
      const cell_idx_t z = cell->z();

      if( outerBB_.contains(x,y,z) )
         addDomain( x, y, z, domainSubFlag );
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceDomain( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   forceDomain( domain_, x, y, z );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceDomain( const flag_t domainSubFlag,
                                                                          const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   clear( x, y, z );
   setDomain( domainSubFlag, x, y, z );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceDomain( const CellInterval & cells )
{
   forceDomain( domain_, cells );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceDomain( const flag_t domainSubFlag, const CellInterval & cells )
{
   clear( cells );
   setDomain( domainSubFlag, cells );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceDomain( const CellIterator & begin, const CellIterator & end )
{
   forceDomain( domain_, begin, end );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceDomain( const flag_t domainSubFlag, const CellIterator & begin, const CellIterator & end )
{
   clear( begin, end );
   setDomain( domainSubFlag, begin, end );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::fillWithDomain( const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_GREATER_EQUAL( flagField_->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT( field::isFlag( domain_ ) );
   fillWithDomain( domain_, numberOfGhostLayersToInclude );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::fillWithDomain( const flag_t domainSubFlag, const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_EQUAL( domain_ & domainSubFlag, domainSubFlag );
   WALBERLA_ASSERT( field::isFlag( domainSubFlag ) );
   WALBERLA_ASSERT_GREATER_EQUAL( flagField_->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   CellInterval cells = getGhostLayerCellInterval( numberOfGhostLayersToInclude );
   fillWithDomain( domainSubFlag, cells );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::fillWithDomain( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT( field::isFlag( domain_ ) );
   fillWithDomain( domain_, x, y, z );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::fillWithDomain( const flag_t domainSubFlag,
                                                                             const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT_EQUAL( domain_ & domainSubFlag, domainSubFlag );
   WALBERLA_ASSERT( field::isFlag( domainSubFlag ) );

   if( outerBB_.contains(x,y,z) && isEmpty(x,y,z) )
      addDomain( x, y, z, domainSubFlag );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::fillWithDomain( const CellInterval & cells )
{
   WALBERLA_ASSERT( field::isFlag( domain_ ) );
   fillWithDomain( domain_, cells );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::fillWithDomain( const flag_t domainSubFlag, const CellInterval & cells )
{
   WALBERLA_ASSERT_EQUAL( domain_ & domainSubFlag, domainSubFlag );
   WALBERLA_ASSERT( field::isFlag( domainSubFlag ) );

   CellInterval localCells( outerBB_ );
   localCells.intersect( cells );

   if( localCells.empty() )
      return;

   for( cell_idx_t z = localCells.zMin(); z <= localCells.zMax(); ++z ) {
      for( cell_idx_t y = localCells.yMin(); y <= localCells.yMax(); ++y ) {
         for( cell_idx_t x = localCells.xMin(); x <= localCells.xMax(); ++x )
         {
            if( isEmpty(x,y,z) )
               addDomain( x, y, z, domainSubFlag );
         }
      }
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::fillWithDomain( const CellIterator & begin, const CellIterator & end )
{
   fillWithDomain( domain_, begin, end );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::fillWithDomain( const flag_t domainSubFlag,
                                                                             const CellIterator & begin, const CellIterator & end )
{
   WALBERLA_ASSERT_EQUAL( domain_ & domainSubFlag, domainSubFlag );
   WALBERLA_ASSERT( field::isFlag( domainSubFlag ) );

   for( auto cell = begin; cell != end; ++cell )
   {
      const cell_idx_t x = cell->x();
      const cell_idx_t y = cell->y();
      const cell_idx_t z = cell->z();

      if( outerBB_.contains(x,y,z) && isEmpty(x,y,z) )
         addDomain( x, y, z, domainSubFlag );
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline shared_ptr<BoundaryConfiguration> BoundaryHandling< FlagField_T, Stencil, Boundaries... >::createBoundaryConfiguration( const BoundaryUID & uid,
                                                                                                                       const Config::BlockHandle & config ) const
{
   return createBoundaryConfiguration( boundaryConditions_, uid, config );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setBoundary( const FlagUID & flag,
                                                                          const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                          const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   setBoundary( flagField_->getFlag( flag ), x, y, z, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setBoundary( const flag_t flag,
                                                                          const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                          const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT_EQUAL( flag & boundary_, flag );
   WALBERLA_ASSERT( field::isFlag( flag ) );
   WALBERLA_ASSERT_EQUAL( numberOfMatchingBoundaryConditions( flag ), uint_t(1) );

   if( outerBB_.contains(x,y,z) )
      setBoundary( boundaryConditions_, flag, x, y, z, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setBoundary( const FlagUID & flag, const CellInterval & cells,
                                                                          const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   setBoundary( flagField_->getFlag( flag ), cells, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setBoundary( const flag_t flag, const CellInterval & cells,
                                                                          const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT_EQUAL( flag & boundary_, flag );
   WALBERLA_ASSERT( field::isFlag( flag ) );
   WALBERLA_ASSERT_EQUAL( numberOfMatchingBoundaryConditions( flag ), uint_t(1) );

   CellInterval localCells( outerBB_ );
   localCells.intersect( cells );

   if( localCells.empty() )
      return;

   setBoundary( boundaryConditions_, flag, localCells, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setBoundary( const FlagUID & flag, const CellIterator & begin, const CellIterator & end,
                                                                          const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   setBoundary( flagField_->getFlag( flag ), begin, end, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setBoundary( const flag_t flag, const CellIterator & begin, const CellIterator & end,
                                                                          const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT_EQUAL( flag & boundary_, flag );
   WALBERLA_ASSERT( field::isFlag( flag ) );
   WALBERLA_ASSERT_EQUAL( numberOfMatchingBoundaryConditions( flag ), uint_t(1) );

   CellVector localCells;
   for( auto cell = begin; cell != end; ++cell )
      if( outerBB_.contains(*cell) ) localCells.push_back( *cell );

   setBoundary( boundaryConditions_, flag, localCells, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceBoundary( const FlagUID & flag,
                                                                            const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                            const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   forceBoundary( flagField_->getFlag( flag ), x, y, z, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceBoundary( const flag_t flag,
                                                                            const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                            const BoundaryConfiguration & parameter )
{
   clear(x,y,z);
   setBoundary( flag, x, y, z, parameter );
}




template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceBoundary( const FlagUID & flag, const CellInterval & cells,
                                                                            const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   forceBoundary( flagField_->getFlag( flag ), cells, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceBoundary( const flag_t flag, const CellInterval & cells,
                                                                            const BoundaryConfiguration & parameter )
{
   clear( cells );
   setBoundary( flag, cells, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceBoundary( const FlagUID & flag, const CellIterator & begin, const CellIterator & end,
                                                                            const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   forceBoundary( flagField_->getFlag( flag ), begin, end, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceBoundary( const flag_t flag, const CellIterator & begin, const CellIterator & end,
                                                                            const BoundaryConfiguration & parameter )
{
   for( auto cell = begin; cell != end; ++cell )
   {
      clear( cell->x(), cell->y(), cell->z() );
      setBoundary( flag, cell->x(), cell->y(), cell->z(), parameter );
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeDomain( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   if( outerBB_.contains(x,y,z) )
   {
      flagField_->removeMask( x, y, z, domain_ );
      flagField_->removeFlag( x, y, z, nearBoundary_ );
      dirty_ = true;
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeDomain( const flag_t mask,
                                                                           const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   if( outerBB_.contains(x,y,z) )
   {
      flagField_->removeMask( x, y, z, static_cast< flag_t >( domain_ & mask ) );

      if( isEmpty(x,y,z) )
         flagField_->removeFlag( x, y, z, nearBoundary_ );

      dirty_ = true;
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeDomain( const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_GREATER_EQUAL( flagField_->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   CellInterval cells = getGhostLayerCellInterval( numberOfGhostLayersToInclude );
   removeDomain( cells );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeDomain( const CellInterval & cells )
{
   CellInterval localCells( outerBB_ );
   localCells.intersect( cells );

   if( localCells.empty() )
      return;

   for( auto cell = flagField_->beginSliceXYZ( localCells ); cell != flagField_->end(); ++cell )
   {
      removeMask( cell, domain_ );
      field::removeFlag( cell, nearBoundary_ );
   }

   dirty_ = true;
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeDomain( const flag_t mask, const CellInterval & cells )
{
   const flag_t dMask = mask & domain_;
   if( dMask == 0 )
      return;

   CellInterval localCells( outerBB_ );
   localCells.intersect( cells );

   if( localCells.empty() )
      return;

   for( auto cell = flagField_->beginSliceXYZ( localCells ); cell != flagField_->end(); ++cell )
   {
      removeMask( cell, dMask );

      if( isEmpty(cell) )
         field::removeFlag( cell, nearBoundary_ );
   }

   dirty_ = true;
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeDomain( const CellIterator & begin, const CellIterator & end )
{
   for( auto cell = begin; cell != end; ++cell )
      removeDomain( cell->x(), cell->y(), cell->z() );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeDomain( const flag_t mask, const CellIterator & begin, const CellIterator & end )
{
   const flag_t dMask = mask & domain_;
   if( dMask == 0 )
      return;

   for( auto cell = begin; cell != end; ++cell )
   {
      const cell_idx_t x = cell->x();
      const cell_idx_t y = cell->y();
      const cell_idx_t z = cell->z();

      if( outerBB_.contains(x,y,z) )
      {
         flagField_->removeMask( x, y, z, dMask );

         if( isEmpty(x,y,z) )
            flagField_->removeFlag( x, y, z, nearBoundary_ );
      }
   }

   dirty_ = true;
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeBoundary( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   if( outerBB_.contains(x,y,z) && flagField_->isPartOfMaskSet( x, y, z, boundary_ ) )
      removeBoundary( boundaryConditions_, x, y, z );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeBoundary( const FlagUID & flag,
                                                                             const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   if( !flagField_->flagExists( flag ) )
      return;

   removeBoundary( flagField_->getFlag( flag ), x, y, z );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeBoundary( const flag_t mask,
                                                                             const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   const flag_t bMask = mask & boundary_;
   if( bMask == 0 )
      return;

   if( outerBB_.contains(x,y,z) && flagField_->isPartOfMaskSet( x, y, z, bMask ) )
      removeBoundary( boundaryConditions_, x, y, z );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeBoundary( const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_GREATER_EQUAL( flagField_->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   CellInterval cells = getGhostLayerCellInterval( numberOfGhostLayersToInclude );
   removeBoundary( cells );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeBoundary( const CellInterval & cells )
{
   CellInterval localCells( outerBB_ );
   localCells.intersect( cells );

   if( localCells.empty() )
      return;

   CellInterval innerHull( localCells.xMin()+1, localCells.yMin()+1, localCells.zMin()+1,
                           localCells.xMax()-1, localCells.yMax()-1, localCells.zMax()-1 );

   if( !innerHull.empty() )
   {
      for( auto cell = flagField_->beginSliceXYZ( innerHull ); cell != flagField_->end(); ++cell )
      {
         if( isPartOfMaskSet( cell, boundary_ ) )
            removeBoundary( boundaryConditions_, cell.x(), cell.y(), cell.z(), false );
         field::removeFlag( cell, nearBoundary_ );
      }
      dirty_ = true;
   }

   for( auto cell = flagField_->beginSliceXYZ( localCells ); cell != flagField_->end(); ++cell )
   {
      const cell_idx_t x = cell.x();
      const cell_idx_t y = cell.y();
      const cell_idx_t z = cell.z();

      if( !innerHull.contains(x,y,z) )
      {
         if( isPartOfMaskSet( cell, boundary_ ) )
            removeBoundary( boundaryConditions_, x, y, z );

         if( isFlagSet( cell, nearBoundary_ ) )
         {
            bool remove = true;
            for( auto d = Stencil::beginNoCenter(); d != Stencil::end() && remove; ++d ) {
               if( isPartOfMaskSet( cell.neighbor( *d ), boundary_ ) )
                  remove = false;
            }
            if( remove )
            {
               WALBERLA_ASSERT( isPartOfMaskSet( cell, domain_ ) );
               field::removeFlag( cell, nearBoundary_ );
               dirty_ = true;
            }
         }
      }
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeBoundary( const FlagUID & flag, const CellInterval & cells )
{
   if( !flagField_->flagExists( flag ) )
      return;

   removeBoundary( flagField_->getFlag( flag ), cells );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeBoundary( const flag_t mask, const CellInterval & cells )
{
   const flag_t bMask = mask & boundary_;
   if( bMask == 0 )
      return;

   CellInterval localCells( outerBB_ );
   localCells.intersect( cells );

   for( auto cell = flagField_->beginSliceXYZ( localCells ); cell != flagField_->end(); ++cell )
      if( isPartOfMaskSet( cell, bMask ) )
         removeBoundary( boundaryConditions_, cell.x(), cell.y(), cell.z() );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeBoundary( const CellIterator & begin, const CellIterator & end )
{
   removeBoundary( boundary_, begin, end );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeBoundary( const FlagUID & flag, const CellIterator & begin, const CellIterator & end )
{
   if( !flagField_->flagExists( flag ) )
      return;

   removeBoundary( flagField_->getFlag( flag ), begin, end );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeBoundary( const flag_t mask, const CellIterator & begin, const CellIterator & end )
{
   const flag_t bMask = mask & boundary_;
   if( bMask == 0 )
      return;

   for( auto cell = begin; cell != end; ++cell )
   {
      const cell_idx_t x = cell->x();
      const cell_idx_t y = cell->y();
      const cell_idx_t z = cell->z();

      if( outerBB_.contains(x,y,z) && flagField_->isPartOfMaskSet( x, y, z, bMask) )
         removeBoundary( boundaryConditions_, x, y, z );
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setFlag( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                      const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   setFlag( flagField_->getFlag( flag ), x, y, z, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setFlag( const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                      const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( field::isFlag(flag) );

   if( ( flag & boundary_ ) == flag )
      setBoundary( flag, x, y, z, parameter );
   else if( ( flag & domain_ ) == flag )
      setDomain( flag, x, y, z );
   else if( outerBB_.contains(x,y,z) )
   {
      WALBERLA_ASSERT( !flagField_->isFlagSet( x, y, z, flag ) );
      flagField_->addFlag( x, y, z, flag );
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setFlag( const FlagUID & flag, const CellInterval & cells,
                                                                      const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   setFlag( flagField_->getFlag( flag ), cells, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setFlag( const flag_t flag, const CellInterval & cells,
                                                                      const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( field::isFlag(flag) );

   if( ( flag & boundary_ ) == flag )
      setBoundary( flag, cells, parameter );
   else if( ( flag & domain_ ) == flag )
      setDomain( flag, cells );
   else
   {
      CellInterval localCells( outerBB_ );
      localCells.intersect( cells );

      if( !localCells.empty() )
         for( auto cell = flagField_->beginSliceXYZ( localCells ); cell != flagField_->end(); ++cell )
         {
            WALBERLA_ASSERT( !field::isFlagSet( cell, flag ) );
            field::addFlag( cell, flag );
         }
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setFlag( const FlagUID & flag, const CellIterator & begin, const CellIterator & end,
                                                                      const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   setFlag( flagField_->getFlag( flag ), begin, end, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setFlag( const flag_t flag, const CellIterator & begin, const CellIterator & end,
                                                                      const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( field::isFlag(flag) );

   if( ( flag & boundary_ ) == flag )
      setBoundary( flag, begin, end, parameter );
   else if( ( flag & domain_ ) == flag )
      setDomain( flag, begin, end );
   else
   {
      for( auto cell = begin; cell != end; ++cell )
         if( outerBB_.contains( cell->x(), cell->y(), cell->z() ) )
         {
            WALBERLA_ASSERT( !flagField_->isFlagSet( cell->x(), cell->y(), cell->z(), flag ) );
            flagField_->addFlag( cell->x(), cell->y(), cell->z(), flag );
         }
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceFlag( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                        const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   forceFlag( flagField_->getFlag( flag ), x, y, z, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceFlag( const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                        const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( field::isFlag(flag) );

   if( ( flag & boundary_ ) == flag )
      forceBoundary( flag, x, y, z, parameter );
   else if( ( flag & domain_ ) == flag )
      forceDomain( flag, x, y, z );
   else if( outerBB_.contains(x,y,z) )
      flagField_->addFlag( x, y, z, flag );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceFlag( const FlagUID & flag, const CellInterval & cells,
                                                                        const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   forceFlag( flagField_->getFlag( flag ), cells, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceFlag( const flag_t flag, const CellInterval & cells,
                                                                        const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( field::isFlag(flag) );

   if( ( flag & boundary_ ) == flag )
      forceBoundary( flag, cells, parameter );
   else if( ( flag & domain_ ) == flag )
      forceDomain( flag, cells );
   else
   {
      CellInterval localCells( outerBB_ );
      localCells.intersect( cells );

      if( !localCells.empty() )
         for( auto cell = flagField_->beginSliceXYZ( localCells ); cell != flagField_->end(); ++cell )
            field::addFlag( cell, flag );
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceFlag( const FlagUID & flag, const CellIterator & begin, const CellIterator & end,
                                                                        const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   forceFlag( flagField_->getFlag( flag ), begin, end, parameter );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::forceFlag( const flag_t flag, const CellIterator & begin, const CellIterator & end,
                                                                        const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( field::isFlag(flag) );

   if( ( flag & boundary_ ) == flag )
      forceBoundary( flag, begin, end, parameter );
   else if( ( flag & domain_ ) == flag )
      forceDomain( flag, begin, end );
   else
   {
      for( auto cell = begin; cell != end; ++cell )
         if( outerBB_.contains( cell->x(), cell->y(), cell->z() ) )
            flagField_->addFlag( cell->x(), cell->y(), cell->z(), flag );
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeFlag( const FlagUID & flag, const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   removeFlag( flagField_->getFlag( flag ), numberOfGhostLayersToInclude );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeFlag( const flag_t flag, const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_GREATER_EQUAL( flagField_->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   CellInterval cells = getGhostLayerCellInterval( numberOfGhostLayersToInclude );
   removeFlag( flag, cells );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeFlag( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   removeFlag( flagField_->getFlag( flag ), x, y, z );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeFlag( const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   if( ( flag & boundary_ ) == flag )
      removeBoundary( flag, x, y, z );
   else if( ( flag & domain_ ) == flag )
      removeDomain( flag, x, y, z );
   else if( outerBB_.contains(x,y,z) )
      flagField_->removeFlag( x, y, z, flag );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeFlag( const FlagUID & flag, const CellInterval & cells )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   removeFlag( flagField_->getFlag( flag ), cells );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeFlag( const flag_t flag, const CellInterval & cells )
{
   if( ( flag & boundary_ ) == flag )
      removeBoundary( flag, cells );
   else if( ( flag & domain_ ) == flag )
      removeDomain( flag, cells );
   else
   {
      CellInterval localCells( outerBB_ );
      localCells.intersect( cells );

      if( !localCells.empty() )
         for( auto cell = flagField_->beginSliceXYZ( localCells ); cell != flagField_->end(); ++cell )
            field::removeFlag( cell, flag );
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeFlag( const FlagUID & flag, const CellIterator & begin, const CellIterator & end )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   removeFlag( flagField_->getFlag( flag ), begin, end );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeFlag( const flag_t flag, const CellIterator & begin, const CellIterator & end )
{
   if( ( flag & boundary_ ) == flag )
      removeBoundary( flag, begin, end );
   else if( ( flag & domain_ ) == flag )
      removeDomain( flag, begin, end );
   else
   {
      for( auto cell = begin; cell != end; ++cell )
         if( outerBB_.contains( cell->x(), cell->y(), cell->z() ) )
            flagField_->removeFlag( cell->x(), cell->y(), cell->z(), flag );
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::clear( const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_GREATER_EQUAL( flagField_->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   CellInterval cells = getGhostLayerCellInterval( numberOfGhostLayersToInclude );
   clear( cells );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::clear( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   if( outerBB_.contains(x,y,z) )
   {
      flagField_->removeMask( x, y, z, domain_ );
      flagField_->removeFlag( x, y, z, nearBoundary_ );

      dirty_ = true;

      if( flagField_->isPartOfMaskSet( x, y, z, boundary_ ) )
         removeBoundary( boundaryConditions_, x, y, z );
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::clear( const CellInterval & cells )
{
   CellInterval localCells( outerBB_ );
   localCells.intersect( cells );

   if( localCells.empty() )
      return;

   CellInterval innerHull( localCells.xMin()+1, localCells.yMin()+1, localCells.zMin()+1,
                           localCells.xMax()-1, localCells.yMax()-1, localCells.zMax()-1 );

   if( !innerHull.empty() )
   {
      for( auto cell = flagField_->beginSliceXYZ( innerHull ); cell != flagField_->end(); ++cell )
      {
         if( isPartOfMaskSet( cell, boundary_ ) )
            removeBoundary( boundaryConditions_, cell.x(), cell.y(), cell.z(), false );
      }
   }

   for( auto cell = flagField_->beginSliceXYZ( localCells ); cell != flagField_->end(); ++cell )
   {
      const cell_idx_t x = cell.x();
      const cell_idx_t y = cell.y();
      const cell_idx_t z = cell.z();

      removeMask( cell, domain_ );
      field::removeFlag( cell, nearBoundary_ );

      if( !innerHull.contains(x,y,z) && isPartOfMaskSet( cell, boundary_ ) )
         removeBoundary( boundaryConditions_, x, y, z );
   }

   dirty_ = true;
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::clear( const CellIterator & begin, const CellIterator & end )
{
   for( auto cell = begin; cell != end; ++cell )
      clear( cell->x(), cell->y(), cell->z() );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::operator()( const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT_LESS( numberOfGhostLayersToInclude, flagField_->nrOfGhostLayers() );
   WALBERLA_ASSERT( checkConsistency( numberOfGhostLayersToInclude ) );

   CellInterval cells = getGhostLayerCellInterval( numberOfGhostLayersToInclude );

   beforeBoundaryTreatment();

   if( mode_ == OPTIMIZED_SPARSE_TRAVERSAL )
   {
      WALBERLA_ASSERT( innerBB_.contains( cells ) );

      if( dirty_ )
      {
         for( uint_t i = 0; i != rebuildCellDirectionPairs_.size(); ++i )
            rebuildCellDirectionPairs_[i] = true;
         dirty_ = false;
      }

      auto & cellDirectionPairs = cellDirectionPairs_[numberOfGhostLayersToInclude];

      if( rebuildCellDirectionPairs_[numberOfGhostLayersToInclude] )
      {
         cellDirectionPairs.clear();
         cellDirectionPairs.resize( bcMaskMapping_.size() );

         for( auto cell = flagField_->beginSliceXYZ( cells ); cell != flagField_->end(); ++cell )
            if( isFlagSet( cell, nearBoundary_ ) )
               for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
                  if( field::isPartOfMaskSet( cell.neighbor( *d ), boundary_ ) )
                  {
                     uint_t index( uint_t(0) );
                     for( auto mask = bcMaskMapping_.begin(); mask != bcMaskMapping_.end(); ++mask, ++index )
                        if( field::isPartOfMaskSet( cell.neighbor( *d ), *mask ) )
                        {
                           cellDirectionPairs[ index ].push_back( std::make_pair( cell.cell(), *d) );
                           break;
                        }
                  }

         rebuildCellDirectionPairs_[numberOfGhostLayersToInclude] = false;
      }

      WALBERLA_ASSERT( checkFlagField( numberOfGhostLayersToInclude ) );

      if( !cellDirectionPairs.empty() )
         treatDirection( boundaryConditions_, uint_t(0), cellDirectionPairs );
   }
   else
   {
      (*this)( cells );
   }

   afterBoundaryTreatment();
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT( innerBB_.contains(x,y,z) );

   if( isNearBoundary(x,y,z) )
   {
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         const cell_idx_t nx = x + cell_idx_c( d.cx() );
         const cell_idx_t ny = y + cell_idx_c( d.cy() );
         const cell_idx_t nz = z + cell_idx_c( d.cz() );

         if( isBoundary(nx,ny,nz) )
            treatDirection( boundaryConditions_, x, y, z, *d, nx, ny, nz );
      }
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::operator()( const CellInterval & cells )
{
   CellInterval localCells( innerBB_ );
   localCells.intersect( cells );

   if( localCells.empty() )
      return;

   WALBERLA_ASSERT( checkConsistency( localCells ) );

   #if defined(_OPENMP) && !(defined(_MSC_VER) && _MSC_VER < 1925)
   const int zMin = int_c( localCells.zMin() );
   const int zMax = int_c( localCells.zMax() );
   #pragma omp parallel for schedule(static) if(threadSafeBCs_)
   for( int iz = zMin; iz <= zMax; ++iz ) {
      cell_idx_t z = cell_idx_c( iz );
   #else
   for( cell_idx_t z = localCells.zMin(); z <= localCells.zMax(); ++z ) {
   #endif
      for( cell_idx_t y = localCells.yMin(); y <= localCells.yMax(); ++y ) {
         for( cell_idx_t x = localCells.xMin(); x <= localCells.xMax(); ++x )
         {
            (*this)(x,y,z);
         }
      }
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename CellIterator >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::operator()( const CellIterator & begin, const CellIterator & end )
{
   for( auto cell = begin; cell != end; ++cell )
   {
      const cell_idx_t x = cell.x();
      const cell_idx_t y = cell.y();
      const cell_idx_t z = cell.z();

      if( innerBB_.contains(x,y,z) )
         (*this)(x,y,z);
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename Buffer_T >
void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::pack( Buffer_T & buffer, const CellInterval & interval, const bool assumeIdenticalFlagMapping ) const
{
   WALBERLA_UNUSED( assumeIdenticalFlagMapping );

#ifdef NDEBUG
   if( !assumeIdenticalFlagMapping )
#endif
      buffer << getFlagMapping();

#ifndef NDEBUG
   const uint_t numberOfCells = uint_c( interval.max()[0] + 1 - interval.min()[0] ) * uint_c( interval.max()[1] + 1 - interval.min()[1] ) *
                                uint_c( interval.max()[2] + 1 - interval.min()[2] );
   buffer << numberOfCells;
#endif

   for( auto z = interval.min()[2]; z <= interval.max()[2]; ++z ) {
      for( auto y = interval.min()[1]; y <= interval.max()[1]; ++y ) {
         for( auto x = interval.min()[0]; x <= interval.max()[0]; ++x )
         {
            flag_t mask = flagField_->get(x,y,z);
            field::removeFlag( mask, nearBoundary_ );
            buffer << mask;
            pack( buffer, mask, x, y, z );
         }
      }
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename Buffer_T >
void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::unpack( Buffer_T & buffer, const CellInterval & interval, const bool assumeIdenticalFlagMapping )
{
   bool identicalFlagMapping = false;
   std::vector< flag_t > flagMapping = getNeighborFlagMapping( buffer, assumeIdenticalFlagMapping, identicalFlagMapping ); // neighbor-flag_t -> flag_t

   clear( interval );

   flag_t handledFlags = boundary_ | domain_;

#ifndef NDEBUG
   uint_t numberOfCells;
   buffer >> numberOfCells;
   WALBERLA_ASSERT_EQUAL( numberOfCells, uint_c( interval.max()[0] + 1 - interval.min()[0] ) * uint_c( interval.max()[1] + 1 - interval.min()[1] ) *
                                uint_c( interval.max()[2] + 1 - interval.min()[2] ) );
#endif

   for( auto z = interval.min()[2]; z <= interval.max()[2]; ++z ) {
      for( auto y = interval.min()[1]; y <= interval.max()[1]; ++y ) {
         for( auto x = interval.min()[0]; x <= interval.max()[0]; ++x )
         {
            flag_t mask;
            buffer >> mask;

            if( !identicalFlagMapping )
               translateMask( mask, flagMapping );

            (*flagField_)(x,y,z) = static_cast< flag_t >( mask & ~handledFlags );
            flag_t flag          = static_cast< flag_t >( mask &  handledFlags );
            if( flag )
               unpack( buffer, flag, x, y, z );
         }
      }
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename Buffer_T >
void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::pack( Buffer_T & buffer, stencil::Direction direction, const uint_t numberOfLayers,
                                                            const bool assumeIdenticalFlagMapping ) const
{
   CellInterval interval = getPackingInterval( direction, numberOfLayers );
   pack( buffer, interval, assumeIdenticalFlagMapping );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename Buffer_T >
void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::unpack( Buffer_T & buffer, stencil::Direction direction, const uint_t numberOfLayers,
                                                              const bool assumeIdenticalFlagMapping )
{
   CellInterval interval = getUnpackingInterval( direction, numberOfLayers );
   unpack( buffer, interval, assumeIdenticalFlagMapping );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::toStream( std::ostream & os ) const {

   os << "==================== BoundaryHandling ====================\n\n"
      << "Identifier: " << uid_.getIdentifier() << "\n\n"
      << "Boundary Conditions:\n";

   std::vector< std::string > boundaryConditions;
   getBoundaryConditions( boundaryConditions_, boundaryConditions );

   for( auto bc = boundaryConditions.begin(); bc != boundaryConditions.end(); ++bc )
      os << "- " << *bc << "\n";

   os << "\nFlags/Masks:"
      << "\n- near boundary: "; valueToStream( os, nearBoundary_ );
   os << "\n- boundary: "; valueToStream( os, boundary_ );
   os << "\n- domain: "; valueToStream( os, domain_ );

   os << "\n\nAssociated Flag Field:\n";
   flagField_->printRegistered( os );

   os << "\n==================== BoundaryHandling ====================\n";
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline std::string BoundaryHandling< FlagField_T, Stencil, Boundaries... >::toString() const {

   std::ostringstream oss;
   toStream( oss );
   return oss.str();
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
CellInterval BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getGhostLayerCellInterval( const uint_t numberOfGhostLayersToInclude ) const
{
   CellInterval cells( -cell_idx_c( numberOfGhostLayersToInclude ),
                       -cell_idx_c( numberOfGhostLayersToInclude ),
                       -cell_idx_c( numberOfGhostLayersToInclude ),
                        cell_idx_c( flagField_->xSize() + numberOfGhostLayersToInclude ) - 1,
                        cell_idx_c( flagField_->ySize() + numberOfGhostLayersToInclude ) - 1,
                        cell_idx_c( flagField_->zSize() + numberOfGhostLayersToInclude ) - 1 );
   return cells;
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
typename std::enable_if<(N!=-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setupBoundaryConditions( BoundariesTuple & boundaryConditions )
{
   using BoundaryType = typename std::tuple_element<N, BoundariesTuple>::type;
   BoundaryType & boundaryCondition = std::get<N>( boundaryConditions );

   if( numberOfMatchingBoundaryConditions( boundaryCondition.getUID() ) != 1 )
      WALBERLA_ABORT( "Every boundary condition registered at the same boundary handler must have a unique boundary UID!\n"
                      "The boundary UID \"" << boundaryCondition.getUID() << "\" is not unique for boundary handler \"" << uid_.getIdentifier() << "\"." );

   flag_t mask(0);
   std::vector< FlagUID > uids;
   boundaryCondition.pushFlags( uids );

   for( auto uid = uids.begin(); uid != uids.end(); ++uid )
   {
      if( flagField_->flagExists( *uid ) )
         mask = static_cast<flag_t>( mask | flagField_->getFlag( *uid ) );
      else
         mask = static_cast<flag_t>( mask | flagField_->registerFlag( *uid ) );
   }
   WALBERLA_ASSERT_EQUAL( boundary_ & mask, flag_t(0) ); // every boundary condition must have a unique mask/set of FlagUIDs

   boundaryCondition.setMask( mask );

   boundary_ = static_cast< flag_t >( boundary_ | mask );

   bcMaskMapping_.push_back( mask );

   threadSafeBCs_ &= isThreadSafe< BoundaryType >::value;

   setupBoundaryConditions< BoundariesTuple, N-1 >( boundaryConditions );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline std::vector< BoundaryUID > BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getBoundaryUIDs() const
{
   std::vector< BoundaryUID > uids;
   getBoundaryUIDs( boundaryConditions_, uids );
   return uids;
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getBoundaryUIDs( const BoundariesTuple & boundaryConditions,
                                                                              std::vector< BoundaryUID > & uids ) const
{
   uids.push_back( std::get<N>( boundaryConditions ).getUID() );
   getBoundaryUIDs< BoundariesTuple, N-1 >( boundaryConditions, uids );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), BoundaryUID>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getBoundaryUID( const BoundariesTuple & boundaryConditions,
                                                                                    const flag_t flag ) const
{
   const auto & boundaryCondition = std::get<N>( boundaryConditions );

   if( ( boundaryCondition.getMask() & flag ) == flag )
   {
      return boundaryCondition.getUID();
   }
   else
   {
      return getBoundaryUID< BoundariesTuple, N-1 >( boundaryConditions, flag );
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N==-1), BoundaryUID>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getBoundaryUID( const BoundariesTuple &,
                                                                                    const flag_t flag ) const
{
   if( !flagField_->isRegistered( flag ) )
      WALBERLA_ABORT( "The requested flag with value " << flag << " is not registered at the flag field and is not handled "\
                      "by any boundary condition of boundary handling " << uid_.getIdentifier() << "!" );

   const FlagUID & flagUID = flagField_->getFlagUID( flag );
   WALBERLA_ABORT( "The requested flag " << flagUID.getIdentifier() << " is not handled by any boundary condition of "\
                   "boundary handling " << uid_.getIdentifier() << "!" );

#ifdef __IBMCPP__
   return *(reinterpret_cast< BoundaryUID * >( NULL )); // silencing incorrect IBM compiler warning
#endif
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), bool>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::containsBoundaryCondition( const BoundariesTuple & boundaryConditions,
                                                                                        const BoundaryUID & uid ) const
{
   if( std::get<N>( boundaryConditions ).getUID() == uid )
      return true;
   return containsBoundaryCondition< BoundariesTuple, N-1 >( boundaryConditions, uid );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), typename BoundaryHandling< FlagField_T, Stencil, Boundaries... >::flag_t>::type
   BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getBoundaryMask( const BoundariesTuple & boundaryConditions,
                                                                     const BoundaryUID & uid ) const
{
   const auto & boundaryCondition = std::get<N>( boundaryConditions );

   if( boundaryCondition.getUID() == uid )
      return boundaryCondition.getMask();
   return getBoundaryMask< BoundariesTuple, N-1 >( boundaryConditions, uid );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline uint_t BoundaryHandling< FlagField_T, Stencil, Boundaries... >::numberOfMatchingBoundaryConditions( const BoundaryUID & uid ) const
{
   return numberOfMatchingBoundaryConditions( boundaryConditions_, uid );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), uint_t>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::numberOfMatchingBoundaryConditions( const BoundariesTuple & boundaryConditions,
                                                                                                 const BoundaryUID & uid ) const
{
   return ( ( std::get<N>( boundaryConditions ).getUID() == uid ) ? uint_c(1) : uint_c(0) ) +
          numberOfMatchingBoundaryConditions< BoundariesTuple, N-1 >( boundaryConditions, uid );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), uint_t>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::numberOfMatchingBoundaryConditions( const BoundariesTuple & boundaryConditions,
                                                                                                 const flag_t mask ) const
{
   return ( ( ( std::get<N>( boundaryConditions ).getMask() & mask ) != 0 ) ? uint_c(1) : uint_c(0) ) +
          numberOfMatchingBoundaryConditions< BoundariesTuple, N-1 >( boundaryConditions, mask );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline bool BoundaryHandling< FlagField_T, Stencil, Boundaries... >::checkFlagField( const uint_t numberOfGhostLayersToInclude ) const
{
   if( !flagField_->isRegistered( nearBoundary_ ) )
      return false;

   if( !flagField_->isRegistered( boundary_ ) )
      return false;

   if( !flagField_->isRegistered( domain_ ) )
      return false;

   if( mode_ == OPTIMIZED_SPARSE_TRAVERSAL )
   {
      WALBERLA_ASSERT_LESS( numberOfGhostLayersToInclude, flagField_->nrOfGhostLayers() );

      CellInterval cells = getGhostLayerCellInterval( numberOfGhostLayersToInclude );

      WALBERLA_ASSERT( innerBB_.contains( cells ) );

      CellVector nearBoundaryCells;
      for( auto cellDirectionPairs = cellDirectionPairs_[numberOfGhostLayersToInclude].begin();
               cellDirectionPairs != cellDirectionPairs_[numberOfGhostLayersToInclude].end(); ++cellDirectionPairs )
         for( auto cellDirectionPair = cellDirectionPairs->begin(); cellDirectionPair != cellDirectionPairs->end(); ++cellDirectionPair )
            nearBoundaryCells.push_back( cellDirectionPair->first );

      for( auto cell = nearBoundaryCells.begin(); cell != nearBoundaryCells.end(); ++cell )
         if( !flagField_->isFlagSet( cell->x(), cell->y(), cell->z(), nearBoundary_ ) )
            return false;

      CellSet nearBoundarySet( nearBoundaryCells );

      for( auto cell = flagField_->beginSliceXYZ( cells ); cell != flagField_->end(); ++cell )
         if( isFlagSet( cell, nearBoundary_ ) )
            if( !nearBoundarySet.contains( cell.x(), cell.y(), cell.z() ) )
               return false;
   }

   return true;
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::addDomain( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const flag_t domain )
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );
   WALBERLA_ASSERT_EQUAL( domain_ & domain, domain );
   WALBERLA_ASSERT( field::isFlag( domain ) );
   WALBERLA_ASSERT( isEmpty(x,y,z) );

   flagField_->addFlag( x, y, z, domain );

   // near boundary?

   if( innerBB_.contains(x,y,z) )
   {
      for( auto d = Stencil::beginNoCenter(); d != Stencil::end(); ++d )
      {
         if( flagField_->isPartOfMaskSet( x + d.cx(), y + d.cy(), z + d.cz(), boundary_ ) )
         {
            flagField_->addFlag( x, y, z, nearBoundary_ );
            dirty_ = true;
            break;
         }
      }
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), shared_ptr<BoundaryConfiguration>>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::createBoundaryConfiguration( const BoundariesTuple & boundaryConditions,
                                                                                                                       const BoundaryUID & uid, const Config::BlockHandle & config ) const
{
   using BoundaryType = typename std::tuple_element<N, BoundariesTuple>::type;
   const BoundaryType & boundaryCondition = std::get<N>( boundaryConditions );
   
   if( boundaryCondition.getUID() == uid )
      return BoundaryType::createConfiguration( config );

   return createBoundaryConfiguration< BoundariesTuple, N-1 >( boundaryConditions, uid, config );
}

template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N==-1), shared_ptr<BoundaryConfiguration>>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::createBoundaryConfiguration( const BoundariesTuple &,
                                                                                                                       const BoundaryUID & uid, const Config::BlockHandle & ) const
{
   WALBERLA_CHECK( false, "There is no boundary condition registered at boundary handling " << uid_ << " for a boundary with UID" << uid << "." );
   return make_shared<BoundaryConfiguration>();
}




template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::addNearBoundary( const CellInterval & cells )
{
   WALBERLA_ASSERT( innerBB_.contains( cells ) );

   for( auto cell = flagField_->beginSliceXYZ( cells ); cell != flagField_->end(); ++cell ) {
      if( isPartOfMaskSet( cell, domain_ ) ) {
         for( auto d = Stencil::beginNoCenter(); d != Stencil::end(); ++d ) { // Even if a domain cell is geometrically adjacent to a boundary cell,
            if( isPartOfMaskSet( cell.neighbor( *d ), boundary_ ) )           // it is not automatically "near" this boundary cell with respect to the stencil.
            {
               addFlag( cell, nearBoundary_ );
               dirty_ = true;
               break;
            }
         }
      }
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::addBoundary( const flag_t flag,
                                                                          const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );
   WALBERLA_ASSERT( !flagField_->isPartOfMaskSet( x, y, z, boundary_ ) ); // no existing boundary
   WALBERLA_ASSERT( !flagField_->isPartOfMaskSet( x, y, z, domain_ ) );   // no domain
   WALBERLA_ASSERT( !flagField_->isFlagSet( x, y, z, nearBoundary_ ) );   // no near boundary marker

   // setting boundary flag

   flagField_->addFlag( x, y, z, flag );

   // setting near boundary flag

   CellInterval hull( x-1, y-1, z-1, x+1, y+1, z+1 );
   hull.intersect( innerBB_ );

   addNearBoundary( hull );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setBoundary( BoundariesTuple & boundaryConditions, const flag_t flag,
                                                                          const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                          const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );

   auto & boundaryCondition = std::get<N>( boundaryConditions );
   if( ( boundaryCondition.getMask() & flag ) == flag )
   {
      addBoundary( flag, x, y, z );

      boundaryCondition.registerCell( flag, x, y, z, parameter );
   }
   else
      setBoundary< BoundariesTuple, N-1 >( boundaryConditions, flag, x, y, z, parameter );
}

template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N==-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setBoundary( const BoundariesTuple &, const flag_t flag,
                                                                          const cell_idx_t, const cell_idx_t, const cell_idx_t,
                                                                          const BoundaryConfiguration & ) const
{
   if( flagField_->isRegistered( flag ) )
   {
      WALBERLA_CHECK( false, "You are trying to set a boundary at boundary handling " << uid_ << " with flag " << flag <<
                             " (" << flagField_->getFlagUID( flag ) << ").\nHowever, no boundary condition is registered for this flag!" );
   }

   WALBERLA_CHECK( false, "You are trying to set a boundary at boundary handling " << uid_ << " with flag " << flag <<
                          ".\nHowever, no boundary condition is registered for this flag!" );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setBoundary( BoundariesTuple & boundaryConditions, const flag_t flag,
                                                                          const CellInterval & cells, const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( outerBB_.contains( cells ) );
   WALBERLA_ASSERT( !cells.empty() );

   auto & boundaryCondition = std::get<N>( boundaryConditions );
   if( ( boundaryCondition.getMask() & flag ) == flag )
   {
      // setting boundary flag

      for( auto cell = flagField_->beginSliceXYZ( cells ); cell != flagField_->end(); ++cell )
      {
         WALBERLA_ASSERT( !isPartOfMaskSet( cell, boundary_ ) ); // no existing boundary
         WALBERLA_ASSERT( !isPartOfMaskSet( cell, domain_ ) );   // no domain
         WALBERLA_ASSERT( !isFlagSet( cell, nearBoundary_ ) );   // no near boundary marker

         addFlag( cell, flag );
      }

      // setting near boundary flag

      CellInterval hull( cells.xMin()-1, cells.yMin()-1, cells.zMin()-1, cells.xMax()+1, cells.yMax()+1, cells.zMax()+1 );
      hull.intersect( innerBB_ );

      addNearBoundary( hull );

      boundaryCondition.registerCells( flag, cells, parameter );
   }
   else
      setBoundary< BoundariesTuple, N-1 >( boundaryConditions, flag, cells, parameter );
}

template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N==-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setBoundary( const BoundariesTuple &, const flag_t flag,
                                                                          const CellInterval &, const BoundaryConfiguration & ) const
{
   if( flagField_->isRegistered( flag ) )
   {
      WALBERLA_CHECK( false, "You are trying to set a boundary at boundary handling " << uid_ << " with flag " << flag <<
                             " (" << flagField_->getFlagUID( flag ) << ").\nHowever, no boundary condition is registered for this flag!" );
   }

   WALBERLA_CHECK( false, "You are trying to set a boundary at boundary handling " << uid_ << " with flag " << flag <<
                          ".\nHowever, no boundary condition is registered for this flag!" );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setBoundary( BoundariesTuple & boundaryConditions, const flag_t flag,
                                                                          const CellVector & cells, const BoundaryConfiguration & parameter )
{
   if( cells.empty() )
      return;

   auto & boundaryCondition = std::get<N>( boundaryConditions );
   if( ( boundaryCondition.getMask() & flag ) == flag )
   {
      for( auto cell = cells.begin(); cell != cells.end(); ++cell )
         addBoundary( flag, cell->x(), cell->y(), cell->z() );

      boundaryCondition.registerCells( flag, cells.begin(), cells.end(), parameter );
   }
   else
      setBoundary< BoundariesTuple, N-1 >( boundaryConditions, flag, cells, parameter );
}

template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N==-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::setBoundary( const BoundariesTuple &, const flag_t flag,
                                                                          const CellVector &, const BoundaryConfiguration & ) const
{
   if( flagField_->isRegistered( flag ) )
   {
      WALBERLA_CHECK( false, "You are trying to set a boundary at boundary handling " << uid_ << " with flag " << flag <<
                             " (" << flagField_->getFlagUID( flag ) << ").\nHowever, no boundary condition is registered for this flag!" );
   }

   WALBERLA_CHECK( false, "You are trying to set a boundary at boundary handling " << uid_ << " with flag " << flag <<
                          ".\nHowever, no boundary condition is registered for this flag!" );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::removeBoundary( BoundariesTuple & boundaryConditions,
                                                                             const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                             const bool checkNearBoundaryFlags )
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );
   WALBERLA_ASSERT( flagField_->isPartOfMaskSet( x, y, z, boundary_ ) );

   auto & boundaryCondition = std::get<N>( boundaryConditions );

   if( flagField_->isPartOfMaskSet( x, y, z, boundaryCondition.getMask() ) )
   {
      boundaryCondition.unregisterCell( flagField_->get(x,y,z) & boundaryCondition.getMask(), x, y, z );

      // remove boundary

      flagField_->removeMask( x, y, z, boundary_ );

      WALBERLA_ASSERT( isEmpty(x,y,z) );
      WALBERLA_ASSERT( !isNearBoundary(x,y,z) );

      if( checkNearBoundaryFlags )
      {
         // potentially clear near boundary flags

         CellInterval hull( x-1, y-1, z-1, x+1, y+1, z+1 );
         hull.intersect( innerBB_ );

         bool neighborIsNearBoundary = false;
         for( auto cell = flagField_->beginSliceXYZ( hull ); cell != flagField_->end(); ++cell ) {
            if( isFlagSet( cell, nearBoundary_ ) ) {
               neighborIsNearBoundary = true;

               bool remove = true;
               for( auto d = Stencil::beginNoCenter(); d != Stencil::end() && remove; ++d ) {
                  if( isPartOfMaskSet( cell.neighbor( *d ), boundary_ ) )
                     remove = false;
               }
               if( remove )
               {
                  WALBERLA_ASSERT( isPartOfMaskSet( cell, domain_ ) );
                  field::removeFlag( cell, nearBoundary_ );
               }
            }
         }
         if( neighborIsNearBoundary )
            dirty_ = true;
      }
   }
   else
   {
      removeBoundary< BoundariesTuple, N-1 >( boundaryConditions, x, y, z, checkNearBoundaryFlags );
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::treatDirection( BoundariesTuple & boundaryConditions, const uint_t index,
                                                                             const std::vector< std::vector< std::pair< Cell, stencil::Direction > > > & cellDirectionPairs )
{
   auto & boundaryCondition = std::get<N>( boundaryConditions );

   WALBERLA_ASSERT_LESS( index, cellDirectionPairs.size() );

   const int size = int_c( cellDirectionPairs[index].size() );
   #ifdef _OPENMP
   #pragma omp parallel for schedule(static) if(threadSafeBCs_)
   #endif
   for( int i = 0; i < size; ++i )
   {
      const auto & cell      = cellDirectionPairs[index][uint_c(i)].first;
      const auto & direction = cellDirectionPairs[index][uint_c(i)].second;

      const cell_idx_t x = cell.x();
      const cell_idx_t y = cell.y();
      const cell_idx_t z = cell.z();
      const cell_idx_t nx = x + cell_idx_c( stencil::cx[ direction ] );
      const cell_idx_t ny = y + cell_idx_c( stencil::cy[ direction ] );
      const cell_idx_t nz = z + cell_idx_c( stencil::cz[ direction ] );

      boundaryCondition.treatDirection( x, y, z, direction, nx, ny, nz, flagField_->get(nx,ny,nz) );
   }

   treatDirection< BoundariesTuple, N-1 >( boundaryConditions, index + uint_t(1), cellDirectionPairs );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::treatDirection( BoundariesTuple & boundaryConditions,
                                                                             const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                             const stencil::Direction dir,
                                                                             const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz )
{
   auto & boundaryCondition = std::get<N>( boundaryConditions );

   if( flagField_->isPartOfMaskSet( nx, ny, nz, boundaryCondition.getMask() ) )
   {
      boundaryCondition.treatDirection( x, y, z, dir, nx, ny, nz, flagField_->get(nx,ny,nz) );
   }
   else
   {
      treatDirection< BoundariesTuple, N-1 >( boundaryConditions, x, y, z, dir, nx, ny, nz );
   }
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::beforeBoundaryTreatment( BoundariesTuple & boundaryConditions )
{
   std::get<N>( boundaryConditions ).beforeBoundaryTreatment();

   beforeBoundaryTreatment< BoundariesTuple, N-1 >( boundaryConditions );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::afterBoundaryTreatment( BoundariesTuple & boundaryConditions )
{
   std::get<N>( boundaryConditions ).afterBoundaryTreatment();

   afterBoundaryTreatment< BoundariesTuple, N-1 >( boundaryConditions );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline std::map< std::string, typename BoundaryHandling< FlagField_T, Stencil, Boundaries... >::flag_t >
BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getFlagMapping() const
{
   std::map< std::string, flag_t > mapping;
   const auto uidMapping = flagField_->getMapping();
   for( auto it = uidMapping.begin(); it != uidMapping.end(); ++it )
      mapping[ it->first.getIdentifier() ] = it->second;
   return mapping;
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename Buffer_T >
std::vector< typename BoundaryHandling< FlagField_T, Stencil, Boundaries... >::flag_t >
BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getNeighborFlagMapping( Buffer_T & buffer, const bool assumeIdenticalFlagMapping,
                                                                         bool & identicalFlagMapping ) const
{
   identicalFlagMapping = assumeIdenticalFlagMapping;
   std::vector< flag_t > flagMapping; // neighbor-flag_t -> flag_t

#ifdef NDEBUG
   if( !assumeIdenticalFlagMapping )
   {
#endif
      std::map< std::string, flag_t >       myMapping = getFlagMapping();
      std::map< std::string, flag_t > neighborMapping;
      buffer >> neighborMapping;

#ifndef NDEBUG
      if( assumeIdenticalFlagMapping )
      {
         WALBERLA_ASSERT_EQUAL( myMapping.size(), neighborMapping.size() );
         WALBERLA_ASSERT( std::equal( myMapping.begin(), myMapping.end(), neighborMapping.begin() ) );
      }
      else {
#endif

      if( myMapping.size() == neighborMapping.size() && std::equal( myMapping.begin(), myMapping.end(), neighborMapping.begin() ) )
         identicalFlagMapping = true;
      else
      {
         for( auto neighbor = neighborMapping.begin(); neighbor != neighborMapping.end(); ++neighbor )
         {
            WALBERLA_ASSERT( field::isFlag(neighbor->second) );
            if( !flagField_->flagExists( neighbor->first ) )
               WALBERLA_ABORT( "There exists no flag with identifier \"" << neighbor->first << "\"!" );

            flagMapping.push_back( neighbor->second );
            flagMapping.push_back( flagField_->getFlag( neighbor->first ) );
         }
      }
   }

   return flagMapping;
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::translateMask( flag_t & mask, const std::vector< flag_t > & flagMapping ) const
{
   flag_t neighbor( mask );
   mask = static_cast< flag_t >(0);
   for( uint_t i = 0; i != flagMapping.size(); i +=2 )
      if( field::isFlagSet( neighbor, flagMapping[i] ) )
         field::addFlag( mask, flagMapping[i+1] );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
CellInterval BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getPackingInterval( stencil::Direction direction, const uint_t numberOfLayers ) const
{
   CellInterval interval = getUnpackingInterval( direction, numberOfLayers );

   for( uint_t i = 0; i != 3; ++i )
   {
      const auto offset = cell_idx_c( stencil::c[i][direction] ) * cell_idx_c( numberOfLayers );
      interval.min()[i] -= offset;
      interval.max()[i] -= offset;
   }

   return interval;
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
CellInterval BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getUnpackingInterval( stencil::Direction direction, const uint_t numberOfLayers ) const
{
   WALBERLA_ASSERT_GREATER_EQUAL( numberOfLayers, uint_t(1) );
   WALBERLA_ASSERT( stencil::cx[direction] == 0 || outerBB_.xSize() >= ( uint_c(4) * numberOfLayers ) );
   WALBERLA_ASSERT( stencil::cy[direction] == 0 || outerBB_.ySize() >= ( uint_c(4) * numberOfLayers ) );
   WALBERLA_ASSERT( stencil::cz[direction] == 0 || outerBB_.zSize() >= ( uint_c(4) * numberOfLayers ) );

   CellInterval interval( outerBB_ );

   for( uint_t i = 0; i != 3; ++i )
   {
      const auto c = cell_idx_c( stencil::c[i][direction] );

      if( c == -1 ) interval.max()[i] = interval.min()[i] + cell_idx_c( numberOfLayers - 1 );
      else if( c == 1 ) interval.min()[i] = interval.max()[i] - cell_idx_c( numberOfLayers - 1 );
      else
      {
         WALBERLA_ASSERT_EQUAL( c, 0 );
         interval.min()[i] += cell_idx_c( numberOfLayers );
         interval.max()[i] -= cell_idx_c( numberOfLayers );
      }
   }

   return interval;
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename Buffer_T >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::pack( Buffer_T & buffer, const flag_t mask,
                                                                   const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   if( field::isPartOfMaskSet( mask, boundary_ ) )
      pack( boundaryConditions_, buffer, mask, x, y, z );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename Buffer_T, typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::pack( const BoundariesTuple & boundaryConditions, Buffer_T & buffer,
                                                                   const flag_t mask, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   const auto & boundaryCondition = std::get<N>( boundaryConditions );
   if( field::isPartOfMaskSet( mask, boundaryCondition.getMask() ) )
   {
      boundaryCondition.packCell( buffer, x, y, z );
   }
   else
      pack< Buffer_T, BoundariesTuple, N-1 >( boundaryConditions, buffer, mask, x, y, z );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename Buffer_T >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::unpack( Buffer_T & buffer, const flag_t flag,
                                                                     const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT( field::isFlag(flag) );

   if( ( flag & boundary_ ) == flag )
      unpackBoundary( buffer, flag, x, y, z );
   else if( ( flag & domain_ ) == flag )
      setDomain( flag, x, y, z );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename Buffer_T >
inline void BoundaryHandling< FlagField_T, Stencil, Boundaries... >::unpackBoundary( Buffer_T & buffer, const flag_t flag,
                                                                             const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT_EQUAL( flag & boundary_, flag );
   WALBERLA_ASSERT( field::isFlag( flag ) );
   WALBERLA_ASSERT_EQUAL( numberOfMatchingBoundaryConditions( flag ), uint_t(1) );
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );

   unpackBoundary( boundaryConditions_, buffer, flag, x, y, z );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename Buffer_T, typename BoundariesTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::unpackBoundary( BoundariesTuple & boundaryConditions, Buffer_T & buffer,
                                                                             const flag_t flag,
                                                                             const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );

   auto & boundaryCondition = std::get<N>( boundaryConditions );
   if( ( boundaryCondition.getMask() & flag ) == flag )
   {
      addBoundary( flag, x, y, z );

      boundaryCondition.registerCell( buffer, flag, x, y, z );
   }
   else
      unpackBoundary< Buffer_T, BoundariesTuple, N-1 >( boundaryConditions, buffer, flag, x, y, z );
}



template< typename FlagField_T, typename Stencil, typename... Boundaries >
template< typename BoundariesTuple, int N >
typename std::enable_if<(N!=-1), void>::type BoundaryHandling< FlagField_T, Stencil, Boundaries... >::getBoundaryConditions( const BoundariesTuple & boundaryConditions,
                                                                             std::vector< std::string > & bcs ) const
{
   const auto & boundaryCondition = std::get<N>( boundaryConditions );

   std::ostringstream oss;
   oss << boundaryCondition.getUID().getIdentifier() << " (";

   std::vector< FlagUID > uids;
   boundaryCondition.pushFlags( uids );
   for( auto uid = uids.begin(); uid != uids.end(); ++uid )
   {
      oss << uid->getIdentifier();
      auto next = uid;
      if( ++next != uids.end() ) oss << " & ";
   }

   oss << " => ";
   valueToStream( oss, boundaryCondition.getMask() );
   oss << ")";

   bcs.push_back( oss.str() );

   getBoundaryConditions< BoundariesTuple, N-1 >( boundaryConditions, bcs );
}



} // namespace boundary



using boundary::BoundaryHandling;

template< typename FlagField_T, typename Stencil, typename... Boundaries >
inline std::ostream & operator<<( std::ostream & os, const BoundaryHandling< FlagField_T, Stencil, Boundaries... > & boundaryHandling ) {

   boundaryHandling.toStream( os );
   return os;
}



} // namespace walberla
