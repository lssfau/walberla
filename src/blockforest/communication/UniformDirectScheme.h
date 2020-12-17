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
//! \file UniformDirectScheme.h
//! \ingroup blockforest
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"
#include "core/Set.h"
#include "core/mpi/Datatype.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/MPIWrapper.h"
#include "core/uid/SUID.h"

#include <vector>
#include <map>
#include <functional>

#include "communication/UniformMPIDatatypeInfo.h"

namespace walberla {
namespace blockforest {
namespace communication {

//*******************************************************************************************************************
/*! Communication for a single field using MPI datatypes ( no intermediate buffers )
*
*/
//*******************************************************************************************************************
template< typename Stencil_T >
class UniformDirectScheme
{
public:
   typedef Stencil_T Stencil;
   typedef walberla::communication::UniformMPIDatatypeInfo UniformMPIDatatypeInfo;
   typedef walberla::communication::UniformMPIDatatypeInfo CommunicationItemInfo;

   //**Construction & Destruction***************************************************************************************
   /*! \name Construction & Destruction */
   //@{
   explicit UniformDirectScheme( const weak_ptr<StructuredBlockForest> & bf,
                                 const shared_ptr<UniformMPIDatatypeInfo> & dataInfo = shared_ptr<UniformMPIDatatypeInfo>(),
                                 const int tag = 778 ) // waLBerla = 119+97+76+66+101+114+108+97
      : blockForest_( bf ),
        setupRequired_( true ),
        communicationRunning_( false ),
        requiredBlockSelectors_( Set<SUID>::emptySet() ),
        incompatibleBlockSelectors_( Set<SUID>::emptySet() ),
        tag_( tag )
   {
      if ( dataInfo )
         dataInfos_.push_back( dataInfo );
   }

   UniformDirectScheme( const weak_ptr<StructuredBlockForest> & bf,
                        const Set<SUID> & requiredBlockSelectors,
                        const Set<SUID> & incompatibleBlockSelectors,
                        const shared_ptr<UniformMPIDatatypeInfo> & dataInfo = shared_ptr<UniformMPIDatatypeInfo>(),
                        const int tag = 778 ) // waLBerla = 119+97+76+66+101+114+108+97
      : blockForest_( bf ),
        setupRequired_( true ),
        communicationRunning_( false ),
        requiredBlockSelectors_( requiredBlockSelectors ),
        incompatibleBlockSelectors_( incompatibleBlockSelectors ),
        tag_( tag )
   {
      if ( dataInfo )
         dataInfos_.push_back( dataInfo );
   }

   //@}
   //*******************************************************************************************************************



   //** Registration of data to communicate ****************************************************************************
   /*! \name Registration of data to communicate */
   //@{
   void addDataToCommunicate(  const shared_ptr<UniformMPIDatatypeInfo> & dataInfo );
   //@}
   //*******************************************************************************************************************


   //** Synchronous Communication **************************************************************************************
   /*! \name Synchronous Communication */
   //@{
   inline void operator() () { communicate(); }
   inline void communicate();
   //@}
   //*******************************************************************************************************************


   //** Asynchronous Communication *************************************************************************************
   /*! \name Asynchronous Communication */
   //@{
   void startCommunication();
   void wait();

   std::function<void()> getStartCommunicateFunctor() { return std::bind( &UniformDirectScheme::startCommunication, this ); }
   std::function<void()> getWaitFunctor()             { return std::bind( &UniformDirectScheme::wait,               this ); }
   //@}
   //*******************************************************************************************************************

protected:
   void setup();

   struct CommInfo
   {
      uint_t             dataIdx;       ///< index into dataInfos_
      BlockID            localBlockId;
      BlockID            remoteBlockId;
      uint_t             remoteProcess;
      stencil::Direction dir;

      static bool sortByLocal( const CommInfo & lhs, const CommInfo & rhs );
      static bool sortByRemote( const CommInfo & lhs, const CommInfo & rhs );
      inline friend std::ostream & operator<<( std::ostream & os, const CommInfo & ci )
      {
         os << '{' << ci.localBlockId.getTreeId() << ',' << ci.remoteBlockId.getTreeId() << ','
                   << ci.remoteProcess << ',' << stencil::dirToString[ci.dir] << '}';
         return os;
      }
   };

   weak_ptr<StructuredBlockForest> blockForest_;

   bool setupRequired_;          //< this is set in the beginning or when new communication item was added
   bool communicationRunning_;   //< this is true between startCommunication() and wait()

   Set<SUID> requiredBlockSelectors_;
   Set<SUID> incompatibleBlockSelectors_;

   std::vector<CommInfo> sendInfos_;
   std::vector<CommInfo> recvInfos_;

   std::vector< MPI_Request  >              mpiRequests_;
   std::vector< shared_ptr<mpi::Datatype> > mpiDatatypes_;

   std::vector< shared_ptr<UniformMPIDatatypeInfo> > dataInfos_;

   int tag_;

}; // class UniformDirectScheme



} // namespace communication
} // namespace blockforest
} // namespace walberla

#include "UniformDirectScheme.impl.h"
