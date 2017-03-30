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
//! \file CommonSchemeFunctions.h
//! \ingroup gather
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Common Functions for Gather Schemes
//
//======================================================================================================================

#pragma once

#include "GatherPackInfo.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"
#include "domain_decomposition/BlockStorage.h"

#include <vector>


/// \cond internal

namespace walberla {
namespace gather {
namespace internal {


   inline void packData(  BlockStorage & blocks,
                          std::vector<shared_ptr<GatherPackInfo> > & packInfos,
                          mpi::SendBuffer & sendBuffer )
   {
      // Loop over all blocks and pack data
      for ( auto blockIter = blocks.begin(); blockIter != blocks.end(); ++ blockIter )
         for( size_t  s =0; s < packInfos.size() ; ++s )
         {
            // Packed Message Example:
            // | 0 | PackInfo0_data | 2 | 0 | PackInfo2_data | 1 | PackInfo1_data |
            // |        Block1                               |    Block2          |  ...
            // On Block1 PackInfo1 does not send anything - since this can happen also the PackInfo Nr has to be
            // included in the message
            // the information about the sending block is lost, so the second line is not available at the unpacking stage
            size_t bufferPos = sendBuffer.size();
            sendBuffer << s;

            size_t bufferSizeBefore = sendBuffer.size();
            packInfos[s]->packData( &( *blockIter ), sendBuffer );
            if( sendBuffer.size() == bufferSizeBefore )
               sendBuffer.rewind( bufferPos ); // undo the "buffer<<s"
         }
   }


   inline void unpackData(  std::vector<shared_ptr<GatherPackInfo> > & packInfos,
                            mpi::RecvBuffer & recvBuffer )
   {
      while ( !recvBuffer.isEmpty() )
      {
         size_t s;
         recvBuffer >> s;
         packInfos[s]->unpackData( recvBuffer );
      }

   }

} // namespace internal
} // namespace gather
} // namespace walberla

/// \endcond

