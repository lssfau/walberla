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
//! \file MPIGatherScheme.h
//! \ingroup gather
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Communicating gathered data using MPI
//
//======================================================================================================================

#pragma once

#include "GatherPackInfo.h"
#include "core/debug/Debug.h"
#include "core/mpi/BufferSystem.h"

#include <set>
#include <vector>


namespace walberla {

namespace domain_decomposition {
class BlockStorage;
}

namespace gather {


/**
 * Collects / Gathers data from multiple blocks to a single process
 *
 * \ingroup gather
 *
 * Usage / Restrictions:
 * - each block sends a fixed amount of data to a collection process
 * - this amount has to be the same for all timesteps (i.e. calls of communicate() )
 * - to use the Scheme implement the GatherPackInfo interface or use one of the existing implementations
 * - The collect operation is done every time the communicate() method is called. If the result
 *   is not needed immediately, but only at the end of the simulation consider using the FileCollectorScheme.
 *
 *
 * Implementation:
 * - when communicate() is called first, the amount of data that each process sends is determined
 * - then an MPI_Allgather operation over all processes has to be done, to determine which processes participate
 *   this potentially very expensive operation is done only once in the setup phases
 * - a MPI communicator is created for all processes participating in the gather operation (i.e. that packed
 *   something ), and the amount of data that each process sends is sent to the gathering process
 * - subsequent calls of communicate() use that communicator for an MPI_Gatherv operation
 *
 */
class MPIGatherScheme
{
   public:

      //** Construction & Destruction        ***************************************************************************
      /*! \name Construction & Destruction  */
      //@{
      MPIGatherScheme( domain_decomposition::BlockStorage & blockStorage, int gatherRank = 0, uint_t everyNTimestep = 1 );
      ~MPIGatherScheme();

      //@}

      /**
       * Registering PackInfo's
       *
       * - The ownership of the passed pack info is transferred to the MPICollectorScheme,
       *   i.e. the pack info is deleted by the scheme
       * - when calling addPackInfo after communicate(), the expensive setupPhase has to be done
       *   again -> better first call all addPackInfo() then start calling communicate()
       */
      void addPackInfo( const shared_ptr<GatherPackInfo>  & pi );


      /**
       * Performs the gather operation.
       * Collects all data on sending processes according to information given in the pack infos
       * and unpacks them on the process that holds the root block (1,1,1)
       */
      void communicate ();

      /// Similar to communicate but only executes everyNTimestep ( see constructor )
      void operator() () {
         static uint_t timestep = 0;
         WALBERLA_ASSERT_UNEQUAL( everyNTimestep_, 0 );
         if ( timestep % everyNTimestep_ == 0)
            communicate();
         timestep++;
      }


   private:
      using PackInfoVector = std::vector<shared_ptr<GatherPackInfo>>;

      domain_decomposition::BlockStorage  & blocks_;
      PackInfoVector                        packInfos_;               ///< all registered PackInfos
      int                                   gatherRankInGlobalComm_;  ///< gather processes rank in mpiManager->comm()



      void runSetupPhase( );
      void setupGatherCommunicator( bool thisProcessParticipates, MPI_Comm & commOut, int & newRank );


      //** Communication parameters        *****************************************************************************
      /*! \name Communication parameters initialized by runSetupPhase */
      //@{

      bool             setupPhaseDone_;      ///< true if runSetupPhase() was called
      MPI_Comm         gatherCommunicator_;  ///< communicator containing only participating processes
      int              gatherRank_;          ///< rank in gatherCommunicator_ that gathers the data
      std::vector<int> displacementVector_;  ///< encodes how much each participating process sends (see MPI_Gatherv) ( only on gather process )
      std::vector<int> sendBytesPerProcess_; ///< For each process in gatherCommunicator_ the number of bytes to send ( only on gather process )
      int              gatherMsgSize_;       ///< total size of gather message ( only on gather process )
      int              bytesToSend_;         ///< number of bytes sent by this process ( on all processes )

      uint_t           everyNTimestep_;
      //@}
      //****************************************************************************************************************
};

} // namespace gather
} // namespace walberla


