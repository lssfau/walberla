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
//! \file Tokenizing.cpp
//! \ingroup mpi
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "Tokenizing.h"
#include "core/mpi/MPIManager.h"


namespace walberla {
namespace mpi {


   Tokenizing::Tokenizing( uint_t nrOfTokens, MPI_Comm comm )
      : nrOfTokens_ ( nrOfTokens ), communicator_ ( comm )
   {}

   void Tokenizing::pre() const
   {
      WALBERLA_MPI_SECTION()
      {
         if( nrOfTokens_ > uint_t(0) )
         {
            int rank = 0;
            MPI_Comm_rank( communicator_, &rank );

            WALBERLA_ASSERT_GREATER( nrOfTokens_, 0 );

            if ( rank >= int_c(nrOfTokens_) )
            {
               int token = 23;
               MPI_Recv( static_cast< void* >( &token ), 1, MPITrait<int>::type(), rank - int_c(nrOfTokens_),
                         mpiTag, communicator_, MPI_STATUS_IGNORE );
            }
         }
      }
   }

   void Tokenizing::post() const
   {
      WALBERLA_MPI_SECTION()
      {
         if( nrOfTokens_ > uint_t(0) )
         {
            int rank  = 0;
            int ranks = 0;

            MPI_Comm_rank( communicator_, &rank  );
            MPI_Comm_size( communicator_, &ranks );

            WALBERLA_ASSERT_GREATER( nrOfTokens_, 0 );

            if ( rank + int_c(nrOfTokens_) < ranks )
            {
               int token = 42;
               MPI_Send( static_cast< void* >( &token ), 1, MPITrait<int>::type(), rank + int_c(nrOfTokens_),
                         mpiTag, communicator_ );
            }
         }
      }
   }


} // namespace mpi
} // namespace walberla
