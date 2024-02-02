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
//! \file MPIOperation.h
//! \ingroup core
//! \author Michael Zikeli <michael.zikeli@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Abort.h"
#include "core/mpi/MPIWrapper.h"

namespace walberla::mpi{

//*******************************************************************************************************************
/*! RAII class for MPI operators that commits and frees them
*
*/
//*******************************************************************************************************************
class MPIOperation
{
 public:
   MPIOperation() = delete;

   explicit MPIOperation( MPI_User_function* fct ) : mpiOperation_( MPI_OP_NULL )
   {
      WALBERLA_MPI_SECTION() {
         if ( MPI_Op_create( fct, true, &mpiOperation_) != MPI_SUCCESS )
         {
            WALBERLA_ABORT("MPI_Op_create for " << typeid(mpiOperation_).name() << " failed." );
         }
      } // WALBERLA_MPI_SECTIION
   }

   ~MPIOperation()
   {
      WALBERLA_MPI_SECTION() { MPI_Op_free( & mpiOperation_ ); }
   }

   operator MPI_Op() const
   {
      return mpiOperation_;
   }

 protected:
   MPI_Op mpiOperation_;
};


} // namespace walberla::mpi