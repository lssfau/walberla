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
//! \file Tokenizing.h
//! \ingroup mpi
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/mpi/MPIWrapper.h"

namespace walberla {
namespace mpi {



   //*******************************************************************************************************************
   /*! MPI tokenizing ensures that not more than N processes execute the same code portion simultaneously
   *
   *  Usage:
   *    \code
           mpi::Tokenizing tokenizing ( 10 );
           tokenizing.pre();
           criticalFunction(); // maximal 10 processes execute that function simultaneously
           tokenizing.post()
   *    \endcode
   *
   */
   //*******************************************************************************************************************
   class Tokenizing
   {
   public:

      Tokenizing( uint_t nrOfTokens, MPI_Comm comm = MPI_COMM_WORLD );

      void pre() const;
      void post() const;

   protected:

      uint_t    nrOfTokens_;
      MPI_Comm  communicator_;

      static const int mpiTag = 513; // Token = 84 + 111 + 107 + 101 + 110 = 513
   };



   //*******************************************************************************************************************
   /*! Object that starts tokenizing in constructor and ends tokenizing when going out of scope
   *
   *  Usage:
   *    \code
          {
            mpi::TokenizedScope tokenizing ( 10 );
            criticalFunction(); // maximal 10 processes execute that function simultaneously
          }
   *    \endcode
   */
   //*******************************************************************************************************************
   class TokenizedScope
   {
   public:

      TokenizedScope( uint_t nrOfTokens, MPI_Comm comm = MPI_COMM_WORLD ) : t_ ( nrOfTokens, comm ) { t_.pre(); }
      ~TokenizedScope() { t_.post(); }

   private:

      Tokenizing t_;
   };



} // namespace mpi
} // namespace walberla
