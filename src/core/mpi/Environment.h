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
//! \file Environment.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "MPIManager.h"


namespace walberla {
namespace mpi {


/*******************************************************************************************************************//**
*  RAII Object to initialize and finalize MPI
*
*  Typical use case is to create an mpi::Environment at the beginning of your main() function:
\code
int main( int argc, char ** argv )
{
   mpi::Environment mpiEnv( argc, argv );
}
\endcode
*
* This takes care of calling MPI_Init and MPI_Finalize at the correct places
*
* You probably want to use walberla::Environment instead, which not only handles MPI initialization but
* also initializes other waLBerla singletons!
***********************************************************************************************************************/
class Environment
{
public:
   Environment( int & argc, char ** & argv, bool abortOnException = true ) {
      MPIManager::instance()->initializeMPI( &argc, &argv, abortOnException );
   }
   ~Environment() {
      MPIManager::instance()->finalizeMPI();
   }
};


} // namespace mpi
} // namespace walberla
