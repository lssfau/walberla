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
//! \file ModuleInit.cpp
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "ModuleInit.h"
#include "waLBerlaDefinitions.h"
#include "core/mpi/MPIManager.h"
#include "core/Abort.h"
#include "core/debug/CheckFunctions.h"

// Workaround for OpenMPI library: it dynamically loads plugins which causes trouble when walberla itself is a shared lib
#if defined(OPEN_MPI) && !defined(_WIN32) && OMPI_MAJOR_VERSION < 3
#define OPEN_MPI_WORKAROUND
#include <dlfcn.h>
#endif



namespace walberla {
namespace python_coupling {

void initWalberlaForPythonModule()
{

#ifdef OPEN_MPI_WORKAROUND
   int mode = RTLD_NOW | RTLD_GLOBAL;
#ifdef RTLD_NOLOAD
  mode |= RTLD_NOLOAD;
#endif
   void *symbol = dlsym(RTLD_DEFAULT, "MPI_Init");
   WALBERLA_CHECK_NOT_NULLPTR(symbol, "Unable to find OpenMPI symbol");
   Dl_info symbol_info;
   dladdr(symbol, &symbol_info);
   void * handle = dlopen(symbol_info.dli_fname, mode);
   WALBERLA_CHECK_NOT_NULLPTR(handle, "Unable to load libmpi into global symbol space");
#endif

   Abort::instance()->resetAbortFunction( Abort::exceptionAbort );

   mpi::MPIManager::instance()->initializeMPI( nullptr,nullptr );
   mpi::MPIManager::instance()->useWorldComm();
}


} // namespace python_coupling
} // namespace walberla


