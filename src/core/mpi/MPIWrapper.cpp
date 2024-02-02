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
//! \file MPIWrapper.cpp
//! \ingroup core
//! \author Michael Zikeli <michael.zikeli@fau.de>
//
//======================================================================================================================

#include "MPIWrapper.h"

#include <set>

#include "MPIManager.h"

#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
namespace walberla
{

namespace mpi
{
namespace
{
/// These functions than can be used by self defined mpi operations, e.g. by using CustomMPIOperation.
using float16_t = walberla::float16;
// The signature of MPI_User_function looks like this
// typedef void MPI_User_function( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype);

void sum(void* mpiRHSArray, void* mpiLHSArray, int* len, MPI_Datatype*)
{
   // cast mpi type to target c++ type
   auto* rhs = (float16_t*) mpiRHSArray;
   auto* lhs = (float16_t*) mpiLHSArray;
   for (int i = 0; i < *len; ++i)
   {
      *lhs += *rhs;
   }
}

void min(void* mpiRHSArray, void* mpiLHSArray, int* len, MPI_Datatype*)
{
   // cast mpi type to target c++ type
   auto* rhs = (float16_t*) mpiRHSArray;
   auto* lhs = (float16_t*) mpiLHSArray;
   for (int i = 0; i < *len; ++i)
   {
      *lhs = (*rhs >= *lhs) ? *lhs : *rhs;
   }
}

void max(void* mpiRHSArray, void* mpiLHSArray, int* len, MPI_Datatype*)
{
   // cast mpi type to target c++ type
   auto* rhs = (float16_t*) mpiRHSArray;
   auto* lhs = (float16_t*) mpiLHSArray;
   for (int i = 0; i < *len; ++i)
   {
      *lhs = (*rhs <= *lhs) ? *lhs : *rhs;
   }
}

MPI_User_function* returnMPIUserFctPointer(const Operation op)
{
   switch (op)
   {
   case SUM:
      return &sum;
   case MIN:
      return &min;
   case MAX:
      return &max;
   default:
      WALBERLA_ABORT("The chosen operation " << typeid(op).name() << " is not implemented for float16 yet.");
   }
}

}
}

/// Here some MPI_Datatypes and MPI_Ops are initialized that are not part of the MPI Standard and therefore have to be
/// define yourself. This is done in the MPIManager, since they need to be freed before MPIFinalize is called and this
/// way it is easiest to keep track of them.
///     For more information about this feature compare MR !647 (
///     https://i10git.cs.fau.de/walberla/walberla/-/merge_requests/647 )

/*!
 *  \brief Specialization of MPITrait for float16
 *
 *  The initialization of the self defined MPI_Datatype and MPI_Op is done in the MPIManager so that it can be freed
 * before MPI is finalized.
 */
MPI_Datatype MPITrait< walberla::float16 >::type()
{

#ifdef WALBERLA_BUILD_WITH_MPI
   // Since this type should be created only once, a static variable is used as safeguard.
   static bool initializedType = false;
   if (!initializedType)
   {
      // Since float16 consists of two Bytes, a continuous datatype with size of two byte is created.
      mpi::MPIManager::instance()->commitCustomType< walberla::float16, const int >(2);
      initializedType = true;
   }
   return mpi::MPIManager::instance()->getCustomType< walberla::float16 >();
#else
   return mpistubs::MPI_FLOAT16;
#endif
}

MPI_Op MPITrait< walberla::float16 >::operation(const mpi::Operation& op)
{
   WALBERLA_MPI_SECTION()
   {
      // mpi::Operation is an enum and not an enum class, thus, it is not sufficient to make a just a bool variable as
      // safeguard, since all operations are of type mpi::Operation and only the first one would pass the safeguard.
      // Therefore, a set is created and each operation that is called the first time, will be initialized.
      static std::set< mpi::Operation > operationInitializationRegister;
      const bool needsInitialization = std::get< 1 >(operationInitializationRegister.emplace(op));
      if (needsInitialization)
      {
         mpi::MPIManager::instance()->commitCustomOperation< walberla::float16 >(
            op, mpi::returnMPIUserFctPointer(op));
      }
      return MPIManager::instance()->getCustomOperation< walberla::float16 >(op);
   }
   WALBERLA_NON_MPI_SECTION() { WALBERLA_ABORT("If MPI is not used, a custom operator should never be called."); }
}

} // namespace walberla
#endif // WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
