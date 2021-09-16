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
//! \file MPIExport.cpp
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "python_coupling/PythonWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "python_coupling/helper/PythonIterableToStdVector.h"

#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"
#include "core/mpi/Gather.h"
#include "core/mpi/Broadcast.h"

#include <vector>
#include "pybind11/stl.h"

namespace py = pybind11;




namespace walberla {
namespace python_coupling {

   typedef std::vector<int64_t>      IntStdVector;
   typedef std::vector<real_t>       RealStdVector;
   typedef std::vector<std::string>  StringStdVector;


   //===================================================================================================================
   //
   //  MPIManager
   //
   //===================================================================================================================


   static int  rank()              { return MPIManager::instance()->rank();               }
   static int  worldRank()         { return MPIManager::instance()->worldRank();          }
   static int  numProcesses()      { return MPIManager::instance()->numProcesses();       }
   static bool hasCartesianSetup() { return MPIManager::instance()->hasCartesianSetup();  }
   static bool rankValid()         { return MPIManager::instance()->rankValid();          }



   //===================================================================================================================
   //
   //  Broadcast
   //
   //===================================================================================================================

   static py::object broadcast_string( py::object value, int sendRank ) //NOLINT
   {
      if ( py::isinstance<std::string>(value) )
      {
         std::string extractedValue = py::cast< std::string >(value);
         mpi::broadcastObject( extractedValue , sendRank );
         return py::cast( extractedValue );
      }
      StringStdVector extractedValue = pythonIterableToStdVector< StringStdVector::value_type >( value );
      mpi::broadcastObject( extractedValue, sendRank );
      return py::cast( extractedValue );
   }

   static py::object broadcast_int( py::object value, int sendRank ) //NOLINT
   {
      if ( py::isinstance<int64_t>(value) )
      {
         int64_t extractedValue = py::cast< int64_t >(value);
         mpi::broadcastObject( extractedValue , sendRank );
         return py::cast( extractedValue );
      }
      IntStdVector extractedValue = pythonIterableToStdVector< IntStdVector::value_type >( value );
      mpi::broadcastObject( extractedValue, sendRank );
      return py::cast( extractedValue );
   }

   static py::object broadcast_real( py::object value, int sendRank ) //NOLINT
   {
      if ( py::isinstance<real_t>(value) )
      {
         real_t extractedValue = py::cast< real_t  >(value);
         mpi::broadcastObject( extractedValue , sendRank);
         return py::cast( extractedValue );
      }
      RealStdVector extractedValue = pythonIterableToStdVector< RealStdVector::value_type >( value );
      mpi::broadcastObject( extractedValue , sendRank);
      return py::cast( extractedValue );
   }


   //===================================================================================================================
   //
   //  Reduce
   //
   //===================================================================================================================


   static py::object reduce_int( py::object value, mpi::Operation op, int recvRank ) //NOLINT
   {
      if ( py::isinstance<int64_t>(value) )
      {
         int64_t extractedValue = py::cast< int64_t >(value);
         mpi::reduceInplace( extractedValue , op, recvRank );
         return py::cast( extractedValue );
      }
      IntStdVector extractedValue = pythonIterableToStdVector< IntStdVector::value_type >( value );
      mpi::reduceInplace( extractedValue, op, recvRank );
      return py::cast( extractedValue );
   }

   static py::object reduce_real( py::object value, mpi::Operation op, int recvRank ) //NOLINT
   {
      if ( py::isinstance<real_t>(value) )
      {
         real_t extractedValue = py::cast< real_t  >(value);
         mpi::reduceInplace( extractedValue , op, recvRank);
         return py::cast( extractedValue );
      }
      RealStdVector extractedValue = pythonIterableToStdVector< RealStdVector::value_type >( value );
      mpi::reduceInplace( extractedValue , op, recvRank);
      return py::cast( extractedValue );
   }


   static py::object allreduce_int( py::object value, mpi::Operation op ) //NOLINT
   {
      if ( py::isinstance<int64_t>(value) )
      {
         int64_t extractedValue = py::cast< int64_t >(value);
         mpi::allReduceInplace( extractedValue , op );
         return py::cast( extractedValue );
      }
      IntStdVector extractedValue = pythonIterableToStdVector< IntStdVector::value_type >( value );
      mpi::allReduceInplace( extractedValue, op );
      return py::cast( extractedValue );
   }

   static py::object allreduce_real( py::object value, mpi::Operation op ) //NOLINT
   {
      if ( py::isinstance<real_t>(value) )
      {
         real_t extractedValue = py::cast< real_t  >(value);
         mpi::allReduceInplace( extractedValue , op );
         return py::cast( extractedValue );
      }
      RealStdVector extractedValue = pythonIterableToStdVector< RealStdVector::value_type >( value );
      mpi::allReduceInplace( extractedValue , op );
      return py::cast( extractedValue );
   }


   //===================================================================================================================
   //
   //  Gather
   //
   //===================================================================================================================

   static IntStdVector gather_int( py::object value, int recvRank ) //NOLINT
   {
      if ( ! py::isinstance<int64_t>(value) )
      {
         throw py::cast_error("Could not gather the given value - unknown type");
      }
      int64_t extractedValue = py::cast< int64_t >(value);
      return mpi::gather( extractedValue , recvRank );
   }

   static RealStdVector gather_real( py::object value, int recvRank ) //NOLINT
   {
      if ( ! py::isinstance<real_t>(value) )
      {
         throw py::cast_error("Could not gather the given value - unknown type");
      }
      real_t extractedValue = py::cast< real_t  >(value);
      return mpi::gather( extractedValue , recvRank);
   }


   static IntStdVector allgather_int( py::object value ) //NOLINT
   {
      if ( ! py::isinstance<int64_t>(value) )
      {
         throw py::cast_error("Could not gather the given value - unknown type");
      }
      int64_t extractedValue = py::cast< int64_t >(value);
      return mpi::allGather( extractedValue );
   }

   static RealStdVector allgather_real( py::object value ) //NOLINT
   {
      if ( ! py::isinstance<real_t>(value) )
      {
         throw py::cast_error("Could not gather the given value - unknown type");
      }
      real_t extractedValue = py::cast< real_t  >(value);
      return mpi::allGather( extractedValue );
   }



   //===================================================================================================================
   //
   //  Export
   //
   //===================================================================================================================

   static void worldBarrier()
   {
      WALBERLA_MPI_WORLD_BARRIER();
   }


   void exportMPI(py::module_ &m)
   {
      py::module_ m2 = m.def_submodule("mpi", "MPI Extension of the waLBerla python bindings");

      m2.def( "rank"             , &rank             );
      m2.def( "worldRank"        , &worldRank        );
      m2.def( "numProcesses"     , &numProcesses     );
      m2.def( "hasCartesianSetup", &hasCartesianSetup);
      m2.def( "rankValid"        , &rankValid        );
      m2.def( "worldBarrier"     , &worldBarrier     );

      py::enum_<mpi::Operation>(m2, "Operation")
              .value("MIN"    ,      mpi::MIN )
              .value("MAX"    ,      mpi::MAX )
              .value("SUM"    ,      mpi::SUM )
              .value("PRODUCT",      mpi::PRODUCT )
              .value("LOGICAL_AND",  mpi::LOGICAL_AND )
              .value("BITWISE_AND",  mpi::BITWISE_AND )
              .value("LOGICAL_OR",   mpi::LOGICAL_OR  )
              .value("BITWISE_OR",   mpi::BITWISE_OR  )
              .value("LOGICAL_XOR",  mpi::LOGICAL_XOR )
              .value("BITWISE_XOR",  mpi::BITWISE_XOR )
              .export_values();

      m2.def( "broadcastInt",   &broadcast_int);
      m2.def( "broadcastReal",  &broadcast_real);
      m2.def( "broadcastString",&broadcast_string);

      m2.def( "reduceInt",     &reduce_int);
      m2.def( "reduceReal",    &reduce_real);
      m2.def( "allreduceInt",  &allreduce_int  );
      m2.def( "allreduceReal", &allreduce_real );

      m2.def( "gatherInt",     &gather_int);
      m2.def( "gatherReal",    &gather_real);
      m2.def( "allgatherInt",  &allgather_int  );
      m2.def( "allgatherReal", &allgather_real );
   }



} // namespace python_coupling
} // namespace walberla


#endif
