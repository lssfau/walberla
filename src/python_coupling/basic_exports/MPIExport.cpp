#include "python_coupling/PythonWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "python_coupling/helper/PythonIterableToStdVector.h"

#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"
#include "core/mpi/Gather.h"
#include "core/mpi/Broadcast.h"

#include <vector>

using namespace boost::python;




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

   static object broadcast_string( object value, int sendRank ) //NOLINT
   {
      if ( extract<std::string>(value).check() )
      {
         std::string extractedValue = extract< std::string >(value);
         mpi::broadcastObject( extractedValue , sendRank );
         return object( extractedValue );
      }
      StringStdVector extractedValue = pythonIterableToStdVector< StringStdVector::value_type >( value );
      mpi::broadcastObject( extractedValue, sendRank );
      return object( extractedValue );
   }

   static object broadcast_int( object value, int sendRank ) //NOLINT
   {
      if ( extract<int64_t>(value).check() )
      {
         int64_t extractedValue = extract< int64_t >(value);
         mpi::broadcastObject( extractedValue , sendRank );
         return object( extractedValue );
      }
      IntStdVector extractedValue = pythonIterableToStdVector< IntStdVector::value_type >( value );
      mpi::broadcastObject( extractedValue, sendRank );
      return object( extractedValue );
   }

   static object broadcast_real( object value, int sendRank ) //NOLINT
   {
      if ( extract<real_t>(value).check() )
      {
         real_t extractedValue = extract< real_t  >(value);
         mpi::broadcastObject( extractedValue , sendRank);
         return object( extractedValue );
      }
      RealStdVector extractedValue = pythonIterableToStdVector< RealStdVector::value_type >( value );
      mpi::broadcastObject( extractedValue , sendRank);
      return object( extractedValue );
   }


   //===================================================================================================================
   //
   //  Reduce
   //
   //===================================================================================================================


   static object reduce_int( object value, mpi::Operation op, int recvRank ) //NOLINT
   {
      if ( extract<int64_t>(value).check() )
      {
         int64_t extractedValue = extract< int64_t >(value);
         mpi::reduceInplace( extractedValue , op, recvRank );
         return object( extractedValue );
      }
      IntStdVector extractedValue = pythonIterableToStdVector< IntStdVector::value_type >( value );
      mpi::reduceInplace( extractedValue, op, recvRank );
      return object( extractedValue );
   }

   static object reduce_real( object value, mpi::Operation op, int recvRank ) //NOLINT
   {
      if ( extract<real_t>(value).check() )
      {
         real_t extractedValue = extract< real_t  >(value);
         mpi::reduceInplace( extractedValue , op, recvRank);
         return object( extractedValue );
      }
      RealStdVector extractedValue = pythonIterableToStdVector< RealStdVector::value_type >( value );
      mpi::reduceInplace( extractedValue , op, recvRank);
      return object( extractedValue );
   }


   static object allreduce_int( object value, mpi::Operation op ) //NOLINT
   {
      if ( extract<int64_t>(value).check() )
      {
         int64_t extractedValue = extract< int64_t >(value);
         mpi::allReduceInplace( extractedValue , op );
         return object( extractedValue );
      }
      IntStdVector extractedValue = pythonIterableToStdVector< IntStdVector::value_type >( value );
      mpi::allReduceInplace( extractedValue, op );
      return object( extractedValue );
   }

   static object allreduce_real( object value, mpi::Operation op ) //NOLINT
   {
      if ( extract<real_t>(value).check() )
      {
         real_t extractedValue = extract< real_t  >(value);
         mpi::allReduceInplace( extractedValue , op );
         return object( extractedValue );
      }
      RealStdVector extractedValue = pythonIterableToStdVector< RealStdVector::value_type >( value );
      mpi::allReduceInplace( extractedValue , op );
      return object( extractedValue );
   }


   //===================================================================================================================
   //
   //  Gather
   //
   //===================================================================================================================

   static IntStdVector gather_int( object value, int recvRank ) //NOLINT
   {
      if ( ! extract<int64_t>(value).check() )
      {
         PyErr_SetString( PyExc_RuntimeError, "Could not gather the given value - unknown type");
         throw error_already_set();
      }
      int64_t extractedValue = extract< int64_t >(value);
      return mpi::gather( extractedValue , recvRank );
   }

   static RealStdVector gather_real( object value, int recvRank ) //NOLINT
   {
      if ( ! extract<real_t>(value).check() )
      {
         PyErr_SetString( PyExc_RuntimeError, "Could not gather the given value - unknown type");
         throw error_already_set();
      }
      real_t extractedValue = extract< real_t  >(value);
      return mpi::gather( extractedValue , recvRank);
   }


   static IntStdVector allgather_int( object value ) //NOLINT
   {
      if ( ! extract<int64_t>(value).check() )
      {
         PyErr_SetString( PyExc_RuntimeError, "Could not gather the given value - unknown type");
         throw error_already_set();
      }
      int64_t extractedValue = extract< int64_t >(value);
      return mpi::allGather( extractedValue );
   }

   static RealStdVector allgather_real( object value ) //NOLINT
   {
      if ( ! extract<real_t>(value).check() )
      {
         PyErr_SetString( PyExc_RuntimeError, "Could not gather the given value - unknown type");
         throw error_already_set();
      }
      real_t extractedValue = extract< real_t  >(value);
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


   void exportMPI()
   {
      object mpiModule( handle<>( borrowed(PyImport_AddModule("walberla_cpp.mpi"))));
      scope().attr("mpi") = mpiModule;
      scope mpiScope = mpiModule;

      def( "rank"             , &rank             );
      def( "worldRank"        , &worldRank        );
      def( "numProcesses"     , &numProcesses     );
      def( "hasCartesianSetup", &hasCartesianSetup);
      def( "rankValid"        , &rankValid        );
      def( "worldBarrier"     , &worldBarrier     );

      enum_<mpi::Operation>("Operation")
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

      def( "broadcastInt",   &broadcast_int,    ( arg("sendRank") = 0) );
      def( "broadcastReal",  &broadcast_real,   ( arg("sendRank") = 0) );
      def( "broadcastString",&broadcast_string, ( arg("sendRank") = 0) );

      def( "reduceInt",     &reduce_int,     ( arg("recvRank") = 0 ) );
      def( "reduceReal",    &reduce_real,    ( arg("recvRank") = 0 ) );
      def( "allreduceInt",  &allreduce_int  );
      def( "allreduceReal", &allreduce_real );


      class_< IntStdVector>   ("IntStdVector") .def(vector_indexing_suite<IntStdVector>()  );
      class_< RealStdVector>  ("RealStdVector").def(vector_indexing_suite<RealStdVector>() );
      class_< StringStdVector>("StringStdVector").def(vector_indexing_suite<StringStdVector>() );

      def( "gatherInt",     &gather_int ,   ( arg("recvRank") = 0 ) );
      def( "gatherReal",    &gather_real,   ( arg("recvRank") = 0 ) );
      def( "allgatherInt",  &allgather_int  );
      def( "allgatherReal", &allgather_real );
   }



} // namespace python_coupling
} // namespace walberla


#endif
