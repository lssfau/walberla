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
//! \file BasicExports.cpp
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

// Do not reorder includes - the include order is important
#include "python_coupling/PythonWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "BasicExports.h"
#include "MPIExport.h"
#include "python_coupling/helper/ModuleScope.h"
#include "core/waLBerlaBuildInfo.h"
#include "core/logging/Logging.h"
#include "core/Abort.h"
#include "core/cell/CellInterval.h"
#include "core/math/AABB.h"
#include "core/mpi/MPIIO.h"
#include "core/timing/ReduceType.h"
#include "core/timing/TimingPool.h"
#include "core/timing/TimingTree.h"
#include "communication/UniformPackInfo.h"
#include "communication/UniformMPIDatatypeInfo.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "python_coupling/Manager.h"
#include "python_coupling/helper/BlockStorageExportHelpers.h"
#include "stencil/Directions.h"

#include <boost/version.hpp>

#include <functional>

using namespace boost::python;


namespace walberla {
namespace python_coupling {


template <class T>
struct NumpyIntConversion
{
    NumpyIntConversion()
    {
        converter::registry::push_back( &convertible, &construct, boost::python::type_id<T>() );
    }

    static void* convertible( PyObject* pyObj)
    {
       auto typeName = std::string( Py_TYPE(pyObj)->tp_name );
       if ( typeName.substr(0,9) == "numpy.int" )
          return pyObj;
       return nullptr;
    }

    static void construct( PyObject* pyObj, converter::rvalue_from_python_stage1_data* data )
    {
        handle<> x(borrowed(pyObj));
        object o(x);
        T value = extract<T>(o.attr("__int__")());
        void* storage =( (boost::python::converter::rvalue_from_python_storage<T>*) data)->storage.bytes;
        new (storage) T(value);
        data->convertible = storage;
    }
};

template <class T>
struct NumpyFloatConversion
{
    NumpyFloatConversion()
    {
        converter::registry::push_back( &convertible, &construct, boost::python::type_id<T>() );
    }

    static void* convertible(PyObject* pyObj)
    {
       auto typeName = std::string( Py_TYPE(pyObj)->tp_name );
       if ( typeName.substr(0,11) == "numpy.float" )
          return pyObj;
        return nullptr;
    }

    static void construct(PyObject* pyObj, converter::rvalue_from_python_stage1_data* data)
    {
        handle<> x(borrowed(pyObj));
        object o(x);
#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
        T value = extract<T>(o.attr("__float__")());
#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif
        void* storage =( (boost::python::converter::rvalue_from_python_storage<T>*) data)->storage.bytes;
        new (storage) T(value);
        data->convertible = storage;
    }
};


#if BOOST_VERSION < 106300
// taken from https://github.com/boostorg/python/commit/97e4b34a15978ca9d7c296da2de89b78bba4e0d5
template <class T>
struct exportSharedPtr
{
   exportSharedPtr()
   {
   converter::registry::insert( &convertible, &construct, boost::python::type_id<std::shared_ptr<T> >()
#ifndef BOOST_PYTHON_NO_PY_SIGNATURES
                                , &converter::expected_from_python_type_direct<T>::get_pytype
#endif
                              );
   }

private:
   static void* convertible( PyObject* p )
   {
      if ( p == Py_None )
         return p;

      return converter::get_lvalue_from_python( p, converter::registered<T>::converters );
   }

   static void construct( PyObject* source, converter::rvalue_from_python_stage1_data* data )
   {
      void* const storage = ( (converter::rvalue_from_python_storage< std::shared_ptr<T> >*) data )->storage.bytes;
      // Deal with the "None" case.
      if ( data->convertible == source )
         new (storage) std::shared_ptr<T>();
      else
      {
         std::shared_ptr<void> hold_convertible_ref_count( (void*)0, converter::shared_ptr_deleter( handle<>( borrowed( source ) ) ) );
         // use aliasing constructor
         new (storage) std::shared_ptr<T>( hold_convertible_ref_count, static_cast<T*>(data->convertible) );
      }

      data->convertible = storage;
   }
};
#endif

//======================================================================================================================
//
//  Helper Functions
//
//======================================================================================================================

void checkForThreeSequence( const object & o, const char * message )
{
   // strange construct required because also the len function throws , if object has no len
   try {
      if ( len(o ) != 3 )  throw error_already_set();
   }
   catch( error_already_set & ) {
      PyErr_SetString(PyExc_RuntimeError, message);
      throw error_already_set();
   }
}

//======================================================================================================================
//
//  Vector3
//
//======================================================================================================================


template<typename T>
struct Vector3_to_PythonTuple
{
   static PyObject* convert( Vector3<T> const& v )
   {
      auto resultTuple = boost::python::make_tuple(v[0], v[1], v[2] );
      return boost::python::incref ( boost::python::object ( resultTuple ).ptr () );
   }
};

template<typename T>
struct PythonTuple_to_Vector3
{
   PythonTuple_to_Vector3()
   {
     boost::python::converter::registry::push_back(
       &convertible,
       &construct,
       boost::python::type_id<Vector3<T> >());
   }

   static void* convertible(PyObject* obj)
   {
      using namespace boost::python;

      if ( ! ( PySequence_Check(obj) && PySequence_Size( obj ) == 3 ))
         return nullptr;

      object element0 ( handle<>( borrowed( PySequence_GetItem(obj,0) )));
      object element1 ( handle<>( borrowed( PySequence_GetItem(obj,1) )));
      object element2 ( handle<>( borrowed( PySequence_GetItem(obj,2) )));

      if (  extract<T>( element0 ).check() &&
            extract<T>( element1 ).check() &&
            extract<T>( element2 ).check() )
         return obj;
      else
         return nullptr;
   }

   static void construct( PyObject* obj, boost::python::converter::rvalue_from_python_stage1_data* data )
   {
      using namespace boost::python;

      object element0 ( handle<>( borrowed( PySequence_GetItem(obj,0) )));
      object element1 ( handle<>( borrowed( PySequence_GetItem(obj,1) )));
      object element2 ( handle<>( borrowed( PySequence_GetItem(obj,2) )));


      // Grab pointer to memory into which to construct the new Vector3
      void* storage = ( (boost::python::converter::rvalue_from_python_storage<Vector3<T> >*) data )->storage.bytes;

      new (storage) Vector3<T> ( extract<T>( element0 ),
                                 extract<T>( element1 ),
                                 extract<T>( element2 ) );

      // Stash the memory chunk pointer for later use by boost.python
      data->convertible = storage;
   }
};

void exportVector3()
{
   // To Python
   boost::python::to_python_converter< Vector3<bool      >, Vector3_to_PythonTuple<bool      > >();

   boost::python::to_python_converter< Vector3<real_t    >, Vector3_to_PythonTuple<real_t    > >();

   boost::python::to_python_converter< Vector3<uint8_t   >, Vector3_to_PythonTuple<uint8_t   > >();
   boost::python::to_python_converter< Vector3<uint16_t  >, Vector3_to_PythonTuple<uint16_t  > >();
   boost::python::to_python_converter< Vector3<uint32_t  >, Vector3_to_PythonTuple<uint32_t  > >();
   boost::python::to_python_converter< Vector3<uint64_t  >, Vector3_to_PythonTuple<uint64_t  > >();

   boost::python::to_python_converter< Vector3<cell_idx_t>, Vector3_to_PythonTuple<cell_idx_t> >();

   // From Python
   PythonTuple_to_Vector3<bool   >();

   PythonTuple_to_Vector3<real_t   >();

   PythonTuple_to_Vector3<uint8_t  >();
   PythonTuple_to_Vector3<uint16_t >();
   PythonTuple_to_Vector3<uint32_t >();
   PythonTuple_to_Vector3<uint64_t >();

   PythonTuple_to_Vector3<cell_idx_t>();
}


//======================================================================================================================
//
//  Cell
//
//======================================================================================================================


struct Cell_to_PythonTuple
{
   static PyObject* convert( Cell const& c )
   {
      auto resultTuple = boost::python::make_tuple(c[0], c[1], c[2] );
      return boost::python::incref ( boost::python::object ( resultTuple ).ptr () );
   }
};

struct PythonTuple_to_Cell
{
   PythonTuple_to_Cell()
   {
     boost::python::converter::registry::push_back(
       &convertible,
       &construct,
       boost::python::type_id< Cell >());
   }

   static void* convertible( PyObject* obj )
   {
      using namespace boost::python;

      if ( ! ( PySequence_Check(obj) && PySequence_Size( obj ) == 3 ))
         return nullptr;

      object element0 ( handle<>( borrowed( PySequence_GetItem(obj,0) )));
      object element1 ( handle<>( borrowed( PySequence_GetItem(obj,1) )));
      object element2 ( handle<>( borrowed( PySequence_GetItem(obj,2) )));

      if (  extract<cell_idx_t>( element0 ).check() &&
            extract<cell_idx_t>( element1 ).check() &&
            extract<cell_idx_t>( element2 ).check() )
         return obj;
      else
         return nullptr;
   }

   static void construct( PyObject* obj, boost::python::converter::rvalue_from_python_stage1_data* data )
   {
      using namespace boost::python;

      object element0 ( handle<>( borrowed( PySequence_GetItem(obj,0) )));
      object element1 ( handle<>( borrowed( PySequence_GetItem(obj,1) )));
      object element2 ( handle<>( borrowed( PySequence_GetItem(obj,2) )));


      // Grab pointer to memory into which to construct the new Vector3
      void* storage = ( (boost::python::converter::rvalue_from_python_storage<Cell>*) data )->storage.bytes;

      new (storage) Cell( extract<cell_idx_t>( element0 ),
                          extract<cell_idx_t>( element1 ),
                          extract<cell_idx_t>( element2 ) );

      // Stash the memory chunk pointer for later use by boost.python
      data->convertible = storage;
   }
};

void exportCell()
{
   // To Python
   boost::python::to_python_converter< Cell, Cell_to_PythonTuple >();
   // From Python
   PythonTuple_to_Cell();
}

//======================================================================================================================
//
//  CellInterval
//
//======================================================================================================================


void cellInterval_setMin( CellInterval & ci, const Cell & min ) {
   ci.min() = min;
}
void cellInterval_setMax( CellInterval & ci, const Cell & max ) {
   ci.max() = max;
}
void cellInterval_shift( CellInterval & ci, cell_idx_t xShift, cell_idx_t yShift, cell_idx_t zShift ) {
   ci.shift( xShift, yShift, zShift );
}

boost::python::tuple cellInterval_size( CellInterval & ci ) {
   return boost::python::make_tuple( ci.xSize(), ci.ySize(), ci.zSize() );
}

CellInterval cellInterval_getIntersection( CellInterval & ci1, CellInterval & ci2 )
{
   CellInterval result ( ci1 );
   result.intersect( ci2 );
   return result;
}

CellInterval cellInterval_getShifted( CellInterval & ci1, cell_idx_t xShift, cell_idx_t yShift, cell_idx_t zShift )
{
   CellInterval result ( ci1 );
   result.shift( xShift, yShift, zShift  );
   return result;
}

CellInterval cellInterval_getExpanded1( CellInterval & ci1, cell_idx_t expandVal )
{
   CellInterval result ( ci1 );
   result.expand( expandVal  );
   return result;
}

CellInterval cellInterval_getExpanded2( CellInterval & ci1, cell_idx_t xExpand, cell_idx_t yExpand, cell_idx_t zExpand )
{
   CellInterval result ( ci1 );
   result.expand( Cell(xExpand, yExpand, zExpand)  );
   return result;
}

void exportCellInterval()
{
   const Cell & ( CellInterval::*p_getMin )( ) const = &CellInterval::min;
   const Cell & ( CellInterval::*p_getMax )( ) const = &CellInterval::max;

   bool ( CellInterval::*p_contains1) ( const Cell         & ) const = &CellInterval::contains;
   bool ( CellInterval::*p_contains2) ( const CellInterval & ) const = &CellInterval::contains;

   void ( CellInterval::*p_expand1) ( const cell_idx_t ) = &CellInterval::expand;
   void ( CellInterval::*p_expand2) ( const Cell &     ) = &CellInterval::expand;
   
   bool          ( CellInterval::*p_overlaps ) ( const CellInterval & ) const = &CellInterval::overlaps;

   class_<CellInterval>("CellInterval")
      .def( init<const Cell&, const Cell&>() )
      .def( init<cell_idx_t, cell_idx_t, cell_idx_t, cell_idx_t, cell_idx_t, cell_idx_t>() )
      .add_property( "min", make_function( p_getMin, return_value_policy<copy_const_reference>() ), &cellInterval_setMin )
      .add_property( "max", make_function( p_getMax, return_value_policy<copy_const_reference>() ), &cellInterval_setMax )
      .add_property( "size", &cellInterval_size  )
      .def( "empty", &CellInterval::empty )
      .def( "positiveIndicesOnly", &CellInterval::positiveIndicesOnly )
      .def( "contains",        p_contains1 )
      .def( "contains",        p_contains2 )
      .def( "overlaps",        p_overlaps )
      .def( "shift",           &cellInterval_shift )
      .def( "getShifted",      &cellInterval_getShifted )
      .def( "expand",          p_expand1 )
      .def( "expand",          p_expand2 )
      .def( "getExpanded",     &cellInterval_getExpanded1 )
      .def( "getExpanded",     &cellInterval_getExpanded2 )
      .def( "intersect",       &CellInterval::intersect )
      .def( "getIntersection", &cellInterval_getIntersection )
      .def("__eq__",           &CellInterval::operator==)
      .def("__ne__",           &CellInterval::operator!=)
      .add_property( "numCells",  &CellInterval::numCells  )
      .def( self_ns::str(self) )
      ;
}

//======================================================================================================================
//
//  AABB
//
//======================================================================================================================


tuple aabb_getMin( const AABB & domainBB ) {
   return boost::python::make_tuple( domainBB.xMin(), domainBB.yMin(), domainBB.zMin() );
}
tuple aabb_getMax( const AABB & domainBB ) {
   return boost::python::make_tuple( domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
}

void aabb_setMin( AABB & domainBB, object min )
{
   checkForThreeSequence(min, "Error assigning minimum of AABB - Sequence of length 3 required" );
   real_t min0 = extract<real_t>( min[0] );
   real_t min1 = extract<real_t>( min[1] );
   real_t min2 = extract<real_t>( min[2] );

   domainBB = AABB( domainBB.max(), AABB::vector_type ( min0, min1, min2 ) );
}

void aabb_setMax( AABB & domainBB, object max )
{
   checkForThreeSequence( max, "Error assigning maximum of AABB - Sequence of length 3 required" );

   real_t max0 = extract<real_t>( max[0] );
   real_t max1 = extract<real_t>( max[1] );
   real_t max2 = extract<real_t>( max[2] );

   domainBB = AABB( domainBB.min(), AABB::vector_type ( max0, max1, max2 ) );
}


void exportAABB()
{
   bool ( AABB::*p_containsBB  )( const AABB & bb           )     const = &AABB::contains;
   bool ( AABB::*p_containsVec )( const Vector3<real_t> & point ) const = &AABB::contains;


   bool ( AABB::*p_containsClosedInterval1 ) ( const Vector3<real_t> & ) const               = &AABB::containsClosedInterval;
   bool ( AABB::*p_containsClosedInterval2 ) ( const Vector3<real_t> &, const real_t ) const = &AABB::containsClosedInterval;

   AABB ( AABB::*p_getExtended1 ) ( const real_t            ) const = &AABB::getExtended;
   AABB ( AABB::*p_getExtended2 ) ( const Vector3<real_t> & ) const = &AABB::getExtended;

   AABB ( AABB::*p_getScaled1 ) ( const real_t            ) const = &AABB::getScaled;
   AABB ( AABB::*p_getScaled2 ) ( const Vector3<real_t> & ) const = &AABB::getScaled;

   AABB ( AABB::*p_getMerged1 ) ( const Vector3<real_t> & ) const = &AABB::getMerged;
   AABB ( AABB::*p_getMerged2 ) ( const AABB            & ) const = &AABB::getMerged;

   bool ( AABB::*p_intersects1  )( const AABB & bb            ) const = &AABB::intersects;
   bool ( AABB::*p_intersects2  )( const AABB & bb, real_t dx ) const = &AABB::intersects;

   bool ( AABB::*p_intersectsClosed1  )( const AABB & bb            ) const = &AABB::intersectsClosedInterval;
   bool ( AABB::*p_intersectsClosed2  )( const AABB & bb, real_t dx ) const = &AABB::intersectsClosedInterval;

   void ( AABB::*p_extend1 ) ( const real_t            ) = &AABB::extend;
   void ( AABB::*p_extend2 ) ( const Vector3<real_t> & ) = &AABB::extend;

   void  ( AABB::*p_scale1 ) ( const real_t            )  = &AABB::scale;
   void  ( AABB::*p_scale2 ) ( const Vector3<real_t> & )  = &AABB::scale;

   void  ( AABB::*p_merge1 ) ( const Vector3<real_t> & )  = &AABB::merge;
   void  ( AABB::*p_merge2 ) ( const AABB            & )  = &AABB::merge;

   real_t  ( AABB::*p_sqDistance1 ) ( const AABB & )            const = &AABB::sqDistance;
   real_t  ( AABB::*p_sqDistance2 ) ( const Vector3<real_t> & ) const = &AABB::sqDistance;

   real_t  ( AABB::*p_sqMaxDistance1 ) ( const AABB & )            const = &AABB::sqMaxDistance;
   real_t  ( AABB::*p_sqMaxDistance2 ) ( const Vector3<real_t> & ) const = &AABB::sqMaxDistance;


   class_<AABB>("AABB")
      .def( init<real_t,real_t,real_t,real_t,real_t,real_t>() )
      .def( init<Vector3<real_t>,Vector3<real_t> >() )
      .def("__eq__",     &walberla::math::operator==<real_t, real_t > )
      .def("__ne__",     &walberla::math::operator!=<real_t, real_t > )
      .add_property( "min",  &aabb_getMin, &aabb_setMin )
      .add_property( "max",  &aabb_getMax, &aabb_setMax )
      .add_property( "size", &AABB::sizes )
      .def( "empty",    &AABB::empty )
      .def( "volume", &AABB::volume )
      .def( "center", &AABB::center )
      .def( "contains", p_containsBB )
      .def( "contains", p_containsVec )
      .def( "containsClosedInterval", p_containsClosedInterval1 )
      .def( "containsClosedInterval", p_containsClosedInterval2 )
      .def( "getExtended", p_getExtended1 )
      .def( "getExtended", p_getExtended2 )
      .def( "getTranslated", &AABB::getTranslated )
      .def( "getScaled", p_getScaled1 )
      .def( "getScaled", p_getScaled2 )
      .def( "getMerged", p_getMerged1 )
      .def( "getMerged", p_getMerged2 )
      .def( "intersects", p_intersects1 )
      .def( "intersects", p_intersects2 )
      .def( "intersectsClosedInterval", p_intersectsClosed1 )
      .def( "intersectsClosedInterval", p_intersectsClosed2 )
      .def( "intersectionVolume", &AABB::intersectionVolume )
      .def( "getIntersection", &AABB::getIntersection )
      .def( "isIdentical", &AABB::isIdentical )
      .def( "isEqual",     &AABB::isEqual )
      .def( "sqDistance",       p_sqDistance1 )
      .def( "sqDistance",       p_sqDistance2 )
      .def( "sqMaxDistance",    p_sqMaxDistance1 )
      .def( "sqMaxDistance",    p_sqMaxDistance2 )
      .def( "sqSignedDistance", &AABB::sqSignedDistance)
      .def( "distance"        , &AABB::distance)
      .def( "signedDistance"  , &AABB::signedDistance)
      .def( "maxDistance"     , &AABB::maxDistance)
      .def( "extend", p_extend1 )
      .def( "extend", p_extend2 )
      .def( "translate", &AABB::translate )
      .def( "scale", p_scale1 )
      .def( "scale", p_scale2 )
      .def( "merge", p_merge1 )
      .def( "merge", p_merge2 )
      .def( "intersect", &AABB::intersect )
      .def( self_ns::str(self) )
      ;
}

//======================================================================================================================
//
//  Timing
//
//======================================================================================================================

dict buildDictFromTimingNode(const WcTimingNode & tn)
{
   dict result;

   result["all"] = tn.timer_;
   for ( auto it = tn.tree_.begin(); it != tn.tree_.end(); ++it)
   {
      if (it->second.tree_.empty())
      {
         result[it->first] = it->second.timer_;
      } else
      {
         result[it->first] = buildDictFromTimingNode(it->second);
      }
   }

   return result;
}

dict buildDictFromTimingTree(const WcTimingTree & tt)
{
   return buildDictFromTimingNode( tt.getRawData() );
}

void timingTreeStopWrapper(WcTimingTree & tt, const std::string& name)
{
   if (!tt.isTimerRunning(name))
   {
      PyErr_SetString( PyExc_ValueError, ("Timer '" + name + "' is currently not running!").c_str() );
      throw error_already_set();
   }
   tt.stop(name);
}

void exportTiming()
{
   class_<WcTimer> ("Timer")
      .def( init<>() )
      .def( "start",  &WcTimer::start )
      .def( "stop",   &WcTimer::end   )
      .def( "reset",  &WcTimer::reset )
      .def( "merge",  &WcTimer::merge )
      .add_property( "counter",      &WcTimer::getCounter   )
      .add_property( "total",        &WcTimer::total        )
      .add_property( "sumOfSquares", &WcTimer::sumOfSquares )
      .add_property( "average",      &WcTimer::average      )
      .add_property( "variance",     &WcTimer::variance     )
      .add_property( "min",          &WcTimer::min          )
      .add_property( "max",          &WcTimer::max          )
      .add_property( "last",         &WcTimer::last         )
      ;


   WcTimer & ( WcTimingPool::*pGetItem ) ( const std::string & ) = &WcTimingPool::operator[];

   {
      scope classScope =
      class_<WcTimingPool, shared_ptr<WcTimingPool> > ("TimingPool")
         .def( init<>() )
         .def( self_ns::str(self) )
         .def( "__getitem__",     pGetItem, return_internal_reference<1>() )
         .def( "__contains__",    &WcTimingPool::timerExists )
         .def( "getReduced",      &WcTimingPool::getReduced,  ( arg("targetRank") = 0) )
         .def( "merge",           &WcTimingPool::merge, ( arg("mergeDuplicates") = true) )
         .def( "clear",           &WcTimingPool::clear )
         .def( "unifyRegisteredTimersAcrossProcesses", &WcTimingPool::unifyRegisteredTimersAcrossProcesses )
         .def( "logResultOnRoot", &WcTimingPool::logResultOnRoot, (arg("unifyRegisteredTimers") = false) )
         .def( self_ns::str(self) )
         ;

      enum_<timing::ReduceType>("ReduceType")
          .value("min"  , timing::REDUCE_MIN)
          .value("avg"  , timing::REDUCE_AVG)
          .value("max"  , timing::REDUCE_MAX)
          .value("total", timing::REDUCE_TOTAL)
          .export_values()
          ;
   }

   const WcTimer & ( WcTimingTree::*pTimingTreeGet ) ( const std::string & ) const = &WcTimingTree::operator[];
   class_<WcTimingTree, shared_ptr<WcTimingTree> > ("TimingTree")
         .def( init<>() )
         .def( "__getitem__",  pTimingTreeGet, return_internal_reference<1>() )
         .def( "start",        &WcTimingTree::start )
         .def( "stop",         &timingTreeStopWrapper )
         .def( "getReduced",   &WcTimingTree::getReduced )
         .def( "toDict",       &buildDictFromTimingTree )
         .def( self_ns::str(self) )
    ;

#if BOOST_VERSION < 106300
   exportSharedPtr<WcTimingTree>();
#endif

}




//======================================================================================================================
//
//  IBlock
//
//======================================================================================================================

boost::python::object IBlock_getData( boost::python::object iblockObject, const std::string & stringID ) //NOLINT
{
   IBlock * block = boost::python::extract<IBlock*>( iblockObject );

   //typedef std::pair< IBlock *, std::string > BlockStringPair;
   //static std::map< BlockStringPair, object > cache;

   //auto blockStringPair = std::make_pair( &block, stringID );
   //auto it = cache.find( blockStringPair );
   //if ( it != cache.end() )
   //   return it->second;

   BlockDataID id = blockDataIDFromString( *block, stringID );

   auto manager = python_coupling::Manager::instance();
   boost::python::object res =  manager->pythonObjectFromBlockData( *block, id );

   if ( res == boost::python::object() )
      throw BlockDataNotConvertible();

   boost::python::objects::make_nurse_and_patient( res.ptr(), iblockObject.ptr() );

   // write result to cache
   //cache[blockStringPair] = res;
   //TODO cache has bugs when cache is destroyed, probably since objects are freed after py_finalize is called
   //move cache to Manager?

   return res;
}


boost::python::list IBlock_blockDataList( boost::python::object iblockObject ) //NOLINT
{
   IBlock * block = boost::python::extract<IBlock*>( iblockObject );

   const std::vector<std::string> & stringIds = block->getBlockStorage().getBlockDataIdentifiers();

   boost::python::list resultList;

   for( auto it = stringIds.begin(); it != stringIds.end(); ++it ) {
      try {
         resultList.append( boost::python::make_tuple( *it, IBlock_getData( iblockObject, *it) ) );
      }
      catch( BlockDataNotConvertible & /*e*/ ) {
      }
   }

   return resultList;
}

boost::python::object IBlock_iter(  boost::python::object iblockObject )
{
   boost::python::list resultList = IBlock_blockDataList( iblockObject ); //NOLINT
   return resultList.attr("__iter__");
}

boost::python::tuple IBlock_atDomainMinBorder( IBlock & block )
{
   return boost::python::make_tuple( block.getBlockStorage().atDomainXMinBorder(block),
                                     block.getBlockStorage().atDomainYMinBorder(block),
                                     block.getBlockStorage().atDomainZMinBorder(block) );
}

boost::python::tuple IBlock_atDomainMaxBorder( IBlock & block )
{
   return boost::python::make_tuple( block.getBlockStorage().atDomainXMaxBorder(block),
                                     block.getBlockStorage().atDomainYMaxBorder(block),
                                     block.getBlockStorage().atDomainZMaxBorder(block) );
}

IBlockID::IDType IBlock_getIntegerID( IBlock & block )
{
   return block.getId().getID();
}

bool IBlock_equals( IBlock & block1, IBlock & block2 )
{
   return block1.getId() == block2.getId();
}

std::string IBlock_str( IBlock & b ) {
   std::stringstream out;
   out <<  "Block at " << b.getAABB();
   return out.str();

}

void exportIBlock()
{
   register_exception_translator<NoSuchBlockData>( & NoSuchBlockData::translate );
   register_exception_translator<BlockDataNotConvertible>( & BlockDataNotConvertible::translate );

   class_<IBlock, boost::noncopyable> ("Block", no_init)
         .def         ( "__getitem__",       &IBlock_getData )
         .add_property( "atDomainMinBorder", &IBlock_atDomainMinBorder )
         .add_property( "atDomainMaxBorder", &IBlock_atDomainMaxBorder )
         .add_property( "items",             &IBlock_blockDataList)
         .add_property( "id",                &IBlock_getIntegerID)
         .def         ( "__hash__",          &IBlock_getIntegerID)
         .def         ( "__eq__",            &IBlock_equals)
         .def         ( "__repr__",          &IBlock_str )
         .add_property( "__iter__",          &IBlock_iter  )
         .add_property("aabb", make_function(&IBlock::getAABB, return_value_policy<copy_const_reference>()))
         ;

}

//======================================================================================================================
//
//  Logging & Abort
//
//======================================================================================================================


static void wlb_log_devel              ( const std::string & msg ) { WALBERLA_LOG_DEVEL          ( msg ); }
static void wlb_log_devel_on_root      ( const std::string & msg ) { WALBERLA_LOG_DEVEL_ON_ROOT  ( msg ); }

static void wlb_log_result             ( const std::string & msg ) { WALBERLA_LOG_RESULT         ( msg ); }
static void wlb_log_result_on_root     ( const std::string & msg ) { WALBERLA_LOG_RESULT_ON_ROOT ( msg ); }

static void wlb_log_warning            ( const std::string & msg ) { WALBERLA_LOG_WARNING         ( msg ); }
static void wlb_log_warning_on_root    ( const std::string & msg ) { WALBERLA_LOG_WARNING_ON_ROOT ( msg ); }

#ifdef WALBERLA_LOGLEVEL_INFO
static void wlb_log_info               ( const std::string & msg ) { WALBERLA_LOG_INFO            ( msg ); }
static void wlb_log_info_on_root       ( const std::string & msg ) { WALBERLA_LOG_INFO_ON_ROOT    ( msg ); }
#else
static void wlb_log_info               ( const std::string & ) {}
static void wlb_log_info_on_root       ( const std::string & ) {}
#endif

#ifdef WALBERLA_LOGLEVEL_PROGRESS
static void wlb_log_progress           ( const std::string & msg ) { WALBERLA_LOG_PROGRESS        ( msg ); }
static void wlb_log_progress_on_root   ( const std::string & msg ) { WALBERLA_LOG_PROGRESS_ON_ROOT( msg ); }
#else
static void wlb_log_progress           ( const std::string & ) {}
static void wlb_log_progress_on_root   ( const std::string & ) {}
#endif

#ifdef WALBERLA_LOGLEVEL_DETAIL
static void wlb_log_detail             ( const std::string & msg ) { WALBERLA_LOG_DETAIL          ( msg ); }
static void wlb_log_detail_on_root     ( const std::string & msg ) { WALBERLA_LOG_DETAIL_ON_ROOT  ( msg ); }
#else
static void wlb_log_detail             ( const std::string & ) {}
static void wlb_log_detail_on_root     ( const std::string & ) {}
#endif

static void wlb_abort                  ( const std::string & msg ) { WALBERLA_ABORT_NO_DEBUG_INFO ( msg ); }

void exportLogging()
{
   def ( "log_devel"         ,  wlb_log_devel           );
   def ( "log_devel_on_root" ,  wlb_log_devel_on_root   );
   def ( "log_result",          wlb_log_result          );
   def ( "log_result_on_root",  wlb_log_result_on_root  );
   def ( "log_warning",         wlb_log_warning         );
   def ( "log_warning_on_root", wlb_log_warning_on_root );
   def ( "log_info",            wlb_log_info            );
   def ( "log_info_on_root",    wlb_log_info_on_root    );
   def ( "log_progress",        wlb_log_progress        );
   def ( "log_progress_on_root",wlb_log_progress_on_root);
   def ( "log_detail",          wlb_log_detail          );
   def ( "log_detail_on_root",  wlb_log_detail_on_root  );

   def ( "abort", wlb_abort );
}




//======================================================================================================================
//
//  StructuredBlockStorage
//
//======================================================================================================================


object * blockDataCreationHelper( IBlock * block, StructuredBlockStorage * bs,  object callable ) //NOLINT
{
   object * res = new object( callable( ptr(block), ptr(bs) ) );
   return res;
}

uint_t StructuredBlockStorage_addBlockData( StructuredBlockStorage & s, const std::string & name, object functionPtr ) //NOLINT
{
   BlockDataID res = s.addStructuredBlockData(name)
               << StructuredBlockDataCreator<object>( std::bind( &blockDataCreationHelper, std::placeholders::_1, std::placeholders::_2, functionPtr ) );
   //TODO extend this for moving block data ( packing und unpacking with pickle )
   return res;
}

// Helper function for iteration over StructuredBlockStorage
// boost::python comes with iteration helpers but non of this worked:
//    .def("__iter__"   range(&StructuredBlockStorage::begin, &StructuredBlockStorage::end))
//    .def("__iter__",  range<return_value_policy<copy_non_const_reference> >( beginPtr, endPtr) )
boost::python::object StructuredBlockStorage_iter( boost::python::object structuredBlockStorage ) //NOLINT
{
   shared_ptr<StructuredBlockStorage> s = extract< shared_ptr<StructuredBlockStorage> > ( structuredBlockStorage );

   std::vector< const IBlock* > blocks;
   s->getBlocks( blocks );
   boost::python::list resultList;

   for( auto it = blocks.begin(); it != blocks.end(); ++it ) {
      boost::python::object theObject( ptr( *it ) );
      // Prevent blockstorage from being destroyed when references to blocks exist
      boost::python::objects::make_nurse_and_patient( theObject.ptr(), structuredBlockStorage.ptr() );
      resultList.append( theObject );
   }

   return resultList.attr("__iter__");
}


boost::python::object StructuredBlockStorage_getItem( boost::python::object structuredBlockStorage, uint_t i ) //NOLINT
{
   shared_ptr<StructuredBlockStorage> s = extract< shared_ptr<StructuredBlockStorage> > ( structuredBlockStorage );

   if ( i >= s->size() )
   {
      PyErr_SetString( PyExc_RuntimeError, "Index out of bounds");
      throw error_already_set();
   }

   std::vector< const IBlock* > blocks;
   s->getBlocks( blocks );

   boost::python::object theObject( ptr( blocks[i] ) );
   boost::python::objects::make_nurse_and_patient( theObject.ptr(), structuredBlockStorage.ptr() );
   return theObject;
}

boost::python::list StructuredBlockStorage_blocksOverlappedByAABB( StructuredBlockStorage & s, const AABB & aabb ) {
   std::vector< IBlock*> blocks;
   s.getBlocksOverlappedByAABB( blocks, aabb );

   boost::python::list resultList;
   for( auto it = blocks.begin(); it != blocks.end(); ++it )
      resultList.append( ptr( *it ) );
   return resultList;
}


boost::python::list StructuredBlockStorage_blocksContainedWithinAABB( StructuredBlockStorage & s, const AABB & aabb ) {
   std::vector< IBlock*> blocks;
   s.getBlocksContainedWithinAABB( blocks, aabb );

   boost::python::list resultList;
   for( auto it = blocks.begin(); it != blocks.end(); ++it )
      resultList.append( ptr( *it ) );
   return resultList;
}


object SbS_transformGlobalToLocal ( StructuredBlockStorage & s, IBlock & block, const object & global )
{
   if ( extract<CellInterval>( global ).check() )
   {
      CellInterval ret;
      s.transformGlobalToBlockLocalCellInterval( ret, block, extract<CellInterval>( global ) );
      return object( ret );
   }
   else if ( extract<Cell>( global ).check() )
   {
      Cell ret;
      s.transformGlobalToBlockLocalCell( ret, block, extract<Cell>( global ) );
      return object( ret );
   }

   PyErr_SetString(PyExc_RuntimeError, "Only CellIntervals and cells can be transformed" );
   throw error_already_set();
}


object SbS_transformLocalToGlobal ( StructuredBlockStorage & s, IBlock & block, const object & local )
{
   if ( extract<CellInterval>( local ).check() )
   {
      CellInterval ret;
      s.transformBlockLocalToGlobalCellInterval( ret, block, extract<CellInterval>( local ) );
      return object( ret );
   }
   else if ( extract<Cell>( local ).check() )
   {
      Cell ret;
      s.transformBlockLocalToGlobalCell( ret, block, extract<Cell>( local ) );
      return object( ret );
   }
   PyErr_SetString(PyExc_RuntimeError, "Only CellIntervals and cells can be transformed" );
   throw error_already_set();
}

void SbS_writeBlockData( StructuredBlockStorage & s,const std::string & blockDataId, const std::string & file )
{
   mpi::SendBuffer buffer;
   s.serializeBlockData( blockDataIDFromString(s, blockDataId), buffer);
   mpi::writeMPIIO(file, buffer);
}

void SbS_readBlockData( StructuredBlockStorage & s,const std::string & blockDataId, const std::string & file )
{
   mpi::RecvBuffer buffer;
   mpi::readMPIIO(file, buffer);

   s.deserializeBlockData( blockDataIDFromString(s, blockDataId), buffer );
   if ( ! buffer.isEmpty() ) {
      PyErr_SetString(PyExc_RuntimeError, "Reading failed - file does not contain matching data for this type." );
      throw error_already_set();
   }
}

CellInterval SbS_getBlockCellBB( StructuredBlockStorage & s, const IBlock * block )
{
   return s.getBlockCellBB( *block );
}


Vector3<real_t> SbS_mapToPeriodicDomain1 ( StructuredBlockStorage & s, real_t x, real_t y, real_t z )
{
   Vector3<real_t> res ( x,y,z );
   s.mapToPeriodicDomain( res );
   return res;
}
Vector3<real_t> SbS_mapToPeriodicDomain2 ( StructuredBlockStorage & s, Vector3<real_t> in )
{
   s.mapToPeriodicDomain( in );
   return in;
}
Cell SbS_mapToPeriodicDomain3 ( StructuredBlockStorage & s, Cell in, uint_t level = 0 )
{
   s.mapToPeriodicDomain( in, level );
   return in;
}

object SbS_getBlock1 ( StructuredBlockStorage & s, const real_t x , const real_t y , const real_t z ) {
   return object( ptr( s.getBlock( x,y,z ) ) );

}

object SbS_getBlock2 ( StructuredBlockStorage & s, const Vector3<real_t> & v ) {
   return object( ptr( s.getBlock( v ) ) );
}


tuple SbS_periodic(  StructuredBlockStorage & s )
{
   return make_tuple( s.isXPeriodic(), s.isYPeriodic(), s.isZPeriodic() );
}

bool SbS_atDomainXMinBorder( StructuredBlockStorage & s, const IBlock * b ) { return s.atDomainXMinBorder( *b ); }
bool SbS_atDomainXMaxBorder( StructuredBlockStorage & s, const IBlock * b ) { return s.atDomainXMaxBorder( *b ); }
bool SbS_atDomainYMinBorder( StructuredBlockStorage & s, const IBlock * b ) { return s.atDomainYMinBorder( *b ); }
bool SbS_atDomainYMaxBorder( StructuredBlockStorage & s, const IBlock * b ) { return s.atDomainYMaxBorder( *b ); }
bool SbS_atDomainZMinBorder( StructuredBlockStorage & s, const IBlock * b ) { return s.atDomainZMinBorder( *b ); }
bool SbS_atDomainZMaxBorder( StructuredBlockStorage & s, const IBlock * b ) { return s.atDomainZMaxBorder( *b ); }

void exportStructuredBlockStorage()
{
   bool ( StructuredBlockStorage::*p_blockExists1         ) ( const Vector3< real_t > & ) const = &StructuredBlockStorage::blockExists;
   bool ( StructuredBlockStorage::*p_blockExistsLocally1  ) ( const Vector3< real_t > & ) const = &StructuredBlockStorage::blockExistsLocally;
   bool ( StructuredBlockStorage::*p_blockExistsRemotely1 ) ( const Vector3< real_t > & ) const = &StructuredBlockStorage::blockExistsRemotely;

   bool ( StructuredBlockStorage::*p_blockExists2         ) ( const real_t, const real_t, const real_t ) const = &StructuredBlockStorage::blockExists;
   bool ( StructuredBlockStorage::*p_blockExistsLocally2  ) ( const real_t, const real_t, const real_t ) const = &StructuredBlockStorage::blockExistsLocally;
   bool ( StructuredBlockStorage::*p_blockExistsRemotely2 ) ( const real_t, const real_t, const real_t ) const = &StructuredBlockStorage::blockExistsRemotely;

   class_<StructuredBlockStorage, shared_ptr<StructuredBlockStorage>, boost::noncopyable>("StructuredBlockStorage", no_init )
       .def( "getNumberOfLevels",                       &StructuredBlockStorage::getNumberOfLevels )
       .def( "getDomain",                               &StructuredBlockStorage::getDomain, return_internal_reference<1>() )
       .def( "mapToPeriodicDomain",                     &SbS_mapToPeriodicDomain1                                 )
       .def( "mapToPeriodicDomain",                     &SbS_mapToPeriodicDomain2                                 )
       .def( "mapToPeriodicDomain",                     &SbS_mapToPeriodicDomain3, (arg("level") = 0 )            )
       .def( "addBlockData",                            &StructuredBlockStorage_addBlockData                      )
       .def( "__getitem__",                             &StructuredBlockStorage_getItem                           )
       .def( "__len__",                                 &StructuredBlockStorage::size                             )
       .def( "getBlock",                                SbS_getBlock1                                             )
       .def( "getBlock",                                SbS_getBlock2                                             )
       .def( "containsGlobalBlockInformation",          &StructuredBlockStorage::containsGlobalBlockInformation   )
       .def( "blocksOverlappedByAABB" ,                 &StructuredBlockStorage_blocksOverlappedByAABB            )
       .def( "blocksContainedWithinAABB",               &StructuredBlockStorage_blocksContainedWithinAABB         )
       .def( "blockExists",                             p_blockExists1                                            )
       .def( "blockExists",                             p_blockExists2                                            )
       .def( "blockExistsLocally",                      p_blockExistsLocally1                                     )
       .def( "blockExistsLocally",                      p_blockExistsLocally2                                     )
       .def( "blockExistsRemotely",                     p_blockExistsRemotely1                                    )
       .def( "blockExistsRemotely",                     p_blockExistsRemotely2                                    )
       .def( "atDomainXMinBorder",                      &SbS_atDomainXMinBorder                                   )
       .def( "atDomainXMaxBorder",                      &SbS_atDomainXMaxBorder                                   )
       .def( "atDomainYMinBorder",                      &SbS_atDomainYMinBorder                                   )
       .def( "atDomainYMaxBorder",                      &SbS_atDomainYMaxBorder                                   )
       .def( "atDomainZMinBorder",                      &SbS_atDomainZMinBorder                                   )
       .def( "atDomainZMaxBorder",                      &SbS_atDomainZMaxBorder                                   )
       .def( "dx",                                      &StructuredBlockStorage::dx, ( args("level")=0 )          )
       .def( "dy",                                      &StructuredBlockStorage::dy, ( args("level")=0 )          )
       .def( "dz",                                      &StructuredBlockStorage::dz, ( args("level")=0 )          )
       .def( "getDomainCellBB",                         &StructuredBlockStorage::getDomainCellBB, return_value_policy<copy_const_reference>(),  ( args( "level") = 0 )  )
       .def( "getBlockCellBB",                          &SbS_getBlockCellBB  )
       .def( "transformGlobalToLocal",                  &SbS_transformGlobalToLocal )
       .def( "transformLocalToGlobal",                  &SbS_transformLocalToGlobal )
       .def( "writeBlockData",                          &SbS_writeBlockData )
       .def( "readBlockData",                           &SbS_readBlockData )
       .add_property("__iter__",                        &StructuredBlockStorage_iter                            )
       .add_property( "containsGlobalBlockInformation", &StructuredBlockStorage::containsGlobalBlockInformation )
       .add_property( "periodic",                       &SbS_periodic )
       ;

#if BOOST_VERSION < 106300
   exportSharedPtr<StructuredBlockStorage>();
#endif
}

//======================================================================================================================
//
//  Communication
//
//======================================================================================================================

void exportCommunication()
{
   using communication::UniformPackInfo;
   class_< UniformPackInfo, shared_ptr<UniformPackInfo>, boost::noncopyable> //NOLINT
      ( "UniformPackInfo", no_init );

   using communication::UniformMPIDatatypeInfo;
   class_< UniformMPIDatatypeInfo, shared_ptr<UniformMPIDatatypeInfo>, boost::noncopyable>
      ( "UniformMPIDatatypeInfo", no_init );

#if BOOST_VERSION < 106300
   exportSharedPtr<UniformPackInfo>();
   exportSharedPtr<UniformMPIDatatypeInfo>();
#endif
}

//======================================================================================================================
//
//  Stencil Directions
//
//======================================================================================================================

void exportStencilDirections()
{
   ModuleScope build_info( "stencil");

   enum_<stencil::Direction>("Direction")
       .value("C"  , stencil::C  )
       .value("N"  , stencil::N  )
       .value("S"  , stencil::S  )
       .value("W"  , stencil::W  )
       .value("E"  , stencil::E  )
       .value("T"  , stencil::T  )
       .value("B"  , stencil::B  )
       .value("NW" , stencil::NW )
       .value("NE" , stencil::NE )
       .value("SW" , stencil::SW )
       .value("SE" , stencil::SE )
       .value("TN" , stencil::TN )
       .value("TS" , stencil::TS )
       .value("TW" , stencil::TW )
       .value("TE" , stencil::TE )
       .value("BN" , stencil::BN )
       .value("BS" , stencil::BS )
       .value("BW" , stencil::BW )
       .value("BE" , stencil::BE )
       .value("TNE", stencil::TNE)
       .value("TNW", stencil::TNW)
       .value("TSE", stencil::TSE)
       .value("TSW", stencil::TSW)
       .value("BNE", stencil::BNE)
       .value("BNW", stencil::BNW)
       .value("BSE", stencil::BSE)
       .value("BSW", stencil::BSW)
       .export_values()
       ;
   boost::python::list cx;

   boost::python::list cy;

   boost::python::list cz;

   boost::python::list dirStrings;
   for( uint_t i=0; i < stencil::NR_OF_DIRECTIONS; ++i  )
   {
      cx.append( stencil::cx[i] );
      cy.append( stencil::cy[i] );
      cz.append( stencil::cz[i] );
      dirStrings.append( stencil::dirToString[i] );
   }
   boost::python::list c;
   c.append( cx );
   c.append( cy );
   c.append( cz );

   using boost::python::scope;
   scope().attr("cx") = cx;
   scope().attr("cy") = cy;
   scope().attr("cz") = cz;
   scope().attr("c") = c;
   scope().attr("dirStrings") = dirStrings;
}


//======================================================================================================================
//
//  Build Info
//
//======================================================================================================================


void exportBuildInfo()
{
   ModuleScope build_info( "build_info");
   using boost::python::scope;
   scope().attr("version")         = WALBERLA_GIT_SHA1;
   scope().attr("type" )           = WALBERLA_BUILD_TYPE;
   scope().attr("compiler_flags" ) = WALBERLA_COMPILER_FLAGS;
   scope().attr("build_machine" )  = WALBERLA_BUILD_MACHINE;
   scope().attr("source_dir")      = WALBERLA_SOURCE_DIR;
   scope().attr("build_dir")       = WALBERLA_BUILD_DIR;
}




void exportBasicWalberlaDatastructures()
{

   NumpyIntConversion<uint8_t>();
   NumpyIntConversion<int32_t>();
   NumpyIntConversion<int64_t>();
   NumpyIntConversion<uint32_t>();
   NumpyIntConversion<uint64_t>();
   NumpyIntConversion<size_t>();
   NumpyIntConversion<bool>();
   NumpyFloatConversion<float>();
   NumpyFloatConversion<double>();
   NumpyFloatConversion<long double>();


   exportMPI();

   exportBuildInfo();
   exportVector3();
   exportCell();
   exportCellInterval();
   exportAABB();

   exportTiming();

   exportIBlock();
   exportStructuredBlockStorage();
   exportCommunication();

   exportLogging();
   exportStencilDirections();

   // Add empty callbacks module
   object callbackModule( handle<>( borrowed(PyImport_AddModule("walberla_cpp.callbacks"))));
   scope().attr("callbacks") = callbackModule;

}

} // namespace python_coupling
} // namespace walberla

#endif

