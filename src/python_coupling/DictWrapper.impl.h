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
//! \file DictWrapper.impl.h
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include <functional>

namespace walberla {
namespace python_coupling {


//===================================================================================================================
//
//  Setting
//
//===================================================================================================================


template<typename T>
void DictWrapper::exposePtr( const std::string & name, T * var ) {
   this->d_[name] =  boost::python::ptr( var );
}

template<typename T>
void DictWrapper::exposePtr( const std::string & name, const shared_ptr<T> & var ) {
   this->d_[name] =  boost::python::ptr( var.get() );
}

template<typename T>
void DictWrapper::exposeValue( const std::string & name, const T & var ) {
   this->d_[name] = var;
}



//===================================================================================================================
//
//  Getting default
//
//===================================================================================================================


template<typename T>
T DictWrapper::get( const std::string & name ) {
   return boost::python::extract<T>( d_[name] );
}

template<typename T>
bool DictWrapper::has( const std::string & name )
{
   if(! d_.has_key(name) )
      return false;

   return boost::python::extract<T>( d_[name]).check();
}

template<typename T>
bool DictWrapper::checkedGet( const std::string & name, T output )
{
   if ( ! has<T>(name) )
      return false;

   output = get<T>(name);
   return true;
}





template<>
inline DictWrapper DictWrapper::get( const std::string & name ) {
   auto dictCopy =  boost::python::extract< boost::python::dict >( d_[name] );
   DictWrapper result;
   result.dict() = dictCopy;
   return result;
}

template<>
inline bool DictWrapper::has<DictWrapper >( const std::string & name )
{
   if(! d_.has_key(name) )
      return false;

   return boost::python::extract< boost::python::dict >( d_[name]).check();
}


//===================================================================================================================
//
//  Getting std::functions
//
//===================================================================================================================

// void()
inline void runPythonObject( boost::python::object obj ) {
   obj();
}
template<>
inline std::function<void()> DictWrapper::get( const std::string & name ) {
   boost::python::object obj ( d_[name] );
   return std::bind( &runPythonObject, obj );
}

template<>
inline bool DictWrapper::has<std::function<void()> >( const std::string & name )
{
   if(! d_.has_key(name) )
      return false;

   return PyCallable_Check( boost::python::object(d_[name]).ptr() ) != 0;
}


} // namespace python_coupling
} // namespace walberla




