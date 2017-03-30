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
//! \file WeakPtrWrapper.h
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"
#include "core/DataTypes.h"


/*
   What is this good for?

      Memory Management of basic waLBerla objects: Use BlockStorage and CommunicationSchemes as example:
          - a block storage is held as a shared_ptr
          - create a communication scheme: internally the communication scheme
            only holds a weak pointer to the block storage i.e. communication does not "own" the block storage.
            If the communication would own the block storage and would internally hold a shared_ptr instead of a
            weak_ptr, cycles could be created e.g. when a communication object is passed to the block storage,
            which happens e.g. when registering communication functors at the block storage itself
            -> weak_ptr is necessary to prevent cycles ( another way would be not to store a pointer to
            block storage at all which is hardly practicable  )
         - trying to communicate after a block storage was deleted leads to an error, which is reasonable

      The problem occurs when exporting this to Python:
         - Due to a known and probably never fixed bug, one cannot created weak pointers from shared pointers
           received from Python functions.  Details here:
              - http://stackoverflow.com/questions/8233252/boostpython-and-weak-ptr-stuff-disappearing
              - https://svn.boost.org/trac/boost/ticket/3673
         - What works are plain old C pointers -> when compiling with Python plain C pointers should be used instead
           of the checked weak pointers -> this is the reason for this wrapper class
         - Special attention is necessary when exporting this construct to Python since a block forest may only be
           deleted when no communication object points to it. This is ensured in this case via
           UniformBufferedSchemeWrapper
 */







namespace walberla {


#ifndef WALBERLA_BUILD_WITH_PYTHON



template<typename T>
class weak_ptr_wrapper : public weak_ptr<T>
{
public:
   weak_ptr_wrapper() {}
   weak_ptr_wrapper( const weak_ptr  <T> & r ) : weak_ptr<T>(r) {}
   weak_ptr_wrapper( const shared_ptr<T> & r ) : weak_ptr<T>(r) {}
};


# else

   // Due to a bug in boost::python weak_ptr cannnot be used:
   // http://stackoverflow.com/questions/8233252/boostpython-and-weak-ptr-stuff-disappearing
   template<typename T>
   class weak_ptr_wrapper
   {
   public:
      weak_ptr_wrapper( const shared_ptr<T> & sp )
         : rawPtr_ ( sp.get() )
      {}

      T * lock() { return rawPtr_; }

   protected:
      T * rawPtr_;
   };


#endif // WALBERLA_BUILD_WITH_PYTHON

} // namespace walberla


