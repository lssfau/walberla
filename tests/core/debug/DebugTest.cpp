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
//! \file DebugTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Implements Tests for the file Debug.h
//
//======================================================================================================================

#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"

#include <limits>


int main()
{
   walberla::debug::enterTestMode();

   WALBERLA_ASSERT(true);
   WALBERLA_ASSERT_NULLPTR(reinterpret_cast<void*>(0));
   WALBERLA_ASSERT_NOT_NULLPTR(reinterpret_cast<void*>(1));
   WALBERLA_ASSERT_EQUAL(1, 1);
   WALBERLA_ASSERT_FLOAT_EQUAL(1.0, 1.0 + std::numeric_limits<double>::epsilon());
   WALBERLA_ASSERT_UNEQUAL(1,2);
   WALBERLA_ASSERT_LESS(1,2);
   WALBERLA_ASSERT_GREATER(2,1);
   WALBERLA_ASSERT_LESS_EQUAL(1,1);
   WALBERLA_ASSERT_LESS_EQUAL(1,2);
   WALBERLA_ASSERT_GREATER_EQUAL(1,1);
   WALBERLA_ASSERT_GREATER_EQUAL(2,1);

   WALBERLA_ASSERT(true, "answer: " << 42);
   WALBERLA_ASSERT_NULLPTR(reinterpret_cast<void*>(0), "answer: " << 42);
   WALBERLA_ASSERT_NOT_NULLPTR(reinterpret_cast<void*>(1), "answer: " << 42);
   WALBERLA_ASSERT_EQUAL(1, 1, "answer: " << 42);
   WALBERLA_ASSERT_FLOAT_EQUAL(1.0, 1.0 + std::numeric_limits<double>::epsilon(), "answer: " << 42);
   WALBERLA_ASSERT_UNEQUAL(1,2, "answer: " << 42);
   WALBERLA_ASSERT_LESS(1,2, "answer: " << 42);
   WALBERLA_ASSERT_GREATER(2,1, "answer: " << 42);
   WALBERLA_ASSERT_LESS_EQUAL(1,1, "answer: " << 42);
   WALBERLA_ASSERT_LESS_EQUAL(1,2, "answer: " << 42);
   WALBERLA_ASSERT_GREATER_EQUAL(1,1, "answer: " << 42);
   WALBERLA_ASSERT_GREATER_EQUAL(2,1, "answer: " << 42);

   WALBERLA_DEBUG_SECTION()
   {
      WALBERLA_ASSERT(true);
   }
   else
   {
      WALBERLA_ASSERT(false); //should be OK, asserts don't get evaluated in debug mode
   }

   WALBERLA_ASSERT_SECTION(true)
   {
      WALBERLA_ASSERT(false); // should never be evaluated, since condition is true
   }

   return 0;
}




