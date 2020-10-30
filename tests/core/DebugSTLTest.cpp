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
//! \file GridGeneratorTest.cpp
//! \author Dominik Thoennes <dominik.thoennes@fau.de>
//
//======================================================================================================================

#include <cstdlib>
#include <vector>

#include "waLBerlaDefinitions.h"

// this test is expected to fail
int main(int /*argc*/, char** /*argv*/)
{
   int ret = EXIT_FAILURE;
#ifdef WALBERLA_STL_BOUNDS_CHECKS
   std::vector< int > a(100);
   // this throws an exception if the debug STL is used
   // otherwise main will return 0 and the test fails since it is expected to fail
   if (a[200] != 1337) { ret = EXIT_SUCCESS; }
#endif
   return ret;
}