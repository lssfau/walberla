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
//! \file UniqueID.cpp
//! \ingroup core
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/ConcatIterator.h"

#include <complex>
#include <list>
#include <vector>

using namespace walberla;

int main( int /*argc*/, char** /*argv*/ )
{
   walberla::debug::enterTestMode();

   std::list<std::complex<double> > a(9, std::complex<double>(1., 2.));
   std::list<std::complex<double> > b(5, std::complex<double>(3., 4.));

   for (auto it = ConcatIterator< std::list<std::complex<double> >::iterator >(a.begin(), a.end(), b.begin(), b.end());
        it != ConcatIterator< std::list<std::complex<double> >::iterator >();
        ++it)
   {
      (*it) = std::complex<double>(99., 98.);
   }

   for (auto it = a.begin(); it != a.end(); ++it)
   {
      WALBERLA_CHECK_FLOAT_EQUAL(it.operator->()->real(), 99.0);
      WALBERLA_CHECK_FLOAT_EQUAL((*it).imag(), 98.0);
   }

   for (auto it = b.begin(); it != b.end(); ++it)
   {
      WALBERLA_CHECK_FLOAT_EQUAL(it->real(), 99.0);
      WALBERLA_CHECK_FLOAT_EQUAL((*it).imag(), 98.0);
   }

   for (auto it = ConcatIterator< std::list<std::complex<double> >::iterator >(a.begin(), a.end(), a.end(), a.end());
        it != ConcatIterator< std::list<std::complex<double> >::iterator >();
        ++it)
   {
      WALBERLA_CHECK_FLOAT_EQUAL(it.operator->()->real(), 99.0);
      WALBERLA_CHECK_FLOAT_EQUAL((*it).imag(), 98.0);
   }

   for (auto it = ConcatIterator< std::list<std::complex<double> >::iterator >(b.begin(), b.end(), b.end(), b.end());
        it != ConcatIterator< std::list<std::complex<double> >::iterator >();
        ++it)
   {
      WALBERLA_CHECK_FLOAT_EQUAL(it->real(), 99.0);
      WALBERLA_CHECK_FLOAT_EQUAL((*it).imag(), 98.0);
   }

   return EXIT_SUCCESS;
}
