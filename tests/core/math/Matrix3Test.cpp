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
//! \file Matrix3Test.cpp
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/math/Matrix3.h"

#include <iostream>


using namespace walberla;
using walberla::uint8_t;


void rotationTest()
{
   Matrix3<real_t> rotationMatrix (0.0);
   rotationMatrix(0,0) = 1.0;
   rotationMatrix(1,1) = 1.0;
   rotationMatrix(2,2) = 1.0;

   Matrix3<real_t> diagonalMatrix ( 0.0 );
   diagonalMatrix(0,0) = 2.0;
   diagonalMatrix(1,1) = 4.0;
   diagonalMatrix(2,2) = 6.0;

   Matrix3<real_t> result = rotationMatrix.rotate( diagonalMatrix );

   std::cout << diagonalMatrix  << std::endl;
   std::cout << result << std::endl;

   WALBERLA_CHECK_FLOAT_EQUAL( result(0,0), 2.0 );
   WALBERLA_CHECK_FLOAT_EQUAL( result(1,1), 4.0 );
   WALBERLA_CHECK_FLOAT_EQUAL( result(2,2), 6.0 );


   for( uint_t i = 0; i < 3; ++i )
      for( uint_t j = 0; j < 3; ++j )
         if ( i != j)
            WALBERLA_CHECK_FLOAT_EQUAL( result(i,j), 0.0 );

   //also checking WALBERLA_CHECK_FLOAT_EQUAL for matrices
   Matrix3<real_t> cmp(2,0,0,0,4,0,0,0,6);
   WALBERLA_CHECK_FLOAT_EQUAL( result, cmp );
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( result, cmp, real_t(1e-5) );
}

void RARTTest()
{
   Matrix3<real_t> A ( 1,2,3,4,5,6,7,8,9 );
   Matrix3<real_t> R ( 2,3,4,5,6,7,8,9,1 );
   WALBERLA_CHECK_FLOAT_EQUAL( math::transformMatrixRART(R,A), R * A * R.getTranspose() );
}

int main()
{

   Matrix3<real_t> m1 ( 1.0 );
   Matrix3<real_t> m2 ( 2.0 );

   // the following line gives a compile error when the operator*(Other, Matrix3) is commented in
   m1 * m2;

   rotationTest();
   RARTTest();

   return 0;
}
