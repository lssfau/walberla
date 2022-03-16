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
//! \file UnrollTest.cpp
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "lbm/IntelCompilerOptimization.h"
#include "core/DataTypes.h"

#include <iostream>
#include <sstream>


using namespace std;
using namespace walberla;

int main( int , char**argv )
{
   stringstream ss (argv[1]);

   cell_idx_t xSize = 20;
   ss >> xSize;
   cout << "XSize = " << xSize << endl;
   X_LOOP
   (
      cout << x << endl;
   )

   return 0;
}
