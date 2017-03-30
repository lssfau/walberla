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
//! \file OutputColor.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "OutputColor.h"


namespace walberla {
namespace blockforest {



void writeColor( std::ofstream& outfile, const uint_t index, const uint_t selector ) {

   if( selector == PARAVIEW ) {

      switch( index ) {

         case 0:
            outfile << "0 1 0 1";
            break;
         case 1:
            outfile << "0 0 1 1";
            break;
         case 2:
            outfile << "1 0 1 1";
            break;
         case 3:
            outfile << "1 0 0 1";
            break;
         default:
            outfile << "1 1 0 1";
      }
   }
   else if( selector == GNUPLOT ) {

      switch( index ) {

         case 0:
            outfile << "#00FF00";
            break;
         case 1:
            outfile << "#0000FF";
            break;
         case 2:
            outfile << "#FF00FF";
            break;
         case 3:
            outfile << "#FF0000";
            break;
         default:
            outfile << "#000000";
      }
   }
}



} // namespace blockforest
} // namespace walberla
