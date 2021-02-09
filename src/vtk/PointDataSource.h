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
//! \file PointDataSource.h
//! \ingroup vtk
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "Base64Writer.h"
#include "UtilityFunctions.h"
#include "core/DataTypes.h"
#include "core/math/Vector3.h"

#include <ostream>
#include <string>
#include <vector>


namespace walberla {
namespace vtk {



class PointDataSource {

public:

   struct Attributes
   {
      Attributes( const std::string & t, const std::string & n, const uint_t c ) :
         type(t), name(n), components(c) {}
      std::string type;
      std::string name;
      uint_t      components;
   };

   virtual ~PointDataSource() = default;

   virtual std::vector< Attributes > getAttributes() const = 0;

   virtual void configure() = 0;

   virtual std::vector< Vector3< real_t > > getPoints() = 0;

   /// For the documentation of this function, please refer to the documentation/general description of this class.
   virtual void push( std::ostream& os,  const uint_t data, const uint_t point, const uint_t component ) = 0;
   /// For the documentation of this function, please refer to the documentation/general description of this class.
   virtual void push( Base64Writer& b64, const uint_t data, const uint_t point, const uint_t component ) = 0;

}; // class PointDataSource



} // namespace vtk
} // namespace walberla
