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
//! \file Variable.cpp
//! \ingroup core
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#include "Variable.h"

#include <cmath>
#include <sstream>


namespace walberla {
namespace math {

   Var::Var ( const std::string& name ) :
      name_ (name),
      valid_ (false),
      value_ (FP_NAN)
   {}

   void Var::setValue( const double value ){
      value_ = value;
      valid_ = !std::isnan( value );
   }

   bool Var::operator==( const Var& var) const {
      return name_ == var.name_;
   }

   std::ostream& operator<<( std::ostream& os, const Var & var ){
      return os << var.name_ << " = " << var.value_;
   }

} // namespace math
} // namespace walberla
