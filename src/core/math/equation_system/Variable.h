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
//! \file Variable.h
//! \ingroup core
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#pragma once

#include "FwdVariable.h"


namespace walberla {
namespace math {

   class Var
   {
   private:
      const std::string name_;

      bool   valid_;
      double value_;

   public:
      Var ( const std::string& name );
   private:
      Var& operator=( const Var& ){ return *this; }

   public:
      bool operator==( const Var& var) const;

   public:
            bool         valid() const { return valid_; }
            double       getValue() const { return value_; }
      const std::string& getName()  const { return name_;  }

   public:
      void setValue( const double value );

   public:
      friend std::ostream& operator<<( std::ostream& os, const Var & var );
   };
   // end class Var

} // namespace math
} // namespace walberla
