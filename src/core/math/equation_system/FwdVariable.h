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
//! \file FwdVariable.h
//! \ingroup core
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#pragma once

#include <boost/shared_ptr.hpp>
#include <map>
#include <string>


namespace walberla {
namespace math {

   extern double NAN_VAL;

   class Var;

   typedef boost::shared_ptr<Var> VarPtr;

   typedef std::map<std::string, VarPtr>                 VarMap;
   typedef std::map<std::string, VarPtr>::const_iterator VarMapIt;

} // namespace math
} // namespace walberla
