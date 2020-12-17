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
//! \file
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <core/DataTypes.h>

#include <typeinfo>
#include <type_traits>

namespace walberla {
namespace mesa_pd {
namespace data {

/**
 * @brief Add this as a base class to identify your class as a accessor.
 *
 * Accessors passed via templated arguments might be checked like
 * \code
 * static_assert(std::is_base_of<IAccessor, X>::value, typeid(X).name() + " is not a accessor");
 * \endcode
 */
class IAccessor
{
public:
   virtual ~IAccessor() = default;
};

} //namespace data
} //namespace mesa_pd
} //namespace walberla
