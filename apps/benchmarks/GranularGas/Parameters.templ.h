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
//! \file Parameters.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <core/config/Config.h>
#include <mesa_pd/data/DataTypes.h>

#include <string>

namespace walberla {
namespace mesa_pd {

struct Parameters
{
   {%- for param in parameters %}
   {{param.type}} {{param.name}} = {{param.defValue}};
   {%- endfor %}
};

void loadFromConfig(Parameters& params,
                    const Config::BlockHandle& cfg);

void saveToSQL(const Parameters& params,
               std::map< std::string, walberla::int64_t >& integerProperties,
               std::map< std::string, double >&            realProperties,
               std::map< std::string, std::string >&       stringProperties );

} //namespace mesa_pd
} //namespace walberla
