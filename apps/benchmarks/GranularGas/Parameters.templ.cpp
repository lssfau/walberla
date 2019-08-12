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
//! \file Parameters.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#include "Parameters.h"

#include <core/logging/Logging.h>

namespace walberla {
namespace mesa_pd {

void loadFromConfig(Parameters& params, const Config::BlockHandle& cfg)
{
   {%- for param in parameters %}
   params.{{param.name}} = cfg.getParameter<{{param.type}}>("{{param.name}}", {{param.defValue}} );
   WALBERLA_LOG_INFO_ON_ROOT("{{param.name}}: " << params.{{param.name}});
   {% endfor %}
}

void saveToSQL(const Parameters& params,
               std::map< std::string, walberla::int64_t >& integerProperties,
               std::map< std::string, double >&            realProperties,
               std::map< std::string, std::string >&       stringProperties )
{
   {%- for param in parameters %}
   {%- if param.type=="int64_t" %}
   integerProperties["{{param.name}}"] = params.{{param.name}};
   {%- endif %}
   {%- if param.type=="real_t" %}
   realProperties["{{param.name}}"] = double_c(params.{{param.name}});
   {%- endif %}
   {%- if param.type=="std::string" %}
   stringProperties["{{param.name}}"] = params.{{param.name}};
   {%- endif %}
   {% endfor %}
}

} //namespace mesa_pd
} //namespace walberla
