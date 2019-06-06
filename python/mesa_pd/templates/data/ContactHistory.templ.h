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
//! \file ContactHistory.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
{%- for include in includes %}
#include <{{include}}>
{%- endfor %}
#include <mesa_pd/data/STLOverloads.h>

#include <core/Abort.h>
#include <core/debug/Debug.h>
#include <core/math/AABB.h>
#include <core/mpi/Datatype.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>
#include <core/STLIO.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace data {

class ContactHistory
{
public:
   {%- for prop in properties %}
   const {{prop.type}}& get{{prop.name | capFirst}}() const {return {{prop.name}}_;}
   {{prop.type}}& get{{prop.name | capFirst}}Ref() {return {{prop.name}}_;}
   void set{{prop.name | capFirst}}(const {{prop.type}}& v) { {{prop.name}}_ = v;}
   {% endfor %}
private:
   {%- for prop in properties %}
   {{prop.type}} {{prop.name}}_ {};
   {%- endfor %}
};

inline
std::ostream& operator<<( std::ostream& os, const ContactHistory& ch )
{
   os << "==========  Contact History  ==========" << "\n" <<
   {%- for prop in properties %}
         "{{'%-20s'|format(prop.name)}}: " << ch.get{{prop.name | capFirst}}() << "\n" <<
   {%- endfor %}
         "================================" << std::endl;
   return os;
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla

//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

template< typename T,    // Element type of SendBuffer
          typename G>    // Growth policy of SendBuffer
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const mesa_pd::data::ContactHistory& obj )
{
   buf.addDebugMarker( "ch" );
   {%- for prop in properties %}
   buf << obj.get{{prop.name | capFirst}}();
   {%- endfor %}
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, mesa_pd::data::ContactHistory& objparam )
{
   buf.readDebugMarker( "ch" );
   {%- for prop in properties %}
   buf >> objparam.get{{prop.name | capFirst}}Ref();
   {%- endfor %}
   return buf;
}

} // mpi
} // walberla
