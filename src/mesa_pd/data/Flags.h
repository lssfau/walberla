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
//! \file Flags.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>

#include <ostream>

namespace walberla {
namespace mesa_pd {
namespace data {

namespace particle_flags {

class FlagT
{
public:
   using value_type = uint8_t;

   value_type& getRawDataRef() {return flags_;}
   const value_type& getRawDataRef() const {return flags_;}
private:
   value_type flags_ = 0;
};

/// default value after initialization
static const FlagT::value_type NO_FLAGS          = 0;
/// the particle extends into at least one direction till infinity
static const FlagT::value_type INFINITE          = (1 << 1);
/// spatial integrators will skip this particle, however, forces and torques are reset
static const FlagT::value_type FIXED             = (1 << 2);
/// all communication functions are skipped for this particle
static const FlagT::value_type NON_COMMUNICATING = (1 << 3);
/// this is a ghost particle
static const FlagT::value_type GHOST             = (1 << 4);
/// this particle exists on all processes with the same uid
static const FlagT::value_type GLOBAL            = (1 << 5);

inline void set(FlagT& flags, const FlagT::value_type flag)
{
   flags.getRawDataRef() |= flag;
}

inline void unset(FlagT& flags, const FlagT::value_type flag)
{
   flags.getRawDataRef() &= static_cast<FlagT::value_type>(~flag);
}

inline bool isSet(const FlagT& flags, const FlagT::value_type flag)
{
   return (flags.getRawDataRef() & flag) != 0;
}

inline std::ostream& operator<<( std::ostream& os, const FlagT& flags )
{
   os << flags.getRawDataRef();
   return os;
}

} //namespace particle_flags
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
mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf, const walberla::mesa_pd::data::particle_flags::FlagT& flags )
{
   buf.addDebugMarker( "fl" );
   buf << flags.getRawDataRef();
   return buf;
}

template< typename T>    // Element type  of RecvBuffer
mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf, walberla::mesa_pd::data::particle_flags::FlagT& flags )
{
   buf.readDebugMarker( "fl" );
   buf >> flags.getRawDataRef();
   return buf;
}

template<> // value type
struct BufferSizeTrait< walberla::mesa_pd::data::particle_flags::FlagT > {
   static const bool constantSize = true;
   static const uint_t size = sizeof(walberla::mesa_pd::data::particle_flags::FlagT::value_type) + mpi::BUFFER_DEBUG_OVERHEAD;
};

} // mpi
} // walberla
