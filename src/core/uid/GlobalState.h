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
//! \file GlobalState.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "SUID.h"

#include "core/debug/Debug.h"
#include "core/singleton/Singleton.h"


namespace walberla {
namespace uid {



//**********************************************************************************************************************
/*!
*   \brief The global state of the simulation
*
*   Can be configured using the member function "configure". Must only be configured once while the simulation is
*   initialized. If not configured, the global simulation state defaults to "Set<SUID>::emptySet()".
*
*   The global simulation state can be retrieved via calling the function "globalState()".
*/
//**********************************************************************************************************************

class GlobalState : public singleton::Singleton< GlobalState > {

public:

   WALBERLA_BEFRIEND_SINGLETON;

   /// Must only be called once!
   /// (fails in debug mode if called for a second time - always succeeds in release mode, however only the first call has an effect)
   void configure( const Set<SUID>& state )   { WALBERLA_ASSERT( !configured_ ); if( !configured_ ) { state_ = state; configured_ = true; } }

   void reconfigure( const Set<SUID> & state) { WALBERLA_ASSERT( configured_);  if ( configured_ ) { state_ = state; } }

   const Set<SUID>& get() const { return state_; }

private:

   bool      configured_;
   Set<SUID> state_;

   GlobalState() : configured_( false ) {}

};



/// convenience function for accessing the global simulation state
inline const Set<SUID>& globalState() { return GlobalState::instance()->get(); }



} // namespace uid
} // namespace walberla
