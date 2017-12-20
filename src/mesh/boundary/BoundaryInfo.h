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
//! \file BoundaryInfo.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "boundary/BoundaryUID.h"
#include "boundary/Boundary.h"

namespace walberla {
namespace mesh {

class BoundaryInfo
{
public:
   typedef boundary::BoundaryUID UID;
   typedef shared_ptr< boundary::BoundaryConfiguration > ConfigPtr;

   BoundaryInfo() : config_( boundary::BoundaryConfiguration::nullPtr() ) { }
   BoundaryInfo( const UID & uid ) : uid_( uid ), config_( boundary::BoundaryConfiguration::nullPtr() ) { }
   BoundaryInfo( const UID & uid, const ConfigPtr & config ) : uid_( uid ), config_( config ) { }

   const UID &       getUid()    const { return uid_;    }
   const ConfigPtr & getConfig() const { return config_; }

   void setUid   ( const UID &       uid    ) { uid_    = uid;    }
   void setConfig( const ConfigPtr & config ) { config_ = config; }

private:
   UID       uid_;
   ConfigPtr config_;
};

} // namespace mesh
} // namespace walberla