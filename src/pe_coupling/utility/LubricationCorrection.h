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
//! \file LubricationCorrection.h
//! \ingroup pe_coupling
//! \author Kristina Pickl <kristina.pickl@fau.de>
//! \author Dominik Bartuschat
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once


#include "domain_decomposition/StructuredBlockStorage.h"

#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Sphere.h"

namespace walberla {
namespace pe_coupling {

class LubricationCorrection
{
public:

   // constructor
   LubricationCorrection ( const shared_ptr<StructuredBlockStorage> & blockStorage, const shared_ptr<pe::BodyStorage> & globalBodyStorage,
                           const BlockDataID & bodyStorageID, real_t dynamicViscosity,
                           real_t cutOffDistance = real_t(2) / real_t(3), real_t minimalGapSize = real_t(1e-5) )
      : blockStorage_ ( blockStorage )
      , globalBodyStorage_( globalBodyStorage )
      , bodyStorageID_( bodyStorageID )
      , dynamicViscosity_( dynamicViscosity )
      , cutOffDistance_( cutOffDistance )
      , minimalGapSize_( minimalGapSize )
   { }

   void operator()();

private:

   // helper functions
   void treatLubricationSphrSphr ( const pe::SphereID sphereI, const pe::SphereID sphereJ, const math::AABB & blockAABB );
   void treatLubricationSphrPlane( const pe::SphereID sphereI, const pe::ConstPlaneID planeJ );

   pe::Vec3 compLubricationSphrSphr ( real_t gap, const pe::SphereID sphereI, const pe::SphereID sphereJ ) const;
   pe::Vec3 compLubricationSphrPlane( real_t gap, const pe::SphereID sphereI, const pe::ConstPlaneID planeJ ) const;

   // member variables
   shared_ptr<StructuredBlockStorage> blockStorage_;
   shared_ptr<pe::BodyStorage> globalBodyStorage_;
   const BlockDataID bodyStorageID_;

   real_t dynamicViscosity_;
   real_t cutOffDistance_;
   real_t minimalGapSize_;

}; // class LubricationCorrection

} // pe_coupling
} // walberla
