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
//! \file Thresholds.cpp
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Numerical thresholds for the physics engine
//
//======================================================================================================================

#include "Thresholds.h"

namespace walberla {
namespace pe {

//=================================================================================================
//
//  GLOBAL THRESHOLD VALUES
//
//=================================================================================================

//*************************************************************************************************
//! Threshold for the contact classification.
/*! This threshold separates between separating, resting and colliding contacts. */
real_t collisionThreshold = Thresholds<real_t>::collisionThreshold();
//*************************************************************************************************


//*************************************************************************************************
//! Threshold for the distance between two rigid bodies.
/*! Rigid bodies with a distance smaller than this threshold are in contact. */
real_t contactThreshold = Thresholds<real_t>::contactThreshold();
//*************************************************************************************************


//*************************************************************************************************
//! Threshold for the restriction of the coefficient of restitution.
/*! In case the relative velocity between two colliding rigid bodies is smaller than this
    threshold, a coefficient of restitution of 0 is used to avoid an infinite number of
    collisions during a single time step. */
real_t restitutionThreshold = Thresholds<real_t>::restitutionThreshold();
//*************************************************************************************************


//*************************************************************************************************
//! Threshold for the separation between static and dynamic friction.
/*! This threshold represents the boundary between static and dynamic friction. */
real_t frictionThreshold = Thresholds<real_t>::frictionThreshold();
//*************************************************************************************************


//*************************************************************************************************
//! Threshold for surface points/checks.
/*! Only points with a distance to the surface smaller than this threshold are considered
    surface point. */
real_t surfaceThreshold = Thresholds<real_t>::surfaceThreshold();
//*************************************************************************************************

} // namespace pe
}  // namespace walberla
