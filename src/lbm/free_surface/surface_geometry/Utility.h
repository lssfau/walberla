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
//! \file Utility.h
//! \ingroup surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Helper functions for surface geometry computations.
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"

#include "ContactAngle.h"

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * Compute the point p that lays on the plane of the interface surface (see Algorithm 2.1 in dissertation of Thomas
 * Pohl, 2008). p = (0.5, 0.5, 0.5) + offset * normal.
 **********************************************************************************************************************/
Vector3< real_t > getInterfacePoint(const Vector3< real_t >& normal, real_t fillLevel);

/***********************************************************************************************************************
 * Compute the intersection point of a surface (defined by normal and surfacePoint) with an edge of a cell.
 **********************************************************************************************************************/
real_t getCellEdgeIntersection(const Vector3< real_t >& edgePoint, const Vector3< real_t >& edgeDirection,
                               const Vector3< real_t >& normal, const Vector3< real_t >& surfacePoint);

/***********************************************************************************************************************
 * Compute the fluid volume within an interface cell with respect to
 * - the interface normal
 * - the interface surface offset.
 *
 * see dissertation of T. Pohl, 2008, section 2.5.3, p. 23-26.
 **********************************************************************************************************************/
real_t computeCellFluidVolume(const Vector3< real_t >& normal, real_t offset);

/***********************************************************************************************************************
 * Compute an artificial wall point according to the artifical curvature wetting model from the dissertation of Stefan
 * Donath, 2011. The artificial wall point and the artificial normal can be used to alter the curvature computation with
 * local triangulation. The interface curvature is changed such that the correct laplace pressure with respect to the
 * contact angle is assumed near solid cells.
 *
 * see dissertation of T. Pohl, 2008, section 6.3.3
 **********************************************************************************************************************/
bool computeArtificalWallPoint(const Vector3< real_t >& globalInterfacePointLocation,
                               const Vector3< real_t >& globalCellCoordinate, const Vector3< real_t >& normal,
                               const Vector3< real_t >& wallNormal, const Vector3< real_t >& obstacleNormal,
                               const ContactAngle& contactAngle, Vector3< real_t >& artificialWallPointCoord,
                               Vector3< real_t >& artificialWallNormal);
} // namespace free_surface
} // namespace walberla
