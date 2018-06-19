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
//! \file BodyOverlapFunctions.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/AABB.h"
#include "core/math/Vector3.h"


namespace walberla {
namespace geometry {


   //*******************************************************************************************************************
   /*! \brief Body Concept
   *
   *  \file   BodyOverlapFunctions.h
   *
   *  A Body is a 3D object, that is used to initialize fields.
   *  This file defines template functions that can be specialized to define new bodies.
   *
   *  The function that has to be implemented by every body is the contains() function.
   *  It is used by overlapFraction() which uses in the general case cell super-sampling to determine the overlap
   *  of a cell with an arbitrary body.
   *  To speed up this process there are two functions that can be specialized to speed up this super-sampling process:
   *  a fast overlap check between a block and a body, and a fast overlap check between a cell and a body.
   *  The default implementation always returns DONT_KNOW.
   *
   *  A body where the overlap with a cell can be calculated analytically (for example a box) should specialize
   *  the overlapFraction directly, so no super-sampling has to be done.
   *
   *  For an example specializations see Sphere.h, Ellipsoid.h or AABBBody.h
   */
   //*******************************************************************************************************************



   /****************************************************************************************************************//**
   * Returns the overlap fraction between a body and a cell
   *
   * Specialize this function directly for bodies where the overlap can be computed analytically.
   * The default implementation does a cell super-sampling, using the contains() method multiple times.
   *
   * \param body          the body object
   * \param cellMidpoint  midpoint of the cell in global coordinates
   * \param dx            the edge length(s) of the cell, dx or (dx, dy, dz)
   * \param maxDepth      sub sampling depth: the cell edge is divided in half \p maxDepth+1 times.
   *                      Values less than zero result in no subdivisions, making this function behave like contains().
   ********************************************************************************************************************/
   template <typename Body> real_t overlapFraction ( const Body & body, const Vector3<real_t> & cellMidpoint,
                                                     real_t dx, uint_t maxDepth=4 );
   template <typename Body> real_t overlapFraction ( const Body & body, const Vector3<real_t> & cellMidpoint,
                                                     real_t dx, int maxDepth );

   template <typename Body> real_t overlapFraction ( const Body & body, const Vector3<real_t> & cellMidpoint,
                                                     const Vector3<real_t> & dx, uint_t maxDepth=4 );


   /****************************************************************************************************************//**
   * Test if a point in space lies inside a body.
   *
   * This function has to be specialized for every body. There is no default implementation.
   ********************************************************************************************************************/
   template <typename Body> bool contains ( const Body & body, const Vector3<real_t> & point );







   enum FastOverlapResult  {  CONTAINED_INSIDE_BODY, COMPLETELY_OUTSIDE, PARTIAL_OVERLAP, DONT_KNOW };

   /****************************************************************************************************************//**
   * Determines in a fast way ( for example using a bounding box) if a body and a block overlap
   * when no fast computation is possible return DONT_KNOW
   ********************************************************************************************************************/
   template <typename Body> FastOverlapResult fastOverlapCheck ( const Body & body, const AABB & block );


   /****************************************************************************************************************//**
   * Determines in a fast way (bounding box etc) if a body and a block overlap
   * when no fast computation is possible return DONT_KNOW
   ********************************************************************************************************************/
   template <typename Body> FastOverlapResult fastOverlapCheck ( const Body & body, const Vector3<real_t> & cellMidpoint, const Vector3<real_t> & dx );




} // namespace geometry
} // namespace walberla

#include "BodyOverlapFunctions.impl.h"
