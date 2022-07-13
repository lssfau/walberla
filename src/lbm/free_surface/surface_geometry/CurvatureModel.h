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
//! \file CurvatureModel.h
//! \ingroup surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Collection of sweeps required for using a specific curvature model.
//
//======================================================================================================================

#pragma once

namespace walberla
{
namespace free_surface
{
// forward declaration
template< typename LatticeModel_T, typename FlagField_T, typename ScalarField_T, typename VectorField_T >
class SurfaceGeometryHandler;

namespace curvature_model
{
/***********************************************************************************************************************
 * Collection of sweeps for computing the curvature using a finite difference-based (Parker-Youngs) approach according
 * to:
 * dissertation of S. Bogner, 2017 (section 4.4.2.1)
 **********************************************************************************************************************/
template< typename Stencil_T, typename LatticeModel_T, typename FlagField_T, typename ScalarField_T,
          typename VectorField_T >
class FiniteDifferenceMethod
{
 private:
   using SurfaceGeometryHandler_T =
      walberla::free_surface::SurfaceGeometryHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >;

 public:
   void addSweeps(SweepTimeloop& timeloop, const SurfaceGeometryHandler_T& geometryHandler);
}; // class FiniteDifferenceMethod

/***********************************************************************************************************************
 * Collection of sweeps for computing the curvature using local triangulation according to:
 * - dissertation of T. Pohl, 2008 (section 2.5)
 * - dissertation of S. Donath, 2011 (wetting model, section 6.3.3)
 **********************************************************************************************************************/
template< typename Stencil_T, typename LatticeModel_T, typename FlagField_T, typename ScalarField_T,
          typename VectorField_T >
class LocalTriangulation
{
 private:
   using SurfaceGeometryHandler_T =
      walberla::free_surface::SurfaceGeometryHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >;

 public:
   void addSweeps(SweepTimeloop& timeloop, const SurfaceGeometryHandler_T& geometryHandler);
}; // class LocalTriangulation

/***********************************************************************************************************************
 * Collection of sweeps for computing the curvature with a simplistic finite difference method. This approach is not
 * documented in literature and neither thoroughly tested or validated.
 * Use it with caution and preferably for testing purposes only.
 **********************************************************************************************************************/
template< typename Stencil_T, typename LatticeModel_T, typename FlagField_T, typename ScalarField_T,
          typename VectorField_T >
class SimpleFiniteDifferenceMethod
{
 private:
   using SurfaceGeometryHandler_T =
      walberla::free_surface::SurfaceGeometryHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T >;

 public:
   void addSweeps(SweepTimeloop& timeloop, const SurfaceGeometryHandler_T& geometryHandler);
}; // class SimpleFiniteDifferenceMethod

} // namespace curvature_model
} // namespace free_surface
} // namespace walberla

#include "CurvatureModel.impl.h"