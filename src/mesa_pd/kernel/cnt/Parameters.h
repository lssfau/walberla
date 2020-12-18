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
//! \file
//! \author Igor Ostanin <i.ostanin@skoltech.ru>
//! \author Grigorii Drozdov <drozd013@umn.edu>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Constants.h"

namespace walberla {
namespace mesa_pd {
namespace kernel {
namespace cnt {

//====================CNT inertial and elastic parameters========================================================

// Within our model, CNTs are represented with intersecting capsule primitives - cylinders with capped ends.
// Capped ends are necessary for a continuous normal definition
// Each capsule represents the inertial properties of a cylindrical segment.

// For the units system used see PFC 5 implementation docs

/// timestep is fixed to 20 fs.
constexpr auto dT = 20_r;

/// vdW interaction radius w/r to inertial segment radius
constexpr auto CutoffFactor     = 4_r;

// CNT lattice numbers (m,n) CNT
constexpr auto mm = 10_r;
constexpr auto nn = 10_r;

/// CNT Young modulus in GPa based on atomistic simulations
constexpr auto En = 1029_r;
/// CNT shear modulus in GPa based on atomistic simulations
constexpr auto Gs = 459_r;

/// Equilibrium distance of a covalent C-C bond
constexpr auto a_CC = 1.42_r;

/// Equilibrium vdW separation of two CNT surfaces
constexpr auto a_VD = 3.35_r;

/// CNT radius
constexpr auto R_CNT = 6.78_r; //( a_CC / (2 * PII)) * std::sqrt( 3.0 * ( mm * mm + nn * nn + mm * nn ) );

/// external radius of an idealized hollow cylindrical shell
constexpr auto R_1 = R_CNT + 0.5_r * a_VD;
/// internal radius of an idealized hollow cylindrical shell
constexpr auto R_2 = R_CNT - 0.5_r * a_VD;

/// height of a cylindrical segment
constexpr auto T = 2_r * R_CNT;

/// linear density in atoms per A
constexpr auto ro = 4_r * mm / ( math::root_three * a_CC );

/// Atomic mass of Carbon in AMU
constexpr auto M_C = 12.011_r;

/// Mass of the repetative cell in AMU
constexpr auto mass_T = ro * T * M_C * 104.397_r;

/// Volume of a capsule
double vol_capsule = (4_r/3_r) * math::pi * (R_CNT * R_CNT * R_CNT) + math::pi * R_CNT * R_CNT * T;

/// Density of a capsule
double dens_capsule = mass_T / vol_capsule;

/// V-bond parameter
constexpr auto knorm = (En * 0.006242_r) / T;
/// V-bond parameter
constexpr auto kshear = (Gs * 0.006242_r) / T;

constexpr auto margin = 10.e-10_r;

constexpr auto outer_radius     = CutoffFactor * R_CNT;
// real_t inner_radius     = R_CNT; // Capsule
constexpr auto inner_radius     = 1.5811_r * R_CNT; // Sphere

/// volume of a sphere
constexpr auto vol_sphere = (4_r/3_r) * math::pi * (inner_radius * inner_radius * inner_radius);

/// density of a sphere
constexpr auto dens_sphere = mass_T / vol_sphere;

} //namespace cnt
} //namespace kernel
} //namespace mesa_pd
} //namespace walberla
