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
//! \file BodyFromConfig.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Functions to create a body from a configuration block
//
//======================================================================================================================

#pragma once

#include "Cylinder.h"
#include "Ellipsoid.h"
#include "Sphere.h"
#include "Torus.h"
#include "BodyLogic.h"
#include "DynamicBody.h"
#include "AABBBody.h"
#include "core/config/Config.h"
#include "core/math/AABB.h"


namespace walberla {
namespace geometry {


   /****************************************************************************************************************//**
   * Parses a configuration block and returns a Sphere
   *
   * Example block:  Sphere{ midpoint <1,2,3>; radius 5; }
   ********************************************************************************************************************/
   Sphere sphereFromConfig( const Config::BlockHandle & block );


   /****************************************************************************************************************//**
   * Parses a configuration block and returns a Cylinder
   *
   * Example block:  Cylinder{ min <1,2,3>; max <5,4,3>; radius 3; }
   ********************************************************************************************************************/
   Cylinder cylinderFromConfig( const Config::BlockHandle & block );


   /****************************************************************************************************************//**
   * Parses a configuration block and returns a Torus
   *
   * Example block:  Torus{ midpoint <1,2,3>; normal <1,0,0>; radius 5; distance 1; }
   ********************************************************************************************************************/
   Torus torusFromConfig( const Config::BlockHandle & block );


   /****************************************************************************************************************//**
   * Parses a configuration block and returns an Ellipsoid
   *
   * Example block:  Ellipsoid{ midpoint <1,2,3>; radii <5,4,3>; axis1  <1,0,0>; axis2 <0,1,0> }
   *
   * For the semantic of the parameters see constructor of Ellipsoid
   ********************************************************************************************************************/
   Ellipsoid ellipsoidFromConfig ( const Config::BlockHandle & block );


   /****************************************************************************************************************//**
   * Parses a configuration block and returns a box
   *
   * Example block:  box{ min <1,2,3>; max <5,4,3>; }
   ********************************************************************************************************************/
   AABB AABBFromConfig( const Config::BlockHandle & block );

   /****************************************************************************************************************//**
   * Parses a configuration block and returns a difference/union/... of bodies.
   ********************************************************************************************************************/
   BodyLogicalAND<Sphere,AABB>                    sphereSliceFromConfig  ( const Config::BlockHandle & block );
   BodyLogicalAND<Sphere,BodyLogicalNOT<Sphere> > hollowSphereFromConfig ( const Config::BlockHandle & block );

   /****************************************************************************************************************//**
   * Parses a configuration block and returns a body via dynamic polymorphism.
   ********************************************************************************************************************/
   shared_ptr<AbstractBody>                       bodyFromConfig    ( const Config::BlockHandle & block );
   shared_ptr<AbstractBody>                       bodyANDFromConfig ( const Config::BlockHandle & block );
   shared_ptr<AbstractBody>                       bodyORFromConfig  ( const Config::BlockHandle & block );
   shared_ptr<AbstractBody>                       bodyXORFromConfig ( const Config::BlockHandle & block );
   shared_ptr<AbstractBody>                       bodyNOTFromConfig ( const Config::BlockHandle & block );

} // namespace geometry
} // namespace walberla

