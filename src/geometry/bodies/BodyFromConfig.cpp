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
//! \file BodyFromConfig.cpp
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "BodyFromConfig.h"

#include <memory>
#include "core/Abort.h"
#include "core/StringUtility.h"


namespace walberla {
namespace geometry {


Sphere sphereFromConfig( const Config::BlockHandle & block )
{
   if( ! block.isDefined( "midpoint" ) )
      WALBERLA_ABORT(  "Missing parameter 'midpoint' for sphere defined in block "  << block.getKey()  );
   if ( ! block.isDefined( "radius" ) )
      WALBERLA_ABORT( "Missing parameter 'radius' for sphere defined in block " << block.getKey() );

   return Sphere( block.getParameter<Vector3<real_t> > ( "midpoint" ),
                  block.getParameter<real_t          > ( "radius"   ) );
}


Cylinder cylinderFromConfig( const Config::BlockHandle & block )
{
   if( ! block.isDefined( "min" ) )
      WALBERLA_ABORT( "Missing parameter 'min' for cylinder defined in block " << block.getKey() );
   if( ! block.isDefined( "max" ) )
      WALBERLA_ABORT( "Missing parameter 'max' for cylinder defined in block " << block.getKey() );
   if ( ! block.isDefined( "radius" ) )
      WALBERLA_ABORT( "Missing parameter 'radius' for cylinder defined in block " << block.getKey() );

   Vector3<real_t> min = block.getParameter<Vector3<real_t> > ( "min" );
   Vector3<real_t> max = block.getParameter<Vector3<real_t> > ( "max" );

   return Cylinder( min, max, block.getParameter<real_t>( "radius" ) );
}


Torus torusFromConfig( const Config::BlockHandle & block )
{
   if( ! block.isDefined( "midpoint" ) )
      WALBERLA_ABORT( "Missing parameter 'midpoint' for torus defined in block " << block.getKey() );
   if( ! block.isDefined( "normal" ) )
      WALBERLA_ABORT( "Missing parameter 'normal' for torus defined in block " << block.getKey() );
   if ( ! block.isDefined( "radius" ) )
      WALBERLA_ABORT( "Missing parameter 'radius' for torus defined in block " << block.getKey() );
   if ( ! block.isDefined( "distance" ) )
      WALBERLA_ABORT( "Missing parameter 'distance' for torus defined in block " << block.getKey() );

   Vector3<real_t> midpoint = block.getParameter<Vector3<real_t> > ( "midpoint" );
   Vector3<real_t> normal   = block.getParameter<Vector3<real_t> > ( "normal"   );

   return Torus( midpoint, normal, block.getParameter<real_t>( "radius" ), block.getParameter<real_t>( "distance" ) );
}


Ellipsoid ellipsoidFromConfig ( const Config::BlockHandle & block )
{
   if( ! block.isDefined( "midpoint" ) )
      WALBERLA_ABORT( "Missing parameter 'midpoint' for ellipsoid defined in block " << block.getKey() );
   if( ! block.isDefined( "radii" ) )
      WALBERLA_ABORT( "Missing parameter 'radii' for ellipsoid defined in block " << block.getKey() );

   Vector3<real_t> radii = block.getParameter<Vector3<real_t> > ( "radii" );

   return Ellipsoid( block.getParameter<Vector3<real_t> >( "midpoint" ),
                     block.getParameter<Vector3<real_t> >( "axis1", Vector3<real_t>( real_t(1),real_t(0),real_t(0) ) ),
                     block.getParameter<Vector3<real_t> >( "axis2", Vector3<real_t>( real_t(0),real_t(1),real_t(0) ) ),
                     radii );
}


AABB AABBFromConfig( const Config::BlockHandle & block )
{
   if( ! block.isDefined( "min" ) )
      WALBERLA_ABORT( "Missing parameter 'min' for box defined in block " << block.getKey() );
   if( ! block.isDefined( "max" ) )
      WALBERLA_ABORT( "Missing parameter 'max' for box defined in block " << block.getKey() );

   Vector3<real_t> min = block.getParameter<Vector3<real_t> > ( "min" );
   Vector3<real_t> max = block.getParameter<Vector3<real_t> > ( "max" );

   return AABB( min[0], min[1], min[2], max[0], max[1], max[2] );
}


BodyLogicalAND<Sphere,AABB> sphereSliceFromConfig ( const Config::BlockHandle & block )
{
   Sphere sphere = sphereFromConfig(block.getOneBlock("Sphere"));
   AABB box = AABBFromConfig  (block.getOneBlock("Box"));
   auto spherePtr = std::make_shared<Sphere>( sphere );
   auto boxPtr    = std::make_shared<AABB>  ( box );
   
   return BodyLogicalAND<Sphere,AABB>(spherePtr, boxPtr);
}

BodyLogicalAND<Sphere,BodyLogicalNOT<Sphere> > hollowSphereFromConfig ( const Config::BlockHandle & block )
{
   if ( ! block.isDefined( "midpoint" ) )
      WALBERLA_ABORT(  "Missing parameter 'midpoint' for sphere defined in block "  << block.getKey()  );
   if ( ! block.isDefined( "inner_radius" ) )
      WALBERLA_ABORT( "Missing parameter 'inner_radius' for sphere defined in block " << block.getKey() );
   if ( ! block.isDefined( "outer_radius" ) )
      WALBERLA_ABORT( "Missing parameter 'outer_radius' for sphere defined in block " << block.getKey() );

   auto inner     = std::make_shared<Sphere>(  block.getParameter<Vector3<real_t> > ( "midpoint" ),
                  							              block.getParameter<real_t          > ( "inner_radius"   ) );
   auto outer     = std::make_shared<Sphere>(  block.getParameter<Vector3<real_t> > ( "midpoint" ),
                  			                           block.getParameter<real_t          > ( "outer_radius"   ) );

   auto not_inner = make_shared<BodyLogicalNOT<Sphere> >(inner);
	
   return BodyLogicalAND< Sphere,BodyLogicalNOT<Sphere> >(outer, not_inner);
}

shared_ptr<AbstractBody> bodyFromConfig (const Config::BlockHandle & block )
{
   std::string shape = block.getParameter<std::string>("shape");
   
   if      ( string_icompare( shape, "Sphere"   ) == 0 )
      return make_DynamicBody(sphereFromConfig   ( block ) );
   else if ( string_icompare( shape, "Cylinder" ) == 0 )
      return make_DynamicBody(cylinderFromConfig ( block ) );
   else if ( string_icompare( shape, "Torus"    ) == 0 )
      return make_DynamicBody(torusFromConfig    ( block ) );
   else if ( string_icompare( shape, "Ellipsoid") == 0 )
      return make_DynamicBody(ellipsoidFromConfig( block ) );
   else if ( string_icompare( shape, "Box"      ) == 0 )
      return make_DynamicBody(AABBFromConfig     ( block ) );
   else if ( string_icompare( shape, "SphereSlice") == 0 )
      return make_DynamicBody(sphereSliceFromConfig   ( block ) );
   else if ( string_icompare( shape, "HollowSphere") == 0 )
      return make_DynamicBody(hollowSphereFromConfig  ( block ) );
   else if ( string_icompare( shape, "LogicalAND") == 0 )
      return bodyANDFromConfig ( block );
   else if ( string_icompare( shape, "LogicalOR") == 0 )
      return bodyORFromConfig  ( block );
   else if ( string_icompare( shape, "LogicalXOR") == 0 )
      return bodyXORFromConfig  ( block );
   else if ( string_icompare( shape, "LogicalNOT") == 0 )
      return bodyNOTFromConfig  ( block );

   WALBERLA_ABORT( "Unknown Block " << block.getKey() << "\nAllowed blocks are 'Sphere', 'Cylinder', 'Torus' and 'Ellipsoid'" );

#ifdef __IBMCPP__
   return shared_ptr<AbstractBody>(); // never reached, helps to suppress a warning from the IBM compiler
#endif
}

shared_ptr<AbstractBody> bodyANDFromConfig (const Config::BlockHandle & block )
{
   auto first  = bodyFromConfig( block.getOneBlock("first"));
   auto second = bodyFromConfig( block.getOneBlock("second"));
   auto body = BodyLogicalAND<AbstractBody, AbstractBody>(first,second);
   
   return make_DynamicBody(body);
}
shared_ptr<AbstractBody> bodyORFromConfig  (const Config::BlockHandle & block )
{
   auto first = bodyFromConfig( block.getOneBlock("first"));
   auto second = bodyFromConfig( block.getOneBlock("second"));
   auto body = BodyLogicalOR <AbstractBody, AbstractBody>(first,second);
   
   return make_DynamicBody(body);
}
shared_ptr<AbstractBody> bodyXORFromConfig (const Config::BlockHandle & block )
{
   auto first = bodyFromConfig( block.getOneBlock("first"));
   auto second = bodyFromConfig( block.getOneBlock("second"));
   auto body = BodyLogicalXOR<AbstractBody, AbstractBody>(first,second);
   
   return make_DynamicBody(body);
}
shared_ptr<AbstractBody> bodyNOTFromConfig (const Config::BlockHandle & block )
{
   auto first = bodyFromConfig( block.getOneBlock("first"));
   auto body = BodyLogicalNOT<AbstractBody>(first);
   
   return make_DynamicBody(body);
}


} // namespace geometry
} // namespace walberla
