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
//! \file MeshBodyWallDistanceCallback.h
//! \ingroup mesh
//! \author Brendan Waters  <brendan.waters@sydney.edu.au>
//! \author Girish Kumatagi <girish.h.kumatagi@fau.de>
//
//======================================================================================================================
#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "core/math/Vector3.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "mesh_common/distance_octree/BranchNode.h"

namespace walberla {
namespace mesh {


//**********************************************************************************************************************
/*!
*   \brief  Class which computes the minimum wall distance between a fluid cell and a mesh object using the 
            Möller-Trumbore Fast Minimum Storage Ray/Triangle Intersection Algorithm

            The operator acts as a callback function which can be bound to Code-Generated boundary conditions requiring 
            the wall distance, i.e. Quadratic Bounce Back.  
*/
//**********************************************************************************************************************

template< typename MeshType >
class MeshBodyWallDistance
{

public:
   MeshBodyWallDistance( const std::shared_ptr< mesh::DistanceOctree<MeshType>> & distanceOctree ) : distanceOctree_(distanceOctree) {}

   real_t operator()(const Cell& fluid, const Cell& boundary, const shared_ptr< StructuredBlockForest >& SbF, IBlock& block) const
   { 
      const uint_t level = SbF->getLevel( block );

      #ifndef NDEBUG
      const real_t dx = SbF->dx( level );

      WALBERLA_ASSERT_FLOAT_EQUAL_2(dx, SbF->dy( level ));
      WALBERLA_ASSERT_FLOAT_EQUAL_2(dx, SbF->dz( level ));
      #endif

      Cell fluidGlobalCell;
      SbF->transformBlockLocalToGlobalCell(fluidGlobalCell, block, fluid);
      const Vector3< real_t > fluidCC = SbF->getCellCenter(fluidGlobalCell, level);

      Cell boundaryGlobalCell;
      SbF->transformBlockLocalToGlobalCell(boundaryGlobalCell, block, boundary);
      const Vector3< real_t > BoundaryCC = SbF->getCellCenter(boundaryGlobalCell, level);

      const typename MeshType::Point ray_origin {fluidCC[0], fluidCC[1], fluidCC[2]};

      typename MeshType::Point ray_direction (  BoundaryCC[0]-fluidCC[0],
                                                BoundaryCC[1]-fluidCC[1],
                                                BoundaryCC[2]-fluidCC[2] );
      
      const real_t ray_length = ray_direction.norm();

      // Returns wall distance from fluidCC to nearest mesh triangle, i.e. |xf - xw|                               
      auto q = distanceOctree_->getRayDistanceToMeshObject(ray_origin, ray_direction.normalize());

      // Interpolated boundary conditions require normalised distance: q = |xf - xw|/|xf - xb|, therefore:
      q /= ray_length;
      
      WALBERLA_ASSERT_GREATER_EQUAL(q, real_c(0.0));
      WALBERLA_ASSERT_LESS_EQUAL(q, real_c(1.0));

      WALBERLA_LOG_DETAIL_ON_ROOT("Normalised wall distance: " << q << " fluid cell center: " << fluidCC  << " boundary cell center: " << BoundaryCC << " dir: " << ray_direction )

      return q;
   };

private:
   const std::shared_ptr< mesh::DistanceOctree<MeshType>> & distanceOctree_;
};

}  // mesh
}  // walberla