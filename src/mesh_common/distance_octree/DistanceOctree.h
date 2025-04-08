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
//! \file DistanceOctree.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "Node.h"
#include "LeafNode.h"
#include "BranchNode.h"

#include "mesh_common/DistanceComputations.h"

#include "core/DataTypes.h"
#include "core/Abort.h"

#include "core/math/GenericAABB.h"
#include "core/math/Vector3.h"

#include "core/typeToString.h"


#include <vector>
#include <queue>

namespace walberla {
namespace mesh {
namespace distance_octree {


template <typename MeshType>
class DistanceOctree
{
public:
   using Point = typename MeshType::Point;
   using Normal = typename MeshType::Normal;
   using Scalar = typename MeshType::Scalar;  
   using FaceHandle = typename MeshType::FaceHandle; 
   using AABB = typename math::GenericAABB<Scalar>;

   DistanceOctree( const shared_ptr< TriangleDistance<MeshType> > & triDist, uint_t maxDepth = 20u, uint_t minNumTriangles = 25u )
   {
      if( triDist->getMesh().faces_empty() )
         WALBERLA_ABORT( "You cannot build a distance octree on a mesh without triangles!");

      if( maxDepth == 0 || triDist->getMesh().n_faces() < minNumTriangles )
         rootNode_ = walberla::make_shared<const LeafNode<MeshType> >( triDist, std::vector<FaceHandle>( triDist->getMesh().faces_begin(), triDist->getMesh().faces_end() ) );
      else
         rootNode_ = walberla::make_shared<const BranchNode<MeshType> >( triDist, triDist->getMesh().faces_begin(), triDist->getMesh().faces_end(), maxDepth - 1, minNumTriangles );
   }

   Scalar sqSignedDistance( const Point & p ) const
   {
      WALBERLA_ASSERT_NOT_NULLPTR( rootNode_ );
      return rootNode_->sqSignedDistance( p );
   }

   Scalar sqSignedDistance( const Point & p, FaceHandle & closestTriangle ) const
   {
      return rootNode_->sqSignedDistance( p, closestTriangle );
   }

   Scalar sqSignedDistance( const Point & p, Point & closestPoint ) const
   {
      return rootNode_->sqSignedDistance( p, closestPoint );
   }

   Scalar sqSignedDistance( const Point & p, Point & closestPoint, Normal & normal ) const
   {
      return rootNode_->sqSignedDistance( p, closestPoint, normal );
   }

   Scalar sqDistance( const Point & p ) const
   {
      WALBERLA_ASSERT_NOT_NULLPTR( rootNode_ );
      return rootNode_->sqDistance( p );
   }

   Scalar sqDistance( const Point & p, FaceHandle & closestTriangle ) const
   {
      return rootNode_->sqDistance( p, closestTriangle );
   }

   Scalar sqDistance( const Point & p, Point & closestPoint ) const
   {
      return rootNode_->sqDistance( p, closestPoint );
   }

   Scalar sqDistance( const Point & p, Point & closestPoint, Normal & normal ) const
   {
      return rootNode_->sqDistance( p, closestPoint, normal );
   }

   /**
    * \brief This function calculates the distance from a point to the closest intersection point of a mesh object 
    * along a specified ray (origin and direction). This is done via the  Möller-Trumbore Fast Minimum Storage 
    * Ray/Triangle Intersection Algorithm.  For more information, see:
    *    Möller, T., & Trumbore, B. (1997). Fast, Minimum Storage Ray-Triangle Intersection.
    *    Journal of Graphics Tools, 2(1), 21–28. https://doi.org/10.1080/10867651.1997.10487468
    * 
    * \ingroup mesh_common
    * 
    * \param ray_origin    The origin of the ray as a MeshType::Point.
    * \param normalised_ray_direction The direction of the ray as a MeshType::Point in its normalised state, i.e.
    *                      if the direction vector is:
    *                        ray_direction = [1,1,0]
    *                        normalised_ray_direction = [1/sqrt(2), 1/sqrt(2), 0]
    * 
    * \return The distance to the closest intersection point, or `std::numeric_limits<Scalar>::max()` 
    *         if no intersection occurs.
    * 
    * Usage:
    * \code
    *    const MeshType::Point ray_origin { ... };
    *    MeshType::Point ray_direction { ... };
    * 
    *    auto q = distanceOctree_->getRayDistanceToMeshObject(ray_origin, ray_direction.normalize());
    * \endcode
    * 
    * \warning If the ray and a triangle are parallel (within a small epsilon value), or do not intersect,
    *          the function returns std::numeric_limits<Scalar>::max().
    */
   Scalar getRayDistanceToMeshObject(const Point & ray_origin, const Point & normalised_ray_direction) const
   {
      return rootNode_->getRayDistanceToMeshObject( ray_origin, normalised_ray_direction);
   }

   uint_t numTriangles() const { return rootNode_->numTriangles(); }

   void numTrianglesToStream( std::ostream & os ) { rootNode_->numTrianglesToStream(os, 0); }

   uint_t height() const { return rootNode_->height(); }

   const AABB & getAABB() const { return rootNode_->getAABB(); }

   static inline Point             toPoint( const Vector3<Scalar> & p ) { return toOpenMesh( p ); }
   static inline Vector3<Scalar> fromPoint( const Point & p )           { return toWalberla( p ); }

   static inline Vector3<Scalar> fromNormal( const Normal & p ) { return toWalberla( p ); }

   static inline Scalar   toScalar( const real_t & x ) { return numeric_cast<Scalar>( x ); }
   static inline real_t fromScalar( const Scalar & x ) { return numeric_cast<real_t>( x ); }

   void writeVTKOutput( const std::string & filestem ) const;

protected:
   shared_ptr< const Node<MeshType> > rootNode_;
};



//**********************************************************************************************************************
/*! \brief Write the distance octree to a VTK file.
 * 
 * This method should only be called by the root process:
 * \code
     WALBERLA_ROOT_SECTION()
     {
        distanceOctree->writeVTKOutput("distanceOctree");
     }
 * \endcode
 *
 * \param filestem name of the VTK file without extension
 */
//**********************************************************************************************************************
template <typename MeshType>
void DistanceOctree<MeshType>::writeVTKOutput( const std::string & filestem ) const
{
   std::ostringstream oss;
   oss << filestem << ".vtk";

   std::ofstream outfile( oss.str().c_str() );

   outfile << "# vtk DataFile Version 3.0\n"
           << "DistanceOctree\n"
           << "ASCII\n\n"
           << "DATASET UNSTRUCTURED_GRID\n\n";

   uint_t numNodes = 0;

   std::queue<const Node<MeshType> *> nodeQueue;
   nodeQueue.push( rootNode_.get() );
   while( !nodeQueue.empty() )
   {
      const Node<MeshType> * frontNode = nodeQueue.front();
      nodeQueue.pop();

      if( frontNode->numTriangles() == 0 )
         continue;

      ++numNodes;

      for( uint_t i = 0; i < frontNode->numChildren(); ++i )
         nodeQueue.push( frontNode->getChild( i ) );
   }

   outfile << "POINTS " << ( 8 * numNodes ) << " " << typeToString<real_t>() << "\n\n";


   nodeQueue.push( rootNode_.get() );
   while( !nodeQueue.empty() )
   {
      const Node<MeshType> * frontNode = nodeQueue.front();
      nodeQueue.pop();

      if( frontNode->numTriangles() == 0 )
         continue;

      const auto aabb = frontNode->getAABB();

      for( uint_t z = 0; z != 2; ++z ) {
         for( uint_t y = 0; y != 2; ++y ) {
            for( uint_t x = 0; x != 2; ++x ) {
               outfile << ( ( x == 0 ) ? aabb.xMin() : aabb.xMax() ) << " "
                       << ( ( y == 0 ) ? aabb.yMin() : aabb.yMax() ) << " "
                       << ( ( z == 0 ) ? aabb.zMin() : aabb.zMax() ) << "\n";
            }
         }
      }

      for( uint_t i = 0; i < frontNode->numChildren(); ++i )
         nodeQueue.push( frontNode->getChild( i ) );
   }

   outfile << "\n\nCELLS " << numNodes << " " << ( 9 * numNodes ) << "\n\n";

   for( uint_t i = 0, c = 0; i != numNodes; ++i ) {

      outfile << "8";

      for( uint_t j = 0; j != 8; ++j, ++c )
         outfile << " " << c;

      outfile << std::endl;
   }

   outfile << "\n\nCELL_TYPES " << numNodes << "\n\n";

   for( uint_t i = 0; i != numNodes; ++i )
      outfile << "11\n";

   outfile << "\n\nCELL_DATA " << numNodes;

   outfile << "\n\nSCALARS numTriangles unsigned_int 1"
           << "\nLOOKUP_TABLE default\n";

   nodeQueue.push( rootNode_.get() );
   while( !nodeQueue.empty() )
   {
      const Node<MeshType> * frontNode = nodeQueue.front();
      nodeQueue.pop();

      if( frontNode->numTriangles() == 0 )
         continue;

      outfile  << uint32_c( frontNode->numTriangles() ) << "\n";

      for( uint_t i = 0; i < frontNode->numChildren(); ++i )
         nodeQueue.push( frontNode->getChild( i ) );
   }

   outfile << "\n\nSCALARS height unsigned_char 1"
           << "\nLOOKUP_TABLE default\n";
   
   nodeQueue.push( rootNode_.get() );
   while( !nodeQueue.empty() )
   {
      const Node<MeshType> * frontNode = nodeQueue.front();
      nodeQueue.pop();

      if( frontNode->numTriangles() == 0 )
         continue;
      
      WALBERLA_ASSERT_LESS_EQUAL( frontNode->height(), std::numeric_limits<uint8_t>::max() );

      outfile  << uint16_c( frontNode->height() ) << "\n";

      for( uint_t i = 0; i < frontNode->numChildren(); ++i )
         nodeQueue.push( frontNode->getChild( i ) );
   }

   outfile << "\n\nSCALARS depth unsigned_char 1"
           << "\nLOOKUP_TABLE default\n";

   std::queue<uint8_t> depthQueue;
   nodeQueue.push( rootNode_.get() );
   depthQueue.push( 0 );
   while( !nodeQueue.empty() )
   {
      const Node<MeshType> * frontNode = nodeQueue.front();
      uint8_t const depth = depthQueue.front();
      nodeQueue.pop();
      depthQueue.pop();

      if( frontNode->numTriangles() == 0 )
         continue;

      outfile  << uint16_c( depth ) << "\n";

      for( uint_t i = 0; i < frontNode->numChildren(); ++i )
      {
         nodeQueue.push( frontNode->getChild( i ) );
         depthQueue.push( uint8_c( depth + 1 ) );
      }
   }

   outfile << std::endl;
   outfile.close();
}

} //namespace distance_octree

using distance_octree::DistanceOctree;

} // namespace mesh
} // namespace walberla
