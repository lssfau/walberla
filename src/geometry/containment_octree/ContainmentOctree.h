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
//! \file ContainmentOctree.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "Node.h"
#include "IndeterminateLeafNode.h"
#include "InsideLeafNode.h"
#include "OutsideLeafNode.h"
#include "BranchNode.h"

#include "core/DataTypes.h"
#include "core/math/GenericAABB.h"
#include "core/math/KahanSummation.h"
#include "core/typeToString.h"

#include <queue>
#include <vector>

namespace walberla {
namespace geometry {

template< typename DistanceObjectT >
class ContainmentOctree
{
public:
   typedef typename DistanceObjectT::Scalar Scalar;
   typedef typename DistanceObjectT::Point  Point;
   typedef math::GenericAABB<Scalar>        AABB;

   typedef math::KahanAccumulator<Scalar> KahanAccumulator;
   typedef DistanceObjectT                DistanceObject;

   typedef containment_octree::Node<ContainmentOctree>                  Node;
   typedef containment_octree::InsideLeafNode<ContainmentOctree>        InsideLeafNode;
   typedef containment_octree::OutsideLeafNode<ContainmentOctree>       OutsideLeafNode;
   typedef containment_octree::IndeterminateLeafNode<ContainmentOctree> IndeterminateLeafNode;
   typedef containment_octree::BranchNode<ContainmentOctree>            BranchNode;

   inline ContainmentOctree( const shared_ptr<const DistanceObject> & distanceObject, const Scalar epsilon = Scalar(0),
                             const uint_t maxDepth = 6u, const Scalar minAABBVolume = Scalar(0) );
    
   inline ContainmentOctree( const shared_ptr<const DistanceObject> & distanceObject, const AABB & aabb,
                             const Scalar epsilon = Scalar(0), const uint_t maxDepth = 6u,
                             const Scalar minAABBVolume = Scalar(0) );


   bool     contains        ( const Point & p ) const { return aabb_.contains( fromPoint(p) ) && rootNode_->contains( p ); }
   Scalar   sqSignedDistance( const Point & p ) const { return distanceObject_->sqSignedDistance( p ); }

   uint_t height() const { return rootNode_->height(); }
   uint_t numNodes() const { return rootNode_->numNodes(); }
   void numNodes( uint_t & numInside, uint_t & numOutside, uint_t & numIndeterminate, uint_t & numBranch ) const;
   void volumes( Scalar & insideVolume, Scalar & outsideVolume, Scalar & indeterminateVolume ) const;
   size_t memory() const;

   const DistanceObjectT & getDistanceObject() const { return *distanceObject_; }

   AABB getAABB() const { return aabb_; }

   inline void writeVTKOutput( const std::string & filestem ) const;


   static inline Point             toPoint( const Vector3<Scalar> & p ) { return DistanceObjectT::toPoint( p ); }
   static inline Vector3<Scalar> fromPoint( const Point & p )           { return DistanceObjectT::fromPoint( p ); }

   static inline Scalar   toScalar( const real_t & x ) { return DistanceObjectT::toScalar( x ); }
   static inline real_t fromScalar( const Scalar & x ) { return DistanceObjectT::fromScalar( x ); }

private:
   inline void init( Scalar epsilon, const uint_t maxDepth, const Scalar minBoxVolume );

protected:
   shared_ptr<const Node > rootNode_;
   shared_ptr<const DistanceObjectT>          distanceObject_;
   AABB                                       aabb_;
};


template< typename DistanceObjectT >
ContainmentOctree<DistanceObjectT>::ContainmentOctree( const shared_ptr<const DistanceObjectT> & distanceObject, const Scalar epsilon /*= Scalar(0)*/,
                                                       const uint_t maxDepth /*= 6u*/, const Scalar minBoxVolume /*= Scalar(0)*/ )
   : distanceObject_( distanceObject ), aabb_( distanceObject->getAABB() )
{
   init( epsilon, maxDepth, minBoxVolume );
}


template< typename DistanceObjectT >
ContainmentOctree<DistanceObjectT>::ContainmentOctree( const shared_ptr<const DistanceObjectT> & distanceObject, const AABB & aabb,
   const Scalar epsilon /*= Scalar(0)*/, const uint_t maxDepth /*= 6u*/,
   const Scalar minBoxVolume /*= Scalar(0)*/ )
   : distanceObject_(  distanceObject ), aabb_( aabb )
{
   init( epsilon, maxDepth, minBoxVolume );
}

template< typename DistanceObjectT >
void ContainmentOctree<DistanceObjectT>::init( Scalar epsilon, const uint_t maxDepth, const Scalar minBoxVolume )
{
   aabb_.extend( epsilon );
   auto aabbDimensions = aabb_.sizes();
   aabbDimensions *= Scalar(0.5);
   const Scalar circumcircleRadius = aabbDimensions.sqrLength();

   const auto aabbCenter = aabb_.center();

   const Scalar sqSignedDist = distanceObject_->sqSignedDistance( toPoint( aabbCenter ) );

   if( std::fabs( sqSignedDist ) > circumcircleRadius )
   {
      if( sqSignedDist < Scalar(0) )
      {
         rootNode_ = make_shared<InsideLeafNode>();
      }
      else
      {
         rootNode_ = make_shared<OutsideLeafNode>();
      }
   }
   else
   {
      if( maxDepth == 0 || aabb_.volume() < minBoxVolume )
      {
         rootNode_ = make_shared<IndeterminateLeafNode>( distanceObject_, epsilon );
      }
      else
      {
         rootNode_ = make_shared<BranchNode>( distanceObject_, aabb_, epsilon, maxDepth - 1, minBoxVolume );
      }
   }
}


template< typename DistanceObjectT >
void ContainmentOctree<DistanceObjectT>::numNodes( uint_t & numInside, uint_t & numOutside, uint_t & numIndeterminate, uint_t & numBranch ) const
{
   numInside = numOutside = numIndeterminate = numBranch = 0;

   rootNode_->numNodes( numInside, numOutside, numIndeterminate, numBranch );
}


template< typename DistanceObjectT >
void ContainmentOctree<DistanceObjectT>::volumes( Scalar & insideVolume, Scalar & outsideVolume, Scalar & indeterminateVolume ) const
{
   KahanAccumulator insideVolumeKahan, outsideVolumeKahan, indeterminateVolumeKahan;

   rootNode_->volumes( insideVolumeKahan, outsideVolumeKahan, indeterminateVolumeKahan, aabb_.volume() );

   insideVolume = insideVolumeKahan.get();
   outsideVolume = outsideVolumeKahan.get();
   indeterminateVolume = indeterminateVolumeKahan.get();
}


template< typename DistanceObjectT >
void ContainmentOctree<DistanceObjectT>::writeVTKOutput( const std::string & filestem ) const
{
   std::ostringstream oss;
   oss << filestem << ".vtk";

   std::ofstream outfile( oss.str().c_str() );

   outfile << "# vtk DataFile Version 3.0\n"
           << "ContainmentOctree\n"
           << "ASCII\n\n"
           << "DATASET UNSTRUCTURED_GRID\n\n";

   uint_t numberOfNodes = 0;

   std::queue<const Node *> nodeQueue;
   nodeQueue.push( rootNode_.get() );
   while( !nodeQueue.empty() )
   {
      const Node * frontNode = nodeQueue.front();
      nodeQueue.pop();

      ++numberOfNodes;

      for( uint_t i = 0; i < frontNode->numChildren(); ++i )
         nodeQueue.push( frontNode->getChild( i ) );
   }

   outfile << "POINTS " << ( 8 * numberOfNodes ) << " " << typeToString<real_t>() << "\n\n";


   std::queue<AABB> aabbQueue;
   nodeQueue.push( rootNode_.get() );
   aabbQueue.push( aabb_ );
   while( !nodeQueue.empty() )
   {
      const Node * frontNode = nodeQueue.front();
      AABB aabb = aabbQueue.front();
      nodeQueue.pop();
      aabbQueue.pop();

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

      if( frontNode->numChildren() > 0 )
      {
         const BranchNode * branchNode = dynamic_cast<const BranchNode *>( frontNode );
         WALBERLA_ASSERT_NOT_NULLPTR( branchNode );
         const auto &    min = aabb.minCorner();
         const auto &    max = aabb.maxCorner();
         const auto & center = branchNode->center();

         aabbQueue.push( AABB::createFromMinMaxCorner(    min[0],    min[1],    min[2], center[0], center[1], center[2] ) );
         aabbQueue.push( AABB::createFromMinMaxCorner(    min[0],    min[1], center[2], center[0], center[1],    max[2] ) );
         aabbQueue.push( AABB::createFromMinMaxCorner(    min[0], center[1],    min[2], center[0],    max[1], center[2] ) );
         aabbQueue.push( AABB::createFromMinMaxCorner(    min[0], center[1], center[2], center[0],    max[1],    max[2] ) );
         aabbQueue.push( AABB::createFromMinMaxCorner( center[0],    min[1],    min[2],    max[0], center[1], center[2] ) );
         aabbQueue.push( AABB::createFromMinMaxCorner( center[0],    min[1], center[2],    max[0], center[1],    max[2] ) );
         aabbQueue.push( AABB::createFromMinMaxCorner( center[0], center[1],    min[2],    max[0],    max[1], center[2] ) );
         aabbQueue.push( AABB::createFromMinMaxCorner( center[0], center[1], center[2],    max[0],    max[1],    max[2] ) );
      }
   }

   outfile << "\n\nCELLS " << numberOfNodes << " " << ( 9 * numberOfNodes ) << "\n\n";

   for( uint_t i = 0, c = 0; i != numberOfNodes; ++i ) {

      outfile << "8";

      for( uint_t j = 0; j != 8; ++j, ++c )
         outfile << " " << c;

      outfile << std::endl;
   }

   outfile << "\n\nCELL_TYPES " << numberOfNodes << "\n\n";

   for( uint_t i = 0; i != numberOfNodes; ++i )
      outfile << "11\n";

   outfile << "\n\nCELL_DATA " << numberOfNodes;

   outfile << "\n\nSCALARS nodeType unsigned_char 1"
      << "\nLOOKUP_TABLE default\n";

   nodeQueue.push( rootNode_.get() );
   while( !nodeQueue.empty() )
   {
      const Node * frontNode = nodeQueue.front();
      nodeQueue.pop();

      if( dynamic_cast<const OutsideLeafNode *>( frontNode ) != 0 )
         outfile  << "1\n";
      else if( dynamic_cast<const InsideLeafNode *>( frontNode ) != 0 )
         outfile  << "2\n";
      else if( dynamic_cast<const IndeterminateLeafNode *>( frontNode ) != 0 )
         outfile  << "3\n";
      else
      {
         WALBERLA_ASSERT_NOT_NULLPTR( dynamic_cast<const BranchNode *>( frontNode ) );
         outfile  << "0\n";
      }

      for( uint_t i = 0; i < frontNode->numChildren(); ++i )
         nodeQueue.push( frontNode->getChild( i ) );
   }

   outfile << "\n\nSCALARS height unsigned_char 1"
           << "\nLOOKUP_TABLE default\n";

   nodeQueue.push( rootNode_.get() );
   while( !nodeQueue.empty() )
   {
      const Node * frontNode = nodeQueue.front();
      nodeQueue.pop();

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
      const Node * frontNode = nodeQueue.front();
      uint8_t depth = depthQueue.front();
      nodeQueue.pop();
      depthQueue.pop();

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

} // namespace geometry
} // namespace walberla
