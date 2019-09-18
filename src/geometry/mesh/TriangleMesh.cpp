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
//! \file TriangleMesh.cpp
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "TriangleMesh.h"
#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/math/AABB.h"

#include <algorithm>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <vector>


namespace walberla {
namespace geometry {

using std::set;
using std::vector;
using std::map;


//======================================================================================================================
//
//  INSERTING
//
//======================================================================================================================

//**********************************************************************************************************************
/*! \brief Adds a vertex to the mesh and returns its index
 *  No duplicate check is done when inserting, instead use removeDuplicateVertices
 */
//*********************************************************************************************************************/
TriangleMesh::index_t TriangleMesh::addVertex(const vertex_t & v)
{
   vertices_.push_back(v);
   WALBERLA_ASSERT_LESS( vertices_.size(), std::numeric_limits<index_t>::max() );
   return index_t(vertices_.size() -1);
}

TriangleMesh::index_t TriangleMesh::addVertexNormal(const vertex_t & n)
{
   vertexNormals_.push_back(n);
   WALBERLA_ASSERT_LESS( vertexNormals_.size(), std::numeric_limits<index_t>::max() );
   return index_t(vertexNormals_.size() -1);
}

TriangleMesh::index_t TriangleMesh::addVertex(const vertex_t & v, const color_t & c)
{
   vertices_.push_back(v);
   vertexColors_.push_back(c);
   WALBERLA_ASSERT_LESS( vertices_.size(), std::numeric_limits<index_t>::max() );
   WALBERLA_ASSERT_EQUAL( vertices_.size(), vertexColors_.size() );
   return index_t(vertices_.size() -1);
}

//**********************************************************************************************************************
/*! \brief Adds a triangle to the mesh, as parameters use the return values of addVertex()
 */
//*********************************************************************************************************************/
void TriangleMesh::addTriangle( index_t v0Idx, index_t v1Idx, index_t v2Idx )
{
   WALBERLA_ASSERT_LESS(v0Idx, vertices_.size());
   WALBERLA_ASSERT_LESS(v1Idx, vertices_.size());
   WALBERLA_ASSERT_LESS(v2Idx, vertices_.size());

   vertexIndices_.push_back(v0Idx);
   vertexIndices_.push_back(v1Idx);
   vertexIndices_.push_back(v2Idx);
}

void TriangleMesh::addTriangle(index_t v0Idx, index_t v1Idx, index_t v2Idx, index_t n0Idx, index_t n1Idx, index_t n2Idx)
{
   WALBERLA_ASSERT_LESS( n0Idx, vertexNormals_.size() );
   WALBERLA_ASSERT_LESS( n1Idx, vertexNormals_.size() );
   WALBERLA_ASSERT_LESS( n2Idx, vertexNormals_.size() );

   addTriangle( v0Idx, v1Idx, v2Idx );

   normalIndices_.push_back(n0Idx);
   normalIndices_.push_back(n1Idx);
   normalIndices_.push_back(n2Idx);

   WALBERLA_ASSERT_EQUAL( normalIndices_.size(), vertexIndices_.size() );
}



//======================================================================================================================
//
//  MESH OPERATIONS
//
//======================================================================================================================



//**********************************************************************************************************************
/*! \brief Helper class, used in TriangleMesh::removeDuplicateVertices
 *
 * The class can be used as a comparator object. It compares to integers,
 * by interpreting them as indices in an array of vertices (which is given at construction).
 * Then the positions of these vertices is compared.
 * This is used to sort an array of indices, rather than the vertex array itself.
 */
//*********************************************************************************************************************/

class LessThan
{
   public:
      LessThan(const vector< Vector3<real_t> > & v, real_t tol)
         : vec(v), tolerance(tol) {}

      inline bool operator()( const size_t i1, const size_t i2 ) const
      {
         for(uint_t i=0; i<3; ++i)
         {
            if( realIsEqual( vec[i1][i], vec[i2][i], tolerance ) )
               continue;

            return ( vec[i1][i] < vec[i2][i] );
         }
         return false; //vertices are equal
      }

   private:
      const vector< Vector3<real_t> > & vec;
      real_t tolerance;
};



//**********************************************************************************************************************
/*! \brief Removes duplicate vertices i.e. vertices with same positions
 *
 * This is expensive, since the vertex array has to be sorted, scanned for duplicates, and the
 * face indices array has to be rewritten.
 *
 * \warning This function changes vertex indices! Consider this when adding triangles.
 * \param tolerance  two vertices are considered equal if each of their coordinates does not
 *                   differ more than the tolerance value (maximum-norm)
 * \return The number of removed vertices
 */
//*********************************************************************************************************************/
size_t TriangleMesh::removeDuplicateVertices( real_t tolerance )
{
   if(vertices_.size() <= 1)
      return 0;

   // create an index array, containing values from 0 to vertices_.size()
   vector<size_t> vInd;
   vInd.reserve(vertices_.size());
   for(size_t i=0; i<vertices_.size(); ++i)
      vInd.push_back(i);

   // the index vector is sorted according to the vertex positions
   sort( vInd.begin(), vInd.end(), LessThan( vertices_, tolerance ) );


   // now remove duplicate vertices. This is done using an algorithm similar
   // to std::unique(). However an additional map is built, that maps the deleted indices
   // to the equivalent index
   // this map is needed to translate the indices_ vector that stores the face information
   map<size_t,size_t> oldToNewIndex;
   size_t removedVertices = 0;
   vector<size_t>::iterator lastValid = vInd.begin();
   for(vector<size_t>::iterator i = vInd.begin()+1; i != vInd.end(); ++i)
   {
      bool verticesEqual = true;
      for(size_t d=0; d<3; ++d)
         if( !realIsEqual(vertices_[*i][d], vertices_[*lastValid][d], tolerance) ){
             verticesEqual = false;
             break;
         }

      if( verticesEqual )
      {
         //delete vertex at position i
         oldToNewIndex[*i] = numeric_cast< size_t >( lastValid - vInd.begin() );
         removedVertices++;
      }
      else{
         ++lastValid;
         *lastValid = *i;
      }
   }
   ++lastValid;


   vector<vertex_t> newVertices;
   newVertices.reserve( vertices_.size() - removedVertices );

   for(vector<size_t>::iterator i = vInd.begin(); i != lastValid; ++i)
   {
      oldToNewIndex[*i] = newVertices.size();
      newVertices.push_back( vertices_[*i] );
   }

   // exchange old by new vertices
   vertices_.swap( newVertices );


   if( hasVertexNormals() ){
      vector<vertex_t> newVertexNormals;
      newVertexNormals.reserve( vertexNormals_.size() - removedVertices );

      for(vector<size_t>::iterator i = vInd.begin(); i != lastValid; ++i)
         newVertexNormals.push_back( vertexNormals_[*i] );

      // exchange old by new vertices
      vertexNormals_.swap( newVertexNormals );
   }


   // adapt the indices in the triangles
   for(size_t i=0; i < vertexIndices_.size(); ++i) {
      WALBERLA_ASSERT_LESS( oldToNewIndex[ vertexIndices_[i] ], vInd.size() - removedVertices);
      vertexIndices_[i] = index_t( oldToNewIndex[ vertexIndices_[i] ] );
   }

   return removedVertices;
}



/// Used only internally in function "TriangleMesh::split( vector<TriangleMesh>& meshes )"
struct TriangleMeshNode{
   TriangleMesh::index_t vOld;
   bool used;

   set   < TriangleMeshNode*   > conns;
   vector< TriangleMesh::index_t > nOld;
};

//**********************************************************************************************************************
/*! \brief Split mesh into unconnected meshes.
 *
 * All triangle connections are checked. Unconnected meshes are split up.
 *
 * \param meshes empty vector where all unconnected meshes are returned
 */
//*********************************************************************************************************************/
void TriangleMesh::split( vector<TriangleMesh>& meshes ) const
{
   // build up connection graph
   map   < index_t, TriangleMeshNode* > nodes;
   vector<          TriangleMeshNode* > tnode( 3u );

   for( size_t triangle = size_t(0u); triangle < getNumTriangles(); ++triangle )
   {
      for( uint8_t index = uint8_t(0u); index < 3; ++index )
      {
         const index_t vIndex = getVertexIndex( triangle, index );
         if( nodes.find( vIndex ) == nodes.end() ){
            TriangleMeshNode* node = new TriangleMeshNode();
            node->vOld = vIndex;
            node->used = false;
            nodes[vIndex] = node;
         }
         tnode[index] = nodes[vIndex];
         if( hasNormalIndices() )
         {
            const index_t nIndex = getNormalIndex( triangle, index );
            tnode[index]->nOld.push_back( nIndex );
         }
      }
      tnode[0]->conns.insert( tnode[1] );
      tnode[0]->conns.insert( tnode[2] );
      tnode[1]->conns.insert( tnode[2] );
      tnode[1]->conns.insert( tnode[0] );
      tnode[2]->conns.insert( tnode[0] );
      tnode[2]->conns.insert( tnode[1] );
   }

   // split vertices by trinagle connetions
   set< vector< TriangleMeshNode* > > ssnode;
   for( auto nit = nodes.begin(); nit != nodes.end(); ++nit )
   {
      TriangleMeshNode* node = nit->second;
      if( node->used )
         continue;

      vector< TriangleMeshNode* > snode;
      snode.push_back( node );
      node->used = true;

      for( size_t index = size_t(0u); index < snode.size(); ++index )
      {
         for( auto cit = snode[index]->conns.begin(); cit != snode[index]->conns.end(); ++cit )
         {
            TriangleMeshNode* child = *cit;
            if( child->used )
               continue;
            child->used = true;
            snode.push_back( child );
         }
      }
      ssnode.insert( snode );
   }

   // add vertices, colors and normals to new meshes
   meshes.resize( ssnode.size() );
   vector< map< index_t, index_t > > vid( ssnode.size() );
   vector< map< index_t, index_t > > nid( ssnode.size() );
   map< index_t, size_t > ind;
   size_t index = size_t(0u);
   for( auto srcIt = ssnode.begin(); srcIt != ssnode.end(); ++srcIt, ++index)
   {
      for( auto it = srcIt->begin(); it != srcIt->end(); ++it )
      {
         TriangleMeshNode* node = *it;
         index_t vIndex;
         if( hasVertexColors() )
            vIndex = meshes[index].addVertex( vertices_[node->vOld], vertexColors_[node->vOld] );
         else
            vIndex = meshes[index].addVertex( vertices_[node->vOld] );
         vid[index][node->vOld] = vIndex;

         if( hasNormalIndices() ){
            for( auto nOld = node->nOld.begin(); nOld != node->nOld.end(); ++nOld ){
               if( nid[index].find(*nOld) != nid[index].end() )
                  continue;
               const index_t nIndex = meshes[index].addVertexNormal( vertexNormals_[*nOld] );
               nid[index][*nOld] = nIndex;
            }
         } else if( hasVertexNormals() ){
            const index_t nIndex = meshes[index].addVertexNormal( vertexNormals_[node->vOld] );
            nid[index][node->vOld] = nIndex;
         }

         ind[node->vOld] = index;
      }
   }

   // add triangles to new meshes
   for( size_t triangle = size_t(0u); triangle < getNumTriangles(); ++triangle )
   {
      const index_t vIndex0 = getVertexIndex( triangle, uint8_t(0u) );
      const index_t vIndex1 = getVertexIndex( triangle, uint8_t(1u) );
      const index_t vIndex2 = getVertexIndex( triangle, uint8_t(2u) );

      index = ind[vIndex0];

      if( hasNormalIndices() )
      {
         const index_t nIndex0 = getNormalIndex( triangle, uint8_t(0u) );
         const index_t nIndex1 = getNormalIndex( triangle, uint8_t(1u) );
         const index_t nIndex2 = getNormalIndex( triangle, uint8_t(2u) );

         meshes[index].addTriangle(
            vid[index][vIndex0], vid[index][vIndex1], vid[index][vIndex2],
            nid[index][nIndex0], nid[index][nIndex1], nid[index][nIndex2] );
      }
      else
      {
         meshes[index].addTriangle( vid[index][vIndex0], vid[index][vIndex1], vid[index][vIndex2] );
      }
   }

   // clear memory
   for( auto it = nodes.begin(); it != nodes.end(); ++it )
   {
      delete it->second;
   }
}



//**********************************************************************************************************************
/*! \brief Merges a second mesh into the given mesh.
 *
 * All vertices and faces of the second mesh are added, no vertex duplicate checking is done.
 *
 * \param other  mesh that has to merge into this
 * \param offset before adding the vertices, they are moved by the given offset
 */
//*********************************************************************************************************************/
void TriangleMesh::merge(const TriangleMesh & other, const Vector3<real_t> & offset)
{
   index_t oldNumVertices      = index_c( vertices_.size() );
   index_t oldNumVertexNormals = index_c( vertexNormals_.size() );

   // Add vertices
   for(index_t i=0; i < other.getNumVertices(); ++i) {
      vertex_t v = other.getVertex(i);
      v += offset;
      addVertex(v);
   }

   // Add Normals
   std::copy(other.vertexNormals_.begin(), other.vertexNormals_.end(), std::back_inserter(vertexNormals_) );
   std::copy(other.vertexColors_.begin(),  other.vertexColors_.end(),  std::back_inserter(vertexColors_)  );

   // Add faces
   for( auto it = other.vertexIndices_.begin(); it != other.vertexIndices_.end(); ++it )
      vertexIndices_.push_back( index_c( *it + oldNumVertices ) );

   for( auto it = other.normalIndices_.begin(); it != other.normalIndices_.end(); ++it )
      normalIndices_.push_back( index_c( *it + oldNumVertexNormals ) );
}

math::AABB TriangleMesh::getAABB() const
{

   if( vertices_.empty() )
      WALBERLA_ABORT( "You are trying to compute the bounding box of an empty mesh!" );

   return math::AABB( vertices_.begin(), vertices_.end() );
}

void TriangleMesh::translate( const Vector3<real_t> & offset )
{
   for( auto it = vertices_.begin(); it != vertices_.end(); ++it )
      *it += offset;
}

void TriangleMesh::scale( const Vector3<real_t> & scaleFactors )
{
   for( auto it = vertices_.begin(); it != vertices_.end(); ++it )
   {
      (*it)[0] *= scaleFactors[0];
      (*it)[1] *= scaleFactors[1];
      (*it)[2] *= scaleFactors[2];
   }
}

void TriangleMesh::exchangeAxes( uint_t xAxisId, uint_t yAxisId, uint_t zAxisId )
{
   for( auto it = vertices_.begin(); it != vertices_.end(); ++it )
   {
      vertex_t copy = *it;
      (*it)[0] = copy[xAxisId];
      (*it)[1] = copy[yAxisId];
      (*it)[2] = copy[zAxisId];
   }
}


real_t TriangleMesh::volume() const
{
   vertex_t v0;
   vertex_t v1;
   vertex_t v2;
   real_t result(0);

   for(size_t i = 0; i < getNumTriangles(); ++i)
   {
      getTriangle( i, v0, v1, v2 );
      result += ( v0 * ( v1 % v2 ) ) / real_t(6);
   }

   return std::fabs(result);
}

real_t TriangleMesh::surfaceArea() const
{
   vertex_t v0;
   vertex_t v1;
   vertex_t v2;
   real_t result(0);

   for(size_t i = 0; i < getNumTriangles(); ++i)
   {
      getTriangle( i, v0, v1, v2 );
      result += ( ( v1 - v0 ) % ( v2 - v0 ) ).length();
   }

   return result * real_t( 0.5 );
}


} // namespace geometry
} // namespace walberla


