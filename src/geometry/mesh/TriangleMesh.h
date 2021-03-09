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
//! \file TriangleMesh.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Class for storing Triangle Mesh using vertex index array
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/math/AABBFwd.h"
#include "core/math/Vector3.h"

#include <cassert>
#include <vector>


namespace walberla {

namespace geometry {



   //*******************************************************************************************************************
   /*! \brief Class for storing a triangle mesh.
    *
    *  Format:
    *    - vertices are stored in an array
    *    - triangles are stored as an index triple into the vertex array, three consecutive
    *      entries in the indices_ array define a triangle
    *
    *  When adding new vertices no check is done if this vertex already exist.
    *  To remove duplicate vertices (and adapt the indices that store the triangles) use removeDuplicateVertices()
    *
    *  \ingroup postprocessing
    */
   //******************************************************************************************************************/
   class TriangleMesh
   {
   public:
      using index_t = uint32_t;
      template<typename T> static index_t index_c( T x ) { return numeric_cast<index_t>(x); }

      using vertex_t = Vector3<real_t>;
      using normal_t = Vector3<real_t>;
      using color_t = Vector3<float>;

      //** Size Information  *******************************************************************************************
      /*! \name Size Information */
      //@{
      inline size_t   getNumTriangles()     const { return vertexIndices_.size() / 3;        }
      inline size_t   getNumVertexIndices() const { return vertexIndices_.size();            }
      inline size_t   getNumNormalIndices() const { return normalIndices_.size();            }
      inline index_t  getNumVertices()      const { return index_c( vertices_.size() );      }
      inline index_t  getNumNormals()       const { return index_c( vertexNormals_.size() ); }
      math::AABB      getAABB()             const;
      //@}
      //****************************************************************************************************************



      //** Data Access     *********************************************************************************************
      /*! \name Data Access */
      //@{
      inline const vertex_t & getVertex(index_t i)                          const { return vertices_[i];      }
      inline const normal_t & getVertexNormal(index_t i)                    const { return vertexNormals_[i]; }
      inline const color_t  & getVertexColor(index_t i)                     const { return vertexColors_[i]; }
      inline void             clear();
      inline index_t          getVertexIndex(index_t i)                     const;
      inline index_t          getVertexIndex(size_t triangle,uint8_t index) const;
      inline index_t          getNormalIndex(index_t i)                     const;
      inline index_t          getNormalIndex(size_t triangle,uint8_t index) const;

      inline bool             hasVertexNormals()                            const { return !vertexNormals_.empty(); }
      inline bool             hasVertexColors()                             const { return !vertexColors_.empty();  }
      inline bool             hasNormalIndices()                            const { return !normalIndices_.empty(); }

      inline void             getTriangle( size_t triangleIdx, vertex_t & v0, vertex_t & v1, vertex_t & v2 ) const;
      inline void             getTriangle( size_t triangleIdx, vertex_t & v0, vertex_t & v1, vertex_t & v2,
                                                                color_t & c0,  color_t & c1,  color_t & c2) const;
      inline void             getTriangleVertexNormals( size_t triangleIdx,
                                                        normal_t & n0, normal_t & n1, normal_t & n2 ) const;

      const std::vector<index_t>  & getVertexIndices() const { return vertexIndices_; }
      const std::vector<index_t>  & getNormalIndices() const { return normalIndices_; }
      const std::vector<vertex_t> & getVertices()      const { return vertices_;      }
      const std::vector<normal_t> & getVertexNormals() const { return vertexNormals_; }
      const std::vector<color_t>  & getVertexColors()  const { return vertexColors_;  }

            std::vector<index_t>  & getVertexIndices()       { return vertexIndices_; }
            std::vector<index_t>  & getNormalIndices()       { return normalIndices_; }
            std::vector<vertex_t> & getVertices()            { return vertices_;      }
            std::vector<normal_t> & getVertexNormals()       { return vertexNormals_; }
            std::vector<color_t>  & getVertexColors()        { return vertexColors_;  }

      template< typename OutputIterator >
      void getVerticesOfColor( const color_t & color, OutputIterator outIt ) const;

      real_t volume() const;
      real_t surfaceArea() const;
      //@}
      //****************************************************************************************************************


      //** Inserting  **************************************************************************************************
      /*! \name Inserting */
      //@{
      index_t addVertex(const vertex_t & v);
      index_t addVertex(const vertex_t & v, const color_t & c);
      index_t addVertexNormal(const normal_t & n);
      void    addTriangle(index_t v0Idx, index_t v1Idx, index_t v2Idx);
      void    addTriangle(index_t v0Idx, index_t v1Idx, index_t v2Idx, index_t n0Idx, index_t n1Idx, index_t n2Idx);
      //@}
      //****************************************************************************************************************


      //** Mesh Operations *********************************************************************************************
      /*! \name Mesh Operations */
      //@{
      void   merge(const TriangleMesh & other, const Vector3<real_t> & offset = Vector3<real_t>(0.0) );
      size_t removeDuplicateVertices( real_t tolerance = real_t(1e-4) );
      void   translate( const Vector3<real_t> & offset );
      void   scale( real_t scaleFactor ) { scale( Vector3<real_t>( scaleFactor, scaleFactor, scaleFactor ) ); }
      void   scale( const Vector3<real_t> & scaleFactors );
      void   split( std::vector<TriangleMesh>& meshes ) const;
      void   exchangeAxes( uint_t xAxisId, uint_t yAxisId, uint_t zAxisId );
      //@}
      //****************************************************************************************************************


   protected:
      std::vector<index_t>  vertexIndices_;
      std::vector<index_t>  normalIndices_;
      std::vector<vertex_t> vertices_;
      std::vector<normal_t> vertexNormals_;
      std::vector<color_t>  vertexColors_;

   };

   inline TriangleMesh::index_t TriangleMesh::getVertexIndex(index_t i) const
   {
      WALBERLA_ASSERT_LESS( i, vertexIndices_.size() );
      return vertexIndices_[i];
   }

   inline TriangleMesh::index_t TriangleMesh::getVertexIndex(size_t triangle,uint8_t index) const
   {
      WALBERLA_ASSERT_LESS( triangle*3+index, vertexIndices_.size() );
      return vertexIndices_[triangle*3+index];
   }

   inline TriangleMesh::index_t TriangleMesh::getNormalIndex(index_t i) const
   {
      WALBERLA_ASSERT_LESS( i, normalIndices_.size() );
      return normalIndices_[i];
   }

   inline TriangleMesh::index_t TriangleMesh::getNormalIndex(size_t triangle,uint8_t index) const
   {
      WALBERLA_ASSERT_LESS( triangle*3+index, normalIndices_.size() );
      return normalIndices_[triangle*3+index];
   }

   inline void TriangleMesh::getTriangle( size_t triangleIdx, vertex_t & v0, vertex_t & v1, vertex_t & v2 ) const
   {
      v0 = vertices_[ getVertexIndex( triangleIdx, 0 ) ];
      v1 = vertices_[ getVertexIndex( triangleIdx, 1 ) ];
      v2 = vertices_[ getVertexIndex( triangleIdx, 2 ) ];
   }

   inline void TriangleMesh::getTriangle( size_t triangleIdx, vertex_t & v0, vertex_t & v1, vertex_t & v2,
                                                        color_t & c0,  color_t & c1,  color_t & c2 ) const
   {
      v0 = vertices_[ getVertexIndex( triangleIdx, 0 ) ];
      v1 = vertices_[ getVertexIndex( triangleIdx, 1 ) ];
      v2 = vertices_[ getVertexIndex( triangleIdx, 2 ) ];
      c0 = vertexColors_[ getVertexIndex( triangleIdx, 0 ) ];
      c1 = vertexColors_[ getVertexIndex( triangleIdx, 1 ) ];
      c2 = vertexColors_[ getVertexIndex( triangleIdx, 2 ) ];
   }

   inline void TriangleMesh::getTriangleVertexNormals( size_t triangleIdx, normal_t & n0, normal_t & n1, normal_t & n2 ) const
   {
      n0 = vertexNormals_[ getNormalIndex( triangleIdx, 0 ) ];
      n1 = vertexNormals_[ getNormalIndex( triangleIdx, 1 ) ];
      n2 = vertexNormals_[ getNormalIndex( triangleIdx, 2 ) ];
   }

   inline void TriangleMesh::clear()
   {
      vertexIndices_.clear();
      normalIndices_.clear();
      vertices_.clear();
      vertexNormals_.clear();
      vertexColors_.clear();
   }

   template< typename OutputIterator >
   void TriangleMesh::getVerticesOfColor( const color_t & color, OutputIterator outIt ) const
   {
      WALBERLA_ASSERT( hasVertexColors() );
      WALBERLA_ASSERT_EQUAL( vertexColors_.size(), vertices_.size() );

      auto vertexIt = vertices_.begin();
      auto colorIt  = vertexColors_.begin();
      
      while( vertexIt != vertices_.end() )
      {
         if( *colorIt++ == color )
            *outIt++ = *vertexIt;
         
         ++vertexIt;
      }
   }

} // namespace geometry
} // namespace walberla


