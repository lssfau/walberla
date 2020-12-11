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
//! \file ConvexPolyhedron.cpp
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "ConvexPolyhedron.h"

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "core/math/Matrix3.h"
#include "core/debug/Debug.h"

#include "mesh_common/MeshOperations.h"
#include "mesh_common/QHull.h"

#include "pe/Materials.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <cmath>

namespace walberla {
namespace mesh {
namespace pe {

//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================


//*************************************************************************************************
/*!\brief Constructor for the ConvexPolyhedron class.
*
* \param sid Unique system-specific ID for the ConvexPolyhedron.
* \param uid User-specific ID for the ConvexPolyhedron.
* \param gpos Global geometric center of the ConvexPolyhedron.
* \param rpos The relative position within the body frame of a superordinate body.
* \param q The orientation of the ConvexPolyhedron's body frame in the global world frame.
* \param radius The radius of the ConvexPolyhedron \f$ (0..\infty) \f$.
* \param material The material of the ConvexPolyhedron.
* \param global specifies if the ConvexPolyhedron should be created in the global storage
* \param communicating specifies if the ConvexPolyhedron should take part in synchronization (syncNextNeighbour, syncShadowOwner)
* \param infiniteMass specifies if the ConvexPolyhedron has infinite mass and will be treated as an obstacle
*/
ConvexPolyhedron::ConvexPolyhedron( id_t sid, id_t uid, const Vec3& gpos, const Quat& q,
                                    const TriangleMesh & mesh, MaterialID material,
                                    const bool global, const bool communicating, const bool infiniteMass )
   : GeomPrimitive( getStaticTypeID(), sid, uid, material ),  // Initialization of the parent class
     mesh_( mesh )
{
   init( gpos, q, global, communicating, infiniteMass );
}
//*************************************************************************************************


void ConvexPolyhedron::init( const Vec3& gpos,  const Quat& q,
                             const bool global, const bool communicating, const bool infiniteMass )
{
   WALBERLA_ASSERT_FLOAT_EQUAL( (toWalberla( computeCentroid( mesh_ ) ) - Vec3() ).length(), real_t(0) );
   WALBERLA_ASSERT_GREATER_EQUAL( mesh_.n_vertices(), 4 );
   
   mesh_.request_face_normals();
   mesh_.update_face_normals();

   // Calculate the bounding sphere radius first, as setPosition will trigger a call to calcBoundingBox(), which needs it
   real_t maxSqRadius(0);
   for(auto vh : mesh_.vertices())
   {
      real_t sqRadius = mesh_.point( vh ).sqrnorm();
      if( sqRadius > maxSqRadius )
         maxSqRadius = sqRadius;
   }
   boundingSphereRadius_ = std::sqrt( maxSqRadius );

   // Setting the center of the ConvexPolyhedron
   setPosition(gpos);

   // Initializing the instantiated ConvexPolyhedron
   setOrientation(q);
   setGlobal( global );
   if (infiniteMass)
   {
      setMassAndInertiaToInfinity();
   } else
   {
      // sets inverse mass and interatio tensor
      setMassAndInertia( getVolume() * Material::getDensity( getMaterial() ), mesh::computeInertiaTensor( mesh_ ) );
   }
   setCommunicating( communicating );
   setFinite( true );

   // Setting the axis-aligned bounding box

   ConvexPolyhedron::calcBoundingBox();

   octandVertices_[0] = supportVertex( TriangleMesh::Normal( real_t( 1), real_t( 1), real_t( 1) ), *mesh_.vertices_begin() );
   octandVertices_[1] = supportVertex( TriangleMesh::Normal( real_t( 1), real_t( 1), real_t(-1) ), *mesh_.vertices_begin() );
   octandVertices_[2] = supportVertex( TriangleMesh::Normal( real_t( 1), real_t(-1), real_t( 1) ), *mesh_.vertices_begin() );
   octandVertices_[3] = supportVertex( TriangleMesh::Normal( real_t( 1), real_t(-1), real_t(-1) ), *mesh_.vertices_begin() );
   octandVertices_[4] = supportVertex( TriangleMesh::Normal( real_t(-1), real_t( 1), real_t( 1) ), *mesh_.vertices_begin() );
   octandVertices_[5] = supportVertex( TriangleMesh::Normal( real_t(-1), real_t( 1), real_t(-1) ), *mesh_.vertices_begin() );
   octandVertices_[6] = supportVertex( TriangleMesh::Normal( real_t(-1), real_t(-1), real_t( 1) ), *mesh_.vertices_begin() );
   octandVertices_[7] = supportVertex( TriangleMesh::Normal( real_t(-1), real_t(-1), real_t(-1) ), *mesh_.vertices_begin() );
}



//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the ConvexPolyhedron class.
 */
ConvexPolyhedron::~ConvexPolyhedron()
{
   // Logging the destruction of the ConvexPolyhedron
   WALBERLA_LOG_DETAIL( "Destroyed ConvexPolyhedron " << sid_ );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Calculation of the bounding box of the ConvexPolyhedron.
 *
 * \return void
 *
 * This function updates the axis-aligned bounding box of the ConvexPolyhedron primitive according to the
 * current position and orientation of the ConvexPolyhedron. Note that the bounding box is increased in
 * all dimensions by pe::contactThreshold to guarantee that rigid bodies in close proximity of
 * the ConvexPolyhedron are also considered during the collision detection process.
 */
void ConvexPolyhedron::calcBoundingBox()
{
   const Vec3 & pos = getPosition();
   const real_t r = boundingSphereRadius_ + contactThreshold;
   aabb_.initMinMaxCorner( pos[0] - r, pos[1] - r , pos[2] - r,
                           pos[0] + r, pos[1] + r , pos[2] + r  );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the volume of the ConvexPolyhedron.
*
* \return void
*/
real_t ConvexPolyhedron::getVolume() const
{
   return mesh::computeVolume( mesh_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculation of the surface area of the ConvexPolyhedron.
*
* \return void
*/
real_t ConvexPolyhedron::getSurfaceArea() const
{
   return mesh::computeSurfaceArea( mesh_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimates the point which is farthest in direction \a d.
 *
 * \param d The normalized search direction in world-frame coordinates.
 * \return The support point in world-frame coordinates in direction a\ d.
 */
Vec3 ConvexPolyhedron::support( const Vec3& d ) const
{
   if (math::equal(d.length(), real_t(0))) return Vec3(0,0,0);

   TriangleMesh::Normal d_loc = toOpenMesh( vectorFromWFtoBF(d) );
   
   TriangleMesh::VertexHandle startVertex;
   if(d_loc[0] >= real_t( 0 ))
   {
      if(d_loc[1] >= real_t( 0 ))
      {
         startVertex = d_loc[2] >= real_t( 0 ) ? octandVertices_[0] : octandVertices_[1];
      }
      else // d_loc[1] < 0
      {
         startVertex = d_loc[2] >= real_t( 0 ) ? octandVertices_[2] : octandVertices_[3];
      }
   }
   else // d_loc[0] < 0
   {
      if(d_loc[1] >= real_t( 0 ))
      {
         startVertex = d_loc[2] >= real_t( 0 ) ? octandVertices_[4] : octandVertices_[5];
      }
      else // d_loc[1] < 0
      {
         startVertex = d_loc[2] >= real_t( 0 ) ? octandVertices_[6] : octandVertices_[7];
      }
   }

   TriangleMesh::VertexHandle vh = supportVertex( d_loc, startVertex );

   return pointFromBFtoWF( toWalberla( mesh_.point( vh ) ) );

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimates the vertex which is farthest in direction \a d.
*
* \param d The normalized search direction in body-frame coordinates.
* \return The support vertex in direction a\ d.
*/
TriangleMesh::VertexHandle ConvexPolyhedron::supportVertex( const TriangleMesh::Normal & d, const TriangleMesh::VertexHandle startVertex ) const
{
   TriangleMesh::VertexHandle maxScalarProductVertex = startVertex;
   real_t maxScalarProduct = mesh_.point(maxScalarProductVertex) | d;

   bool isExtremum = false;
   while( !isExtremum )
   {
      isExtremum = true;
      for(auto vh : mesh_.vv_range( maxScalarProductVertex ))
      {
         const real_t sp = mesh_.point(vh) | d;
         if(sp > maxScalarProduct)
         {
            isExtremum = false;
            maxScalarProductVertex = vh;
            maxScalarProduct = sp;
            break;
         }
      }
   }

   return maxScalarProductVertex;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimates the point which is farthest in direction \a d.
 *
 * \param d The normalized search direction in world-frame coordinates
 * \return The support point in world-frame coordinates in direction a\ d extended by a vector in
 *         direction \a d of length \a pe::contactThreshold.
 */
Vec3 ConvexPolyhedron::supportContactThreshold( const Vec3& /*d*/ ) const
{
 WALBERLA_ABORT("supportContactThreshold not implemented!");
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies inside the ConvexPolyhedron.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies inside the ConvexPolyhedron, \a false if not.
 */
bool ConvexPolyhedron::containsRelPointImpl( real_t px, real_t py, real_t pz ) const
{
   if( px * px + py * py + pz * pz > boundingSphereRadius_ * boundingSphereRadius_ )
      return false;

   return std::none_of(mesh_.faces().begin(),
                       mesh_.faces().end(),
                       [&](auto fh)
                       {
                          //check if point is on positive side of the face
                          const TriangleMesh::Normal &n = mesh_.normal(fh); // Plane normal
                          const TriangleMesh::Point &pp = mesh_.point(
                                mesh_.to_vertex_handle(mesh_.halfedge_handle(fh))); // Point on plane

                          return (n[0] * (px - pp[0]) + n[1] * (py - pp[1]) + n[2] * (pz - pp[2]) >= real_t(0));
                       });
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks, whether a point in body relative coordinates lies on the surface of the ConvexPolyhedron.
 *
 * \param px The x-component of the relative coordinate.
 * \param py The y-component of the relative coordinate.
 * \param pz The z-component of the relative coordinate.
 * \return \a true if the point lies on the surface of the ConvexPolyhedron, \a false if not.
 *
 * The tolerance level of the check is pe::surfaceThreshold.
 */
bool ConvexPolyhedron::isSurfaceRelPointImpl( real_t /*px*/, real_t /*py*/, real_t /*pz*/ ) const
{
 WALBERLA_ABORT("isSurfaceRelPointImpl not implemented!");
}
//*************************************************************************************************




//=================================================================================================
//
//  OUTPUT FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Output of the current state of a ConvexPolyhedron.
 *
 * \param os Reference to the output stream.
 * \param tab Indentation in front of every line of the ConvexPolyhedron output.
 * \return void
 */
void ConvexPolyhedron::print( std::ostream& os, const char* tab ) const
{
   using std::setw;

   os << tab << " ConvexPolyhedron " << uid_ << "\n"
      << tab << "   #Vertices = " << mesh_.n_vertices() << "\n"
      << tab << "   #Faces    = " << mesh_.n_faces()    << "\n";


   os << tab << "   Fixed: " << isFixed() << " , sleeping: " << !isAwake() << "\n";

   os << tab << "   System ID         = " << getSystemID() << "\n"
      << tab << "   Total mass        = " << getMass() << "\n"
      << tab << "   Material          = " << Material::getName( material_ ) << "\n"
      << tab << "   Owner             = " << MPITrait.getOwner() << "\n"
      << tab << "   Global position   = " << getPosition() << "\n"
      << tab << "   Relative position = " << getRelPosition() << "\n"
      << tab << "   Linear velocity   = " << getLinearVel() << "\n"
      << tab << "   Angular velocity  = " << getAngularVel() << "\n";

   //   if( verboseMode )
   //   {
   Mat3 R = getRotation();
   os << tab << "   Bounding box      = " << getAABB() << "\n"
      << tab << "   Quaternion        = " << getQuaternion() << "\n"
      << tab << "   Rotation matrix   = ( " << setw(9) << R[0] << " , " << setw(9) << R[1] << " , " << setw(9) << R[2] << " )\n"
      << tab << "                       ( " << setw(9) << R[3] << " , " << setw(9) << R[4] << " , " << setw(9) << R[5] << " )\n"
      << tab << "                       ( " << setw(9) << R[6] << " , " << setw(9) << R[7] << " , " << setw(9) << R[8] << " )\n";

   os << std::setiosflags(std::ios::right)
      << tab << "   Moment of inertia = ( " << setw(9) << I_[0] << " , " << setw(9) << I_[1] << " , " << setw(9) << I_[2] << " )\n"
      << tab << "                       ( " << setw(9) << I_[3] << " , " << setw(9) << I_[4] << " , " << setw(9) << I_[5] << " )\n"
      << tab << "                       ( " << setw(9) << I_[6] << " , " << setw(9) << I_[7] << " , " << setw(9) << I_[8] << " )\n"
      << std::resetiosflags(std::ios::right);
   //   }
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global output operator for ConvexPolyhedra.
 *
 * \param os Reference to the output stream.
 * \param s Reference to a constant ConvexPolyhedron object.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, const ConvexPolyhedron& s )
{
   os << "--" << "ConvexPolyhedron PARAMETERS"
      << "-------------------------------------------------------------\n";
   s.print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for ConvexPolyhedron handles.
 *
 * \param os Reference to the output stream.
 * \param s Constant ConvexPolyhedron handle.
 * \return Reference to the output stream.
 */
std::ostream& operator<<( std::ostream& os, ConstConvexPolyhedronID s )
{
   os << "--" << "ConvexPolyhedron PARAMETERS"
      << "-------------------------------------------------------------\n";
   s->print( os, "" );
   os << "--------------------------------------------------------------------------------\n"
      << std::endl;
   return os;
}
//*************************************************************************************************

id_t ConvexPolyhedron::staticTypeID_ = std::numeric_limits<id_t>::max();

} // namespace pe
} // namespace mesh
} // namespace walberla
