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
//! \file EPA.h
//! \author Tobias Scharpff
//! \author Tobias Leemann
//
//  DISCLAIMER: The following source file contains modified code from the SOLID-3.5 library for
//  interference detection as it is published in the book "Collision Detection in Interactive
//  3D Environments" by Gino van den Bergen <info@dtecta.com>. Even though the original source
//  was published under the GPL version 2 not allowing later versions, the original author of the
//  source code permitted the relicensing of the SOLID-3.5 library under the GPL version 3 license.
//
//=================================================================================================

#pragma once

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "GJK.h"
#include <pe/Thresholds.h>
#include <pe/Types.h>

#include <core/math/Constants.h>
#include <core/math/Limits.h>
#include <core/math/Matrix3.h>
#include <core/math/Quaternion.h>

#include <vector>

namespace walberla {
namespace pe {
namespace fcd {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The Expanding-Polytope Algorithm.
 * \ingroup fine_collision_detection
 */
class EPA
{
private :
   //**Type definitions****************************************************************************
   class EPA_Edge;
   class EPA_Triangle;
   class EPA_TriangleComp;

   using EPA_EntryBuffer = std::vector<EPA_Triangle>;
   using EPA_EntryHeap = std::vector<EPA_Triangle *>;
   using EPA_EdgeBuffer = std::vector<EPA_Edge>;
   //**********************************************************************************************

public:
   //**Query functions*****************************************************************************
   /*!\name Query functions */
   //@{
   bool doEPAcontactThreshold( GeomPrimitive &geom1, GeomPrimitive &geom2, const GJK& gjk, Vec3& normal,
                                      Vec3& contactPoint, real_t& penetrationDepth);


   bool doEPAcontactThreshold( GeomPrimitive &geom1, GeomPrimitive &geom2, const GJK& gjk, Vec3& normal,
                                      Vec3& contactPoint, real_t& penetrationDepth, real_t eps_rel);
   
   bool doEPAmargin( GeomPrimitive &geom1, GeomPrimitive &geom2, const GJK& gjk, Vec3& normal,
                            Vec3& contactPoint, real_t& penetrationDepth, real_t margin);

   bool doEPA( GeomPrimitive &geom1, GeomPrimitive &geom2, const GJK& gjk, Vec3& normal,
                      Vec3& contactPoint, real_t& penetrationDepth, real_t margin, real_t eps_rel );



   //@}
   //**********************************************************************************************

   //**Getter/Setter functions*****************************************************************************
   /*!\name Getter and Setter functions */
   //@{
   inline void setMaxSupportPoints( size_t maxSupportPoints) {maxSupportPoints_ = maxSupportPoints;}

   inline size_t getMaxSupportPoints() {return maxSupportPoints_;}

   inline void setMaxTriangles( size_t maxTriangles) {maxTriangles_ = maxTriangles;}

   inline size_t getMaxTriangles() {return maxTriangles_;}

   inline int getNumIterations() const {return numIterations_; }

   inline bool useSphereOptimization() const {return bUseSphereOptimization_; }
   inline void useSphereOptimization(const bool useIt) {bUseSphereOptimization_ = useIt;}

   //@}
   //**********************************************************************************************

private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline void pushSupportMargin(const GeomPrimitive &geom1, const GeomPrimitive &geom2, const Vec3& dir, const real_t margin,
                                 std::vector<Vec3>& epaVolume, std::vector<Vec3>& supportA, std::vector<Vec3>& supportB);

   inline void replaceSupportMargin(const GeomPrimitive &geom1, const GeomPrimitive &geom2, const Vec3& dir, const real_t margin,
                                    std::vector<Vec3>& epaVolume, std::vector<Vec3>& supportA, std::vector<Vec3>& supportB, size_t indexToReplace);

   inline void removeSupportMargin(std::vector<Vec3>& epaVolume, std::vector<Vec3>& supportA, std::vector<Vec3>& supportB);

   inline bool originInTetrahedron            ( const Vec3& A, const Vec3& B, const Vec3& C,
                                                const Vec3& D );
   inline bool originInTetrahedronVolumeMethod( const Vec3& A, const Vec3& B, const Vec3& C,
                                                const Vec3& D );
   inline bool pointInTetrahedron             ( const Vec3& A, const Vec3& B, const Vec3& C,
                                                const Vec3& D, const Vec3& point );
   bool searchTetrahedron              (GeomPrimitive &geom1, GeomPrimitive &geom2, std::vector<Vec3>& epaVolume,
                                               std::vector<Vec3>& supportA, std::vector<Vec3>& supportB, EPA_EntryBuffer& entryBuffer, real_t margin );

   void createInitialTetrahedron       ( size_t top, size_t frontLeft, size_t frontRight,
                                                size_t back, std::vector<Vec3>& epaVolume,
                                                EPA_EntryBuffer& entryBuffer );

   void createInitialSimplex           ( size_t numPoints, GeomPrimitive &geom1, GeomPrimitive &geom2,
                                                std::vector<Vec3>& supportA,
                                                std::vector<Vec3>& supportB,
                                                std::vector<Vec3>& epaVolume,
                                                EPA_EntryBuffer& entryBuffer, real_t margin );
   inline real_t calculateCircle              ( const Vec3& A, const Vec3& B, const Vec3& C,
                                                const Vec3& D, Vec3& center );
   //@}
   //**********************************************************************************************


private:
   //EPA constants
   size_t maxSupportPoints_ = 100;
   size_t maxTriangles_     = 200;

   int numIterations_ = 0;
   bool bUseSphereOptimization_ = false;
};
//*************************************************************************************************




//=================================================================================================
//
//  EPA::EPA_EDGE CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Class storing Information about an Edge of the EPA-Polytope
 */
class EPA::EPA_Edge {
public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   EPA_Edge( EPA_Triangle* triangle, size_t index );
   //@}
   //**********************************************************************************************

   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   EPA_Triangle* getTriangle() const;
   size_t        getIndex()    const;
   size_t        getStart()    const;
   size_t        getEnd()      const;
   //@}
   //**********************************************************************************************

private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   EPA_Triangle* triangle_; //!< the EPA triangle the edge is contained in
   size_t startIdx_; //!< the index of the point the edge starts at (0, 1, 2)
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  EPA::EPA_TRIANGLE CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Class storing Information about a triangular facette (Triangle) of the EPA-Polytope
 *
 * see Collision detection in interactive 3D environments; Gino van den bergen page 155
 */
class EPA::EPA_Triangle {
public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline EPA_Triangle( size_t a, size_t b, size_t c, const std::vector<Vec3>& points );
   //@}
   //**********************************************************************************************

   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline size_t      operator[]( size_t i )                    const;
   inline const Vec3& getClosest()                                     const;
   inline const Vec3& getNormal()                                      const;
   inline Vec3        getClosestPoint(const std::vector<Vec3>& points) const;
   inline real_t      getSqrDist()                                     const;
   inline bool        isObsolete()                                     const;
   inline bool        isClosestInternal()                              const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline bool        link( size_t edge0, EPA_Triangle* tria, size_t edge1 );
   inline void        silhouette( const Vec3& w, EPA_EdgeBuffer& edgeBuffer );
   //@}
   //**********************************************************************************************

private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void silhouette( size_t index, const Vec3& w, EPA_EdgeBuffer& edgeBuffer );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t         indices_[3];     //!< indices of the vertices of the triangle
   bool           obsolete_;       //!< flag to denote whether die triangle is visible from the new support point

   Vec3           closest_;        //!< the point closest to the origin of the affine hull of the triangle
   Vec3           normal_;         //!< normal pointing away from the origin
   real_t         bar_[3];         //!< the barycentric coordinate of closest_
   real_t         sqrDist_;        //!< =key; square distance of closest_ to the origin

   EPA_Triangle*  adjTriangle_[3]; //!< pointer to the triangle adjacent to edge i(=0,1,2)
   size_t         adjEdges_[3];    //!< for each adjoining triangle adjTriangle_[i], the index of the adjoining edge
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//=================================================================================================
//
//  EPA::EPA_TRIANGLECOMP CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief
 * Compare Triangles by their closest points to sort the triangle heap.
 */
class EPA::EPA_TriangleComp {
public:
   //**Binary function call operator***************************************************************
   /*!\name Binary function call operator */
   //@{
   inline bool operator()( const EPA_Triangle *tria1, const EPA_Triangle *tria2 );
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//=================================================================================================
//
//  EPA_EDGE CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief
 * Construct a new Triangle Edge.
 */
inline EPA::EPA_Edge::EPA_Edge( EPA_Triangle* triangle, size_t index )
   : triangle_(triangle)
   , startIdx_(index)
{
}
//*************************************************************************************************




//=================================================================================================
//
//  EPA_EDGE GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \brief Return the triangle this edge belongs to.
 */
inline EPA::EPA_Triangle* EPA::EPA_Edge::getTriangle() const
{
   return triangle_;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Get the Index of this edge in its triangle.
 */
inline size_t EPA::EPA_Edge::getIndex() const
{
   return startIdx_;
}
//*************************************************************************************************



//*************************************************************************************************
/*! \brief Return the start point index  of an edge.
 *
 */
inline size_t EPA::EPA_Edge::getStart() const
{
   return (*triangle_)[startIdx_];
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Return the end point index of an edge.
 */
inline size_t EPA::EPA_Edge::getEnd() const
{
   return (*triangle_)[(startIdx_+1) % 3];
}
//*************************************************************************************************




//=================================================================================================
//
//  EPA::EPA_TRIANGLE GET FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \brief Returns the index of the internal vertex i(=0,1,2) within the EPA scope.
 */
inline size_t EPA::EPA_Triangle::operator[]( size_t i ) const
{
   return indices_[i];
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Returns the point closest to the origin of the affine hull of the triangle, which is also the normal.
 */
inline const Vec3& EPA::EPA_Triangle::getClosest() const
{
   return closest_;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Returns the normal of the triangle. Normal is not normalized!
 */
inline const Vec3& EPA::EPA_Triangle::getNormal() const
{
   return normal_;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Calculates the corresponding closest point from the given points, using barycentric coordinates.
 */
inline Vec3 EPA::EPA_Triangle::getClosestPoint(const std::vector<Vec3>& points) const
{
   return  bar_[0] * points[indices_[0]]
         + bar_[1] * points[indices_[1]]
         + bar_[2] * points[indices_[2]];

}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Returns the squared distance to the closest to the origin of the affine hull of the triangle.
 */
inline real_t EPA::EPA_Triangle::getSqrDist() const
{
   return sqrDist_;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Returns true if the triangle is no longer part of the EPA polygon.
 */
inline bool EPA::EPA_Triangle::isObsolete() const
{
   return obsolete_;
}
//*************************************************************************************************


//*************************************************************************************************
/*! Returns true if the point closest to the origin of the affine hull of the triangle, lies inside the triangle.
 */
inline bool EPA::EPA_Triangle::isClosestInternal() const
{
   real_t tol = real_t(0.0);
   return bar_[0] >= tol
         && bar_[1] >= tol
         && bar_[2] >= tol;
}
//*************************************************************************************************



//=================================================================================================
//
//  EPA::EPA_TRIANGLECOMP BINARY FUNCTION CALL OPERATOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compare two triangles by their distance.
 */
inline bool EPA::EPA_TriangleComp::operator()( const EPA_Triangle *tria1,
                                               const EPA_Triangle *tria2 )
{
   return tria1->getSqrDist() > tria2->getSqrDist();
}
//*************************************************************************************************


//=================================================================================================
//
//  EPA UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \brief Calculates a support point of a body extended by threshold.
 * Adds this support and the base points at bodies A and B to the vector.
 * \param geom1 The body A.
 * \param geom2 The body B.
 * \param dir The support point direction.
 * \param margin Extension of the Body.
 */
inline void EPA::pushSupportMargin(const GeomPrimitive &geom1, const GeomPrimitive &geom2, const Vec3& dir, const real_t margin,
                                   std::vector<Vec3>& epaVolume, std::vector<Vec3>& supportA, std::vector<Vec3>& supportB)
{
   Vec3 ndir;
   if(floatIsEqual(dir.sqrLength(), real_t(1.0))){
      ndir = dir.getNormalized();
   }else{
      ndir = dir;
   }
   Vec3 sA = geom1.support(ndir);
   Vec3 sB = geom2.support(-ndir);
   supportA.push_back(sA);
   supportB.push_back(sB);

   Vec3 support = sA -sB + real_t(2.0) * ndir * margin;
   epaVolume.push_back(support);
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Calculates a support point of a body extended by threshold.
 * Replaces the old value in the vectors at "IndexToReplace" with this support and the base points at bodies A and B.
 * \param geom1 The body A.
 * \param geom2 The body B.
 * \param dir The support point direction.
 * \param margin Extension of the Body.
 */
inline void EPA::replaceSupportMargin(const GeomPrimitive &geom1, const GeomPrimitive &geom2, const Vec3& dir, const real_t margin,
                                      std::vector<Vec3>& epaVolume, std::vector<Vec3>& supportA, std::vector<Vec3>& supportB, size_t indexToReplace)
{
   Vec3 ndir;
   if(floatIsEqual(dir.sqrLength(), real_t(1.0))){
      ndir = dir.getNormalized();
   }else{
      ndir = dir;
   }
   Vec3 sA = geom1.support(ndir);
   Vec3 sB = geom2.support(-ndir);
   Vec3 support = sA -sB + real_t(2.0) * ndir * margin;

   supportA[indexToReplace] = sA;
   supportB[indexToReplace] = sB;
   epaVolume[indexToReplace] = support;
}
//*************************************************************************************************

//*************************************************************************************************
/*! \brief Removes a support point from the volume.
 */
inline void EPA::removeSupportMargin(std::vector<Vec3>& epaVolume, std::vector<Vec3>& supportA, std::vector<Vec3>& supportB)
{
   supportA.pop_back();
   supportB.pop_back();
   epaVolume.pop_back();
}
//*************************************************************************************************


//@}

} // namespace fcd
} // namespace pe
} // namespace walberla
