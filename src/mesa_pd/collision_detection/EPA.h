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

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/collision_detection/GJK.h>
#include <mesa_pd/collision_detection/Support.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace collision_detection {

/*!
 * \brief The Expanding-Polytope Algorithm.
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
   static constexpr real_t contactThreshold = real_t(1e-8);

   //**Query functions*****************************************************************************
   /*!\name Query functions */
   //@{
   bool doEPAcontactThreshold( Support &geom1,
                               Support &geom2,
                               const GJK& gjk,
                               Vec3& normal,
                               Vec3& contactPoint,
                               real_t& penetrationDepth);


   bool doEPAcontactThreshold( Support &geom1,
                               Support &geom2,
                               const GJK& gjk,
                               Vec3& normal,
                               Vec3& contactPoint,
                               real_t& penetrationDepth,
                               real_t eps_rel);
   
   bool doEPAmargin( Support &geom1,
                     Support &geom2,
                     const GJK& gjk,
                     Vec3& normal,
                     Vec3& contactPoint,
                     real_t& penetrationDepth,
                     real_t margin);

   bool doEPA( Support &geom1,
               Support &geom2,
               const GJK& gjk,
               Vec3& normal,
               Vec3& contactPoint,
               real_t& penetrationDepth,
               real_t margin,
               real_t eps_rel );



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
   inline void pushSupportMargin(const Support &geom1,
                                 const Support &geom2,
                                 const Vec3& dir,
                                 const real_t margin,
                                 std::vector<Vec3>& epaVolume,
                                 std::vector<Vec3>& supportA,
                                 std::vector<Vec3>& supportB);

   inline void replaceSupportMargin(const Support &geom1,
                                    const Support &geom2,
                                    const Vec3& dir,
                                    const real_t margin,
                                    std::vector<Vec3>& epaVolume,
                                    std::vector<Vec3>& supportA,
                                    std::vector<Vec3>& supportB,
                                    size_t indexToReplace);

   inline void removeSupportMargin(std::vector<Vec3>& epaVolume, std::vector<Vec3>& supportA, std::vector<Vec3>& supportB);

   inline bool originInTetrahedron            ( const Vec3& A, const Vec3& B, const Vec3& C,
                                                const Vec3& D );
   inline bool originInTetrahedronVolumeMethod( const Vec3& A, const Vec3& B, const Vec3& C,
                                                const Vec3& D );
   inline bool pointInTetrahedron             ( const Vec3& A, const Vec3& B, const Vec3& C,
                                                const Vec3& D, const Vec3& point );
   bool searchTetrahedron              (Support &geom1,
                                        Support &geom2,
                                        std::vector<Vec3>& epaVolume,
                                        std::vector<Vec3>& supportA,
                                        std::vector<Vec3>& supportB,
                                        EPA_EntryBuffer& entryBuffer,
                                        real_t margin );

   void createInitialTetrahedron       ( size_t top,
                                         size_t frontLeft,
                                         size_t frontRight,
                                         size_t back,
                                         std::vector<Vec3>& epaVolume,
                                         EPA_EntryBuffer& entryBuffer );

   void createInitialSimplex           ( size_t numPoints,
                                         Support &geom1,
                                         Support &geom2,
                                         std::vector<Vec3>& supportA,
                                         std::vector<Vec3>& supportB,
                                         std::vector<Vec3>& epaVolume,
                                         EPA_EntryBuffer& entryBuffer,
                                         real_t margin );
   inline real_t calculateCircle              ( const Vec3& A,
                                                const Vec3& B,
                                                const Vec3& C,
                                                const Vec3& D,
                                                Vec3& center );
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
   ///Returns the index of the internal vertex i(=0,1,2) within the EPA scope.
   inline size_t      operator[]( size_t i )                           const {return indices_[i];}
   ///Returns the point closest to the origin of the affine hull of the triangle, which is also the normal.
   inline const Vec3& getClosest()                                     const {return closest_;}
   ///Returns the normal of the triangle. Normal is not normalized!
   inline const Vec3& getNormal()                                      const {return normal_;}
   inline Vec3        getClosestPoint(const std::vector<Vec3>& points) const;
   ///Returns the squared distance to the closest to the origin of the affine hull of the triangle.
   inline real_t      getSqrDist()                                     const {return sqrDist_;}
   ///Returns true if the triangle is no longer part of the EPA polygon.
   inline bool        isObsolete()                                     const {return obsolete_;}
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
   ///Return the triangle this edge belongs to.
   EPA_Triangle* getTriangle() const {return triangle_;}
   ///Get the Index of this edge in its triangle.
   size_t        getIndex()    const {return startIdx_;}
   ///Return the start point index  of an edge.
   size_t        getStart()    const {return (*triangle_)[startIdx_];}
   ///Return the end point index of an edge.
   size_t        getEnd()      const {return (*triangle_)[(startIdx_+1) % 3];}
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
   ///Compare two triangles by their distance.
   inline bool operator()( const EPA_Triangle *tria1,
                           const EPA_Triangle *tria2 )
   { return tria1->getSqrDist() > tria2->getSqrDist(); }
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
inline EPA::EPA_Edge::EPA_Edge( EPA_Triangle* triangle,
                                size_t index )
   : triangle_(triangle)
   , startIdx_(index)
{
}
//*************************************************************************************************


//=================================================================================================
//
//  EPA::EPA_TRIANGLE GET FUNCTIONS
//
//=================================================================================================


//*************************************************************************************************
/*! \brief Calculates the corresponding closest point from the given points, using barycentric coordinates.
 */
inline Vec3 EPA::EPA_Triangle::getClosestPoint(const std::vector<Vec3>& points) const
{
   return   bar_[0] * points[indices_[0]]
         + bar_[1] * points[indices_[1]]
         + bar_[2] * points[indices_[2]];

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
inline void EPA::pushSupportMargin(const Support &geom1,
                                   const Support &geom2,
                                   const Vec3& dir,
                                   const real_t margin,
                                   std::vector<Vec3>& epaVolume,
                                   std::vector<Vec3>& supportA,
                                   std::vector<Vec3>& supportB)
{
   Vec3 ndir;
   if(floatIsEqual(dir.sqrLength(), real_t(1.0))){
      ndir = dir.getNormalizedIfNotZero();
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
inline void EPA::replaceSupportMargin(const Support &geom1,
                                      const Support &geom2,
                                      const Vec3& dir,
                                      const real_t margin,
                                      std::vector<Vec3>& epaVolume,
                                      std::vector<Vec3>& supportA,
                                      std::vector<Vec3>& supportB,
                                      size_t indexToReplace)
{
   Vec3 ndir;
   if(floatIsEqual(dir.sqrLength(), real_t(1.0))){
      ndir = dir.getNormalizedIfNotZero();
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
inline void EPA::removeSupportMargin(std::vector<Vec3>& epaVolume,
                                     std::vector<Vec3>& supportA,
                                     std::vector<Vec3>& supportB)
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
