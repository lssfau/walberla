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

#include "EPA.h"
#include <core/math/Constants.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace collision_detection {

//=================================================================================================
//
//  EPA::EPA_TRIANGLE CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*! \brief Construct a new EPA_Triangle.
 *  \param a First point index
 *  \param b Second point index
 *  \param c Third point index
 *  \param points Vector with all points
 */
EPA::EPA_Triangle::EPA_Triangle( size_t a,
                                 size_t b,
                                 size_t c,
                                 const std::vector<Vec3>& points )
{
   const Vec3& A = points[a];
   const Vec3& B = points[b];
   const Vec3& C = points[c];

   indices_[0] = a;
   indices_[1] = b;
   indices_[2] = c;

   //calculate the closest point to the origin
   //Real-Time Collsion Buch Seite 137
   Vec3 ab = B-A;
   Vec3 ac = C-A;
   //Vec3 bc = C-B;

   normal_ = ab % ac;
   Vec3 nT = normal_;

   //
   real_t vc = nT * (A % B);
   real_t va = nT * (B % C);
   real_t vb = nT * (C % A);
   real_t denom = real_t(1.0) / (va + vb + vc);

   bar_[0] = va * denom;
   bar_[1] = vb * denom;
   bar_[2] = real_t(1.0) - bar_[0] - bar_[1];

   closest_ = bar_[0] * A + bar_[1] * B + bar_[2] * C;

   //sqrDist=key is square distance of v to origin
   sqrDist_ = closest_.sqrLength();

   //adjoined triangles not set yet
   adjTriangle_[0] = adjTriangle_[1] = adjTriangle_[2] = nullptr;
   adjEdges_[0]    = adjEdges_[1]    = adjEdges_[2] = 4;

   obsolete_ = false;
}

//=================================================================================================
//
//  EPA::EPA_TRIANGLE UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \brief Sets the link of this triangles edge0 neighbor to tria and vice versa.
 */
inline bool EPA::EPA_Triangle::link( size_t edge0, EPA_Triangle* tria, size_t edge1 )
{
   WALBERLA_ASSERT_LESS(edge0, 3, "link: invalid edge index");
   WALBERLA_ASSERT_LESS(edge1, 3, "link: invalid edge index");

   adjTriangle_[edge0] = tria;
   adjEdges_[edge0] = edge1;
   tria->adjTriangle_[edge1] = this;
   tria->adjEdges_[edge1] = edge0;

   bool b = indices_[edge0]       == tria->indices_[(edge1+1)%3] &&
            indices_[(edge0+1)%3] == tria->indices_[edge1];
   return b;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Fills edgeBuffer with the CCW contour of triangles not seen from point w which is in normal direction of the triangle.
 */
inline void EPA::EPA_Triangle::silhouette( const Vec3& w, EPA_EdgeBuffer& edgeBuffer )
{
   //std::cerr << "Starting Silhoutette search on Triangle {" << indices_[0] << "," << indices_[1] << "," << indices_[2] << "}" << std::endl;
   edgeBuffer.clear();
   obsolete_ = true;

   adjTriangle_[0]->silhouette(adjEdges_[0], w, edgeBuffer);
   adjTriangle_[1]->silhouette(adjEdges_[1], w, edgeBuffer);
   adjTriangle_[2]->silhouette(adjEdges_[2], w, edgeBuffer);
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Recursive silhuette finding method.
 */
void EPA::EPA_Triangle::silhouette( size_t index, const Vec3& w,
                                    EPA_EdgeBuffer& edgeBuffer )
{
   if (!obsolete_) {
      real_t test = (closest_ * w);
      if (test < sqrDist_) {
         edgeBuffer.push_back(EPA_Edge(this, index));
      }
      else {
         obsolete_ = true; // Facet is visible
         size_t next = (index+1) % 3;
         adjTriangle_[next]->silhouette(adjEdges_[next], w, edgeBuffer);
         next = (next+1) % 3;
         adjTriangle_[next]->silhouette(adjEdges_[next], w, edgeBuffer);
      }
   }
}
//*************************************************************************************************

//=================================================================================================
//
//  EPA QUERY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
//EPA Precision default values for different data types
template< class T > struct EpsilonRelEPA;
template<> struct EpsilonRelEPA<       float > { static const       float value; };
template<> struct EpsilonRelEPA<      double > { static const      double value; };
template<> struct EpsilonRelEPA< long double > { static const long double value; };

const       float EpsilonRelEPA<       float >::value = static_cast<       float >(1e-4);
const      double EpsilonRelEPA<      double >::value = static_cast<      double >(1e-6);
const long double EpsilonRelEPA< long double >::value = static_cast< long double >(1e-6);

//*************************************************************************************************
/*! \brief Does an EPA computation with contactthreshold added. Use Default relative Error.
 */
bool EPA::doEPAcontactThreshold( Support &geom1,
                                 Support &geom2,
                                 const GJK& gjk,
                                 Vec3& retNormal,
                                 Vec3& contactPoint,
                                 real_t& penetrationDepth){

   //Default relative epsilon
   return doEPA(geom1, geom2, gjk, retNormal, contactPoint, penetrationDepth, contactThreshold, EpsilonRelEPA<real_t>::value);
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Does an EPA computation with contactThreshold added. Relative Error can be specified.
 */
bool EPA::doEPAcontactThreshold( Support &geom1,
                                 Support &geom2,
                                 const GJK& gjk,
                                 Vec3& retNormal,
                                 Vec3& contactPoint,
                                 real_t& penetrationDepth,
                                 real_t eps_rel){

   return doEPA(geom1, geom2, gjk, retNormal, contactPoint, penetrationDepth, contactThreshold, eps_rel);
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Does an EPA computation with margin added. Use Default relative Error.
 */
bool EPA::doEPAmargin( Support &geom1,
                       Support &geom2,
                       const GJK& gjk,
                       Vec3& retNormal,
                       Vec3& contactPoint,
                       real_t& penetrationDepth,
                       real_t margin){
   //Default relative epsilon
   return doEPA(geom1, geom2, gjk, retNormal, contactPoint, penetrationDepth, margin, EpsilonRelEPA<real_t>::value);
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Does an epa computation with contact margin added and specified realtive error.
 */
bool EPA::doEPA( Support &geom1,
                 Support &geom2,
                 const GJK& gjk,
                 Vec3& retNormal,
                 Vec3& contactPoint,
                 real_t& penetrationDepth,
                 real_t margin,
                 real_t eps_rel )
{
   //have in mind that we use a support mapping which blows up the objects a wee bit so
   //zero penetraion aka toching contact means that the original bodies have a distance of 2*margin between them

   //Set references to the results of GJK
   size_t     numPoints( static_cast<size_t>( gjk.getSimplexSize() ) );
   std::vector<Vec3> epaVolume( gjk.getSimplex() );
   std::vector<Vec3> supportA ( gjk.getSupportA() );
   std::vector<Vec3> supportB ( gjk.getSupportB() );

   Vec3 support;

   epaVolume.reserve( maxSupportPoints_ );
   supportA.reserve ( maxSupportPoints_ );
   supportB.reserve ( maxSupportPoints_ );

   EPA_EntryBuffer entryBuffer;
   entryBuffer.reserve(maxTriangles_);

   EPA_EntryHeap entryHeap;
   entryHeap.reserve(maxTriangles_);

   EPA_EdgeBuffer edgeBuffer;
   edgeBuffer.reserve(20);

   real_t lowerBoundSqr = math::Limits<real_t>::inf();
   real_t upperBoundSqr = math::Limits<real_t>::inf();

   //create an Initial simplex
   if(numPoints == 1) {
      //If the GJK-Simplex contains only one point, it must be the origin and it must be on the boundary of the CSO.
      //This means the enlarged bodies are in touching contact and the original bodies do not intersect.
      return false;
   }
   else {
      createInitialSimplex(numPoints, geom1, geom2, supportA, supportB, epaVolume, entryBuffer, margin);
   }

   for(EPA_EntryBuffer::iterator it=entryBuffer.begin(); it != entryBuffer.end(); ++it) {
      if(it->isClosestInternal()) {
         entryHeap.push_back(&(*it));
      }
   }

   if(entryHeap.empty()) {
      //unrecoverable error.
      return false;
   }

   std::make_heap(entryHeap.begin(), entryHeap.end(), EPA::EPA_TriangleComp());
   EPA_Triangle* current = nullptr;

   numIterations_ = 0;
   //EPA Main-Loop
   do {
      ++numIterations_;
      std::pop_heap(entryHeap.begin(), entryHeap.end(), EPA::EPA_TriangleComp());
      current = entryHeap.back();
      entryHeap.pop_back();
      if(!current->isObsolete()) {
         lowerBoundSqr = current->getSqrDist();

         if(epaVolume.size() == maxSupportPoints_) {
            WALBERLA_ASSERT(false, "Support point limit reached.");
            break;
         }

         // Compute new support direction
         // if origin is contained in plane, use out-facing normal.
         Vec3 normal;
         if(current->getSqrDist() < real_comparison::Epsilon<real_t>::value*real_comparison::Epsilon<real_t>::value){
            if(current->getNormal().sqrLength() < real_comparison::Epsilon<real_t>::value*real_comparison::Epsilon<real_t>::value){
               break;
            }

            normal = current->getNormal().getNormalizedOrZero();
         }else{
            normal = current->getClosest().getNormalizedOrZero();
         }
         //std::cerr << "Current Closest: " << current->getClosest();
         //std::cerr << "New support direction: " <<  normal << std::endl;

         pushSupportMargin(geom1, geom2, normal, margin, epaVolume, supportA, supportB);
         support = epaVolume.back();

         numPoints++;

         real_t farDist = support * normal; //not yet squared

         WALBERLA_ASSERT_GREATER(farDist, real_t(0.0), "EPA support mapping gave invalid point in expansion direction");
         //std::cerr << "New upper bound: " <<  farDist*farDist << std::endl;
         upperBoundSqr = std::min(upperBoundSqr, farDist*farDist);



         //Try to approximate the new surface with a sphere
         if (bUseSphereOptimization_)
         {
            Vec3 ctr;
            real_t radius2 = calculateCircle(support,
                                             epaVolume[(*current)[0]],
                  epaVolume[(*current)[1]],
                  epaVolume[(*current)[2]],
                  ctr);
            if(radius2 > real_t(0.0)){ //if a Sphere exists
               // std::cerr << "Circle created with center at " << ctr << ". r2=" << radius2 << std::endl;
               real_t center_len = ctr.length();
               real_t circle_dist = (std::sqrt(radius2) - center_len); //Distance from center to the spheres surface

               // Check if the circle matches the bounds given by EPA and limit max error to ca. 5%
               if (circle_dist*circle_dist <= upperBoundSqr &&
                   circle_dist*circle_dist >= lowerBoundSqr &&
                   (circle_dist*circle_dist) < real_t(1.10) * lowerBoundSqr &&
                   !floatIsEqual(center_len, real_t(0.0)) &&
                   circle_dist > real_t(0.0)) // In case of numerical errors, this can be the case
               {
                  const auto ilen = real_t(1.0) / center_len;
                  ctr *= -ilen;
                  pushSupportMargin(geom1, geom2, ctr, margin, epaVolume, supportA, supportB);
                  support = epaVolume.back();
                  // Check if support is in expected direction
                  real_t supp_dist = support.length();
                  if(floatIsEqual((support % ctr).sqrLength()/support.sqrLength(), real_t(0.0)) &&
                     supp_dist*supp_dist <= upperBoundSqr &&
                     supp_dist*supp_dist >= lowerBoundSqr)
                  {
                     //Accept sphere
                     contactPoint = real_t(0.5) * (supportA.back() + supportB.back());
                     penetrationDepth = -supp_dist + real_t(2.0) * margin;
                     retNormal = -ctr;
                     
                     return penetrationDepth < contactThreshold;
                  } else {
                     //Reject sphere
                     removeSupportMargin(epaVolume, supportA, supportB);
                     support = epaVolume.back();
                  }
               }
            }
         }

         //terminating criteria's
         //- we found that the two bounds are close enough
         //- the added support point was already in the epaVolume
         if(upperBoundSqr <= (real_t(1.0)+eps_rel)*(real_t(1.0)+eps_rel)*lowerBoundSqr
               || support == epaVolume[(*current)[0]]
               || support == epaVolume[(*current)[1]]
               || support == epaVolume[(*current)[2]])
         {
            //std::cerr << "Tolerance reached." << std::endl;
            break;
         }

         // Compute the silhouette cast by the new vertex
         // Note that the new vertex is on the positive side
         // of the current triangle, so the current triangle
         // will not be in the convex hull. Start local search
         // from this facet.

         current->silhouette(support, edgeBuffer);
         if(edgeBuffer.size() < 3 ) {
            return false;
         }

         if(entryBuffer.size() == maxSupportPoints_) {
            //"out of memory" so stop here
            //std::cerr << "Memory Limit reached." << std::endl;
            break;
         }

         EPA_EdgeBuffer::const_iterator it = edgeBuffer.begin();
         entryBuffer.push_back(EPA_Triangle(it->getEnd(), it->getStart(), epaVolume.size()-1, epaVolume));

         EPA_Triangle* firstTriangle = &(entryBuffer.back());
         //if it is expanding candidate add to heap
         //std::cerr << "Considering Triangle (" << firstTriangle->getSqrDist() << ") {"  << (*firstTriangle)[0] <<  "," << (*firstTriangle)[1] <<  ","<< (*firstTriangle)[2] << "} ("<< epaVolume[(*firstTriangle)[0]] * firstTriangle->getNormal()<< ")" << std::endl;
         if(epaVolume[(*firstTriangle)[0]] * firstTriangle->getNormal() < real_t(0.0)){
            //the whole triangle is on the wrong side of the origin.
            //This is a numerical error and will produce wrong results, if the search is continued. Stop here.
            break;
         }
         if(firstTriangle->isClosestInternal()
               && firstTriangle->getSqrDist() > lowerBoundSqr
               && firstTriangle->getSqrDist() < upperBoundSqr)
         {
            entryHeap.push_back(firstTriangle);
            std::push_heap(entryHeap.begin(), entryHeap.end(), EPA::EPA_TriangleComp());
         }

         firstTriangle->link(0, it->getTriangle(), it->getIndex());

         EPA_Triangle* lastTriangle = firstTriangle;

         ++it;
         for(; it != edgeBuffer.end(); ++it){
            if(entryBuffer.size() == maxSupportPoints_) {
               //"out of memory" so stop here
               break;
            }

            entryBuffer.push_back(EPA_Triangle(it->getEnd(), it->getStart(), epaVolume.size()-1, epaVolume));
            EPA_Triangle* newTriangle = &(entryBuffer.back());

            //std::cerr << "Considering Triangle (" << newTriangle->getSqrDist() << ") {"  << (*newTriangle)[0] <<  "," << (*newTriangle)[1] <<  ","<< (*newTriangle)[2] << "} ("<< epaVolume[(*newTriangle)[0]] * newTriangle->getNormal() << ")" << std::endl;

            if(epaVolume[(*newTriangle)[0]] * newTriangle->getNormal() < real_t(0.0)){
               //the whole triangle is on the wrong side of the origin.
               //This is an error.
               break;
            }
            //if it is expanding candidate add to heap
            if(newTriangle->isClosestInternal()
                  &&  newTriangle->getSqrDist() > lowerBoundSqr
                  &&  newTriangle->getSqrDist() < upperBoundSqr)
            {
               entryHeap.push_back(newTriangle);
               std::push_heap(entryHeap.begin(), entryHeap.end(), EPA::EPA_TriangleComp());
            }

            if(!newTriangle->link(0, it->getTriangle(), it->getIndex())) {
               break;
            }

            if(!newTriangle->link(2, lastTriangle, 1)) {
               break;
            }

            lastTriangle = newTriangle;
         }

         if(it != edgeBuffer.end()) {
            //For some reason the silhouette couldn't be processed completely
            //so we stop here and take the last result
            break;
         }

         firstTriangle->link(2, lastTriangle, 1);
      }
   } while (!entryHeap.empty() && entryHeap[0]->getSqrDist() <= upperBoundSqr);

   //Normal must be inverted
   retNormal   = -current->getClosest().getNormalizedOrZero();

   //Calculate Witness points
   const Vec3 wittnessA = current->getClosestPoint(supportA);
   const Vec3 wittnessB = current->getClosestPoint(supportB);
   contactPoint = real_t(0.5) * (wittnessA + wittnessB);

   //Penetration Depth
   penetrationDepth = -(current->getClosest().length() - real_t(2.0) * margin);

   /*std::cerr << "normal=" << retNormal <<std::endl;
   std::cerr << "close =" << current->getClosest() << std::endl;
   std::cerr << "diff  =" << wittnesA - wittnesB  <<std::endl;
   std::cerr << "wittnesA    =" << wittnesA <<std::endl;
   std::cerr << "wittnesB    =" << wittnesB <<std::endl;
   std::cerr << "contactPoint=" << contactPoint << std::endl;
   std::cerr << "penDepth=" << penetrationDepth  <<std::endl;
   std::cerr << "lowerBound=" << sqrt(lowerBoundSqr) <<std::endl;
   std::cerr << "curreBound=" << current->getClosest().length() << std::endl;
   std::cerr << "upperBound=" << sqrt(upperBoundSqr) <<std::endl;
   std::cerr << "Heap Size=" << entryHeap.size() << std::endl;
   std::cerr << "entryHeap[0]->getSqrDist()=" << entryHeap[0]->getSqrDist() << std::endl;*/
   //std::cout << "EPA penetration depth: " << penetrationDepth <<  std::endl;

   return penetrationDepth < contactThreshold;
}
//*************************************************************************************************

//=================================================================================================
//
//  EPA UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \brief Create a starting tetrahedron for EPA, from the GJK Simplex.
 */
inline void EPA::createInitialSimplex( size_t numPoints,
                                       Support &geom1,
                                       Support &geom2,
                                       std::vector<Vec3>& supportA,
                                       std::vector<Vec3>& supportB,
                                       std::vector<Vec3>& epaVolume,
                                       EPA_EntryBuffer& entryBuffer,
                                       real_t margin )
{
   switch(numPoints) {
   case 2:
   {
      //simplex is a line segement
      //add 3 points around the this segment
      //the COS is konvex so the resulting hexaheadron should be konvex too

      Vec3 d = epaVolume[1] - epaVolume[0];
      //find coordinate axis e_i which is furthest from paralell to d
      //and therefore d has the smallest abs(d[i])
      real_t abs0 = std::abs(d[0]);
      real_t abs1 = std::abs(d[1]);
      real_t abs2 = std::abs(d[2]);

      Vec3 axis;
      if( abs0 < abs1 && abs0 < abs2) {
         axis = Vec3(real_t(1.0), real_t(0.0), real_t(0.0));
      }
      else if( abs1 < abs0 && abs1 < abs2) {
         axis = Vec3(real_t(0.0), real_t(1.0), real_t(0.0));
      }
      else {
         axis = Vec3(real_t(0.0), real_t(0.0), real_t(1.0));
      }

      Vec3 direction1 = (d % axis).getNormalizedOrZero();
      Quat q(d, (real_t(2.0)/real_t(3.0)) * real_t(walberla::math::pi));
      Mat3 rot = q.toRotationMatrix();
      Vec3 direction2 = (rot*direction1).getNormalizedOrZero();
      Vec3 direction3 = (rot*direction2).getNormalizedOrZero();

      //add point in positive normal direction1
      pushSupportMargin(geom1, geom2, direction1, margin, epaVolume, supportA, supportB);
      //std::cerr << "S1: " << support1 << std::endl;

      //add point in negative normal direction2
      pushSupportMargin(geom1, geom2, direction2, margin, epaVolume, supportA, supportB);
      //std::cerr << "S2: " << support2 << std::endl;

      //add point in negative normal direction
      pushSupportMargin(geom1, geom2, direction3, margin, epaVolume, supportA, supportB);
      //std::cerr << "S3: " << support3 << std::endl;

      //Build the hexahedron as it is convex
      //epaVolume[1] = up
      //epaVolume[0] = down
      //epaVolume[2] = ccw1
      //epaVolume[3] = ccw2
      //epaVolume[4] = ccw3


      //check for containment inside
      if(originInTetrahedron(epaVolume[0], epaVolume[2], epaVolume[3], epaVolume[4]) || originInTetrahedron(epaVolume[1], epaVolume[2], epaVolume[3], epaVolume[4]) ){
         //insert triangle 1
         entryBuffer.push_back(EPA_Triangle(1, 2, 3, epaVolume)); //[0] up->ccw1->ccw2
         //insert triangle 2
         entryBuffer.push_back(EPA_Triangle(1, 3, 4, epaVolume)); //[1] up->ccw2->ccw3
         //insert triangle 3
         entryBuffer.push_back(EPA_Triangle(1, 4, 2, epaVolume)); //[2] up->ccw3->ccw1

         //link these 3 triangles
         entryBuffer[0].link(2, &(entryBuffer[1]), 0); //edge up->ccw1
         entryBuffer[1].link(2, &(entryBuffer[2]), 0); //edge up->ccw2
         entryBuffer[2].link(2, &(entryBuffer[0]), 0); //edge up->ccw3


         //insert triangle 4
         entryBuffer.push_back(EPA_Triangle(0, 2, 4, epaVolume)); //[3] down->ccw1->ccw3
         //insert triangle 5
         entryBuffer.push_back(EPA_Triangle(0, 4, 3, epaVolume)); //[4] down->ccw3->ccw2
         //insert triangle 6
         entryBuffer.push_back(EPA_Triangle(0, 3, 2, epaVolume)); //[5] down->ccw2->ccw1

         //link these 3 triangles
         entryBuffer[3].link(2, &(entryBuffer[4]), 0); //edge down->ccw3
         entryBuffer[4].link(2, &(entryBuffer[5]), 0); //edge down->ccw1
         entryBuffer[5].link(2, &(entryBuffer[3]), 0); //edge down->ccw1

         //link the two pyramids
         entryBuffer[0].link(1, &(entryBuffer[5]), 1); //edge ccw1->ccw2
         entryBuffer[1].link(1, &(entryBuffer[4]), 1); //edge ccw2->ccw3
         entryBuffer[2].link(1, &(entryBuffer[3]), 1); //edge ccw3->ccw1
      }else{
         //Apply iterative search
         removeSupportMargin(epaVolume, supportA, supportB); //remove 5th point.
         //Search starts from the remaining 4 points
         searchTetrahedron(geom1, geom2, epaVolume, supportA, supportB, entryBuffer, margin);
      }

      break;
   }
   case 3:
   {
      //simplex is a triangle, add tow points in positive and negative normal direction

      const Vec3& A = epaVolume[2];  //The Point last added to the simplex
      const Vec3& B = epaVolume[1];  //One Point that was already in the simplex
      const Vec3& C = epaVolume[0];  //One Point that was already in the simplex
      //ABC is a conterclockwise triangle

      const Vec3  AB  = B-A;       //The vector A->B
      const Vec3  AC  = C-A;       //The vector A->C
      const Vec3  ABC = (AB%AC).getNormalizedOrZero();     //The the normal pointing towards the viewer if he sees a CCW triangle ABC

      //add point in positive normal direction
      pushSupportMargin(geom1, geom2, ABC, margin, epaVolume, supportA, supportB);

      //add point in negative normal direction
      pushSupportMargin(geom1, geom2, -ABC, margin, epaVolume, supportA, supportB);
      //Vec3 support2 = epaVolume.back();

      //check if the hexahedron is convex aka check if a partial tetrahedron contains the last point
      if(pointInTetrahedron(epaVolume[3], epaVolume[4], epaVolume[0], epaVolume[2], epaVolume[1])) {
         //epaVolumne[1] is whithin the tetraheadron 3-4-0-2 so this is the epaVolume to take
         createInitialTetrahedron(3,4,0,2, epaVolume, entryBuffer);
      }
      else if(pointInTetrahedron(epaVolume[3], epaVolume[4], epaVolume[1], epaVolume[0], epaVolume[2])) {
         createInitialTetrahedron(3,4,1,0, epaVolume, entryBuffer);
      }
      else if(pointInTetrahedron(epaVolume[3], epaVolume[4], epaVolume[2], epaVolume[1], epaVolume[0])) {
         createInitialTetrahedron(3,4,2,1, epaVolume, entryBuffer);
      }
      else {
         //Build the hexahedron as it is convex
         //insert triangle 1
         entryBuffer.push_back(EPA_Triangle(3, 2, 1, epaVolume)); //[0] support1->A->B
         //insert triangle 2
         entryBuffer.push_back(EPA_Triangle(3, 1, 0, epaVolume)); //[1] support1->B->C
         //insert triangle 3
         entryBuffer.push_back(EPA_Triangle(3, 0, 2, epaVolume)); //[2] support1->C->A

         //link these 3 triangles
         entryBuffer[0].link(2, &(entryBuffer[1]), 0); //edge support1->A
         entryBuffer[1].link(2, &(entryBuffer[2]), 0); //edge support1->B
         entryBuffer[2].link(2, &(entryBuffer[0]), 0); //edge support1->C


         //insert triangle 4
         entryBuffer.push_back(EPA_Triangle(4, 2, 0, epaVolume)); //[3] support2->A->C
         //insert triangle 5
         entryBuffer.push_back(EPA_Triangle(4, 0, 1, epaVolume)); //[4] support2->C->B
         //insert triangle 6
         entryBuffer.push_back(EPA_Triangle(4, 1, 2, epaVolume)); //[5] support2->B->A

         //link these 3 triangles
         entryBuffer[3].link(2, &(entryBuffer[4]), 0); //edge support2->C
         entryBuffer[4].link(2, &(entryBuffer[5]), 0); //edge support2->B
         entryBuffer[5].link(2, &(entryBuffer[3]), 0); //edge support2->A

         //link the two pyramids
         entryBuffer[0].link(1, &(entryBuffer[5]), 1); //edge A->B
         entryBuffer[1].link(1, &(entryBuffer[4]), 1); //edge B->C
         entryBuffer[2].link(1, &(entryBuffer[3]), 1); //edge C->A

      }

      break;
   }
   case 4:
   {
      createInitialTetrahedron(3,2,1,0, epaVolume, entryBuffer);
      break;
   }
   default:
   {
      WALBERLA_ASSERT( false, "invalid number of simplex points in EPA" );
      break;
   }
   }
}

//*************************************************************************************************

//*************************************************************************************************
/*! \brief TODO
 *
 * see Book "collision detection in interactive 3D environments" page161
 * ATTENTION seems to have no consistent behavior on the surface and vertices
 */
inline bool EPA::originInTetrahedron( const Vec3& p0, const Vec3& p1, const Vec3& p2,
                                      const Vec3& p3 )
{
   Vec3 normal0T = (p1 -p0) % (p2-p0);
   if( (normal0T*p0 > real_t(0.0)) == (normal0T*p3 > real_t(0.0)) ) {
      return false;
   }
   Vec3 normal1T = (p2 -p1) % (p3-p1);
   if( (normal1T*p1 > real_t(0.0)) == (normal1T*p0 > real_t(0.0)) ) {
      return false;
   }
   Vec3 normal2T = (p3 -p2) % (p0-p2);
   if( (normal2T*p2 > real_t(0.0)) == (normal2T*p1 > real_t(0.0)) ) {
      return false;
   }
   Vec3 normal3T = (p0 -p3) % (p1-p3);
   return (normal3T*p3 > real_t(0.0)) != (normal3T*p2 > real_t(0.0));
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Retrurns true, if the origin lies in the tetrahedron ABCD.
 */
inline bool EPA::originInTetrahedronVolumeMethod( const Vec3& A, const Vec3& B, const Vec3& C,
                                                  const Vec3& D )
{
   const Vec3& aoT = A;
   if((aoT * (B % C)) <= real_t(0.0)) {
      //if volume of ABC and Origin <0.0 than the origin is on the wrong side of ABC
      //http://mathworld.wolfram.com/Tetrahedron.html volume formula
      return false;
   }
   if((aoT * (C % D)) <= real_t(0.0)) {
      return false;
   }
   if((aoT * (D % B)) <= real_t(0.0)) {
      return false;
   }
   if((B * (D % C)) <= real_t(0.0)) {
      return false;
   }
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Retrurns true, if a point lies in the tetrahedron ABCD.
 *  \param point The point to be checked for containment.
 */
inline bool EPA::pointInTetrahedron( const Vec3& A, const Vec3& B, const Vec3& C, const Vec3& D,
                                     const Vec3& point )
{
   return originInTetrahedronVolumeMethod( A-point, B-point, C-point, D-point );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief TODO
 * top, frontLeft ... are indices
 */
inline void EPA::createInitialTetrahedron( size_t top,
                                           size_t frontLeft,
                                           size_t frontRight,
                                           size_t back,
                                           std::vector<Vec3>& epaVolume,
                                           EPA_EntryBuffer& entryBuffer )
{
   //insert triangle 1
   entryBuffer.push_back(EPA_Triangle(top, frontLeft, frontRight, epaVolume)); //[0] vorne
   //insert triangle 2
   entryBuffer.push_back(EPA_Triangle(top, frontRight, back, epaVolume)); //[1] rechts hinten
   //insert triangle 3
   entryBuffer.push_back(EPA_Triangle(top, back, frontLeft, epaVolume)); //[2] links hinten
   //insert triangle 4
   entryBuffer.push_back(EPA_Triangle(back, frontRight, frontLeft, epaVolume)); //[3] unten

   //make links between the triangles
   entryBuffer[0].link(0, &(entryBuffer[2]), 2); //Kante vorne links
   entryBuffer[0].link(2, &(entryBuffer[1]), 0); //Kante vorne rechts
   entryBuffer[0].link(1, &(entryBuffer[3]), 1); //kante vorne unten

   entryBuffer[1].link(2, &(entryBuffer[2]), 0); //Kante hinten
   entryBuffer[1].link(1, &(entryBuffer[3]), 0); //kante rechts unten

   entryBuffer[2].link(1, &(entryBuffer[3]), 2); //kante links unten

}
//*************************************************************************************************

/*! \brief Search a tetrahedron that contains the origin.
 * Start with four arbitrary support points in epaVolume that form a
 * tetrahedron. (This needn't contain the origin.)
 * This algorithm will search and return an altered tetrahedron
 * containing the origin. Do only use this function if the object/body
 * certainly contains the origin!
 * \return True, if a tetrahedron was found. False if search has been aborted.
 */
inline bool EPA::searchTetrahedron(Support &geom1,
                                   Support &geom2,
                                   std::vector<Vec3>& epaVolume,
                                   std::vector<Vec3>& supportA,
                                   std::vector<Vec3>& supportB,
                                   EPA_EntryBuffer& entryBuffer,
                                   real_t margin )
{
   //Store the point no longer needed (0 if all points are needed, and origin is contained.)
   int loopCount = 0;
   int pointIndexToRemove = -1;
   Vec3 newSearchDirection;
   do{
      loopCount++;
      pointIndexToRemove = -1;
      //Check if opposite tetrahedron point and orign are on the same side
      //of the face. (for all faces)
      Vec3 normal0T = (epaVolume[1] -epaVolume[0]) % (epaVolume[2]-epaVolume[0]);
      real_t dot_val = normal0T*epaVolume[0];
      if( (normal0T*epaVolume[3] < dot_val) == (dot_val < real_t(0.0)) ) {
         pointIndexToRemove = 3;
         newSearchDirection = (normal0T*epaVolume[3] < dot_val) ? normal0T : -normal0T;
      }

      Vec3 normal1T = (epaVolume[2] -epaVolume[1]) % (epaVolume[3]-epaVolume[1]);
      dot_val = normal1T*epaVolume[1];
      if( (normal1T*epaVolume[0] < dot_val) == (dot_val < real_t(0.0)) ) {
         pointIndexToRemove = 0;
         newSearchDirection = (normal1T*epaVolume[0] < dot_val) ? normal1T : -normal1T;
      }

      Vec3 normal2T = (epaVolume[3] -epaVolume[2]) % (epaVolume[0]-epaVolume[2]);
      dot_val = normal2T*epaVolume[2];
      if( (normal2T*epaVolume[1] < dot_val) == (dot_val < real_t(0.0)) ) {
         pointIndexToRemove = 1;
         newSearchDirection = (normal2T*epaVolume[1] < dot_val) ? normal2T : -normal2T;
      }

      Vec3 normal3T = (epaVolume[0] -epaVolume[3]) % (epaVolume[1]-epaVolume[3]);
      dot_val = normal3T*epaVolume[3];
      if( (normal3T*epaVolume[2] < dot_val) == (dot_val < real_t(0.0)) ) {
         pointIndexToRemove = 2;
         newSearchDirection = (normal3T*epaVolume[2] < dot_val) ? normal3T : -normal3T;
      }
      //Origin not contained in tetrahedron.
      if(pointIndexToRemove != -1){
         if(loopCount > 50){
            return false;
         }
         //Get new support point and replace old.
         /*std::cerr << "Search Direction is: "<< newSearchDirection << std::endl;
                   std::cerr << "Projection of unnecc. point " << pointIndexToRemove << ": " << epaVolume[pointIndexToRemove] * newSearchDirection << std::endl;
                   std::cerr << "Projection of other points: " << epaVolume[(pointIndexToRemove+1)%4] * newSearchDirection << std::endl;*/
         newSearchDirection = newSearchDirection.getNormalizedOrZero();
         /*supportA[pointIndexToRemove] = geom1.supportContactThreshold(newSearchDirection);
                   supportB[pointIndexToRemove] = geom2.supportContactThreshold(-newSearchDirection);
                   epaVolume[pointIndexToRemove] = supportA[pointIndexToRemove] - supportB[pointIndexToRemove];*/
         replaceSupportMargin(geom1, geom2, newSearchDirection, margin, epaVolume, supportA, supportB, (size_t)pointIndexToRemove);
         //std::cerr << "Projection of new support point " << epaVolume[pointIndexToRemove] << ": " << epaVolume[pointIndexToRemove] * newSearchDirection << std::endl;

      }
   }
   while(pointIndexToRemove != -1);
   //std::cerr << "Found Tet after " << loopCount << " searches." << std::endl;

   //Build final tetrahedron
   Vec3 check_normal = (epaVolume[1] -epaVolume[0]) % (epaVolume[2]-epaVolume[0]);
   if(check_normal*epaVolume[3] > check_normal*epaVolume[0]){
      //p3 is behind.
      createInitialTetrahedron(1, 0, 2, 3, epaVolume, entryBuffer);
   }else{
      //p3 is in front
      createInitialTetrahedron(1, 3, 2, 0, epaVolume, entryBuffer);
   }
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \brief Calculate a Circle through the for Points A, B, C, D.
 * \param center Contains the center point of the circle after the call
 * \return The squared radius of the circle or a negative value if no such circle exists.
 */
inline real_t EPA::calculateCircle(const Vec3& A, const Vec3& B, const Vec3& C,
                                   const Vec3& D, Vec3& center ){
   const Vec3 n1(A-B);
   const Vec3 n2(A-C);
   const Vec3 n3(A-D);

   // Here we already see if such circle exists.
   const real_t det = n1 * (n2 % n3);
   if(floatIsEqual(det, real_t(0.0))){
      //no circle exists. Leave center untouched, and return -1.0
      return real_t(-1.0);
   }
   const real_t Alen = A.sqrLength();
   const real_t d1 = (Alen - B.sqrLength())*real_t(0.5);
   const real_t d2 = (Alen - C.sqrLength())*real_t(0.5);
   const real_t d3 = (Alen - D.sqrLength())*real_t(0.5);

   //Apply solution formula
   center = (real_t(1.0)/det)*(d1 * (n2 % n3) + d2 * (n3 % n1) + d3 * (n1 % n2));

   return (A - center).sqrLength();
}
//*************************************************************************************************

} //collision_detection
} //mesa_pd
} //walberla
