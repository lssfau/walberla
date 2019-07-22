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
//! \file LubricationCorrection.cpp
//! \ingroup pe_coupling
//! \author Kristina Pickl <kristina.pickl@fau.de>
//! \author Dominik Bartuschat
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "LubricationCorrection.h"

#include "core/logging/all.h"
#include "core/debug/Debug.h"

#include "pe/rigidbody/BodyIterators.h"
#include "pe/utility/Distance.h"

namespace walberla {
namespace pe_coupling {

void LubricationCorrection::operator ()()
{
   WALBERLA_LOG_PROGRESS( "Calculating Lubrication Force" );

   for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
   {
      // loop over all rigid bodies
      for( auto body1It = pe::BodyIterator::begin( *blockIt, bodyStorageID_ ); body1It != pe::BodyIterator::end(); ++body1It )
      {
         // lubrication forces for spheres
         if ( body1It->getTypeID() == pe::Sphere::getStaticTypeID() )
         {
            pe::SphereID sphereI = static_cast<pe::SphereID> ( body1It.getBodyID() );

            auto copyBody1It = body1It;
            // loop over all rigid bodies after current body1 to avoid double forces
            for( auto body2It = (++copyBody1It); body2It != pe::BodyIterator::end(); ++body2It )
            {
               // sphere-sphere lubrication
               if ( body2It->getTypeID() == pe::Sphere::getStaticTypeID() )
               {
                  pe::SphereID sphereJ = static_cast<pe::SphereID>( body2It.getBodyID() );
                  treatLubricationSphrSphr( sphereI, sphereJ, blockIt->getAABB() );
               }
            }
         }
      }

      // lubrication correction for local bodies with global bodies (for example sphere-plane)
      for( auto body1It = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID_ ); body1It != pe::LocalBodyIterator::end(); ++body1It )
      {
         if ( body1It->getTypeID() == pe::Sphere::getStaticTypeID() )
         {
            pe::SphereID sphereI = static_cast<pe::SphereID> ( body1It.getBodyID() );

            for (auto body2It = globalBodyStorage_->begin(); body2It != globalBodyStorage_->end(); ++body2It)
            {
               if ( body2It->getTypeID() == pe::Plane::getStaticTypeID() )
               {
                  // sphere-plane lubrication
                  pe::PlaneID planeJ = static_cast<pe::PlaneID>( body2It.getBodyID() );
                  treatLubricationSphrPlane( sphereI, planeJ );
               } else if ( body2It->getTypeID() == pe::Sphere::getStaticTypeID() )
               {
                  // sphere-sphere lubrication
                  pe::SphereID sphereJ = static_cast<pe::SphereID>( body2It.getBodyID() );
                  treatLubricationSphrSphr( sphereI, sphereJ, blockIt->getAABB() );
               }
            }
         }
      }
   }
}
//*****************************************************************************************************************************************



//////////////////////
// Helper Functions //
//////////////////////

void LubricationCorrection::treatLubricationSphrSphr( const pe::SphereID sphereI, const pe::SphereID sphereJ, const math::AABB & blockAABB )
{

   WALBERLA_ASSERT_UNEQUAL( sphereI->getSystemID(), sphereJ->getSystemID() );

   real_t gap = pe::getSurfaceDistance( sphereI, sphereJ );

   if ( gap > cutOffDistance_ || gap < real_t(0) )
   {
      WALBERLA_LOG_DETAIL("gap " << gap << " larger than cutOff " << cutOffDistance_ << " - ignoring pair");
      return;
   }

   if ( gap < minimalGapSize_ )
   {
      WALBERLA_LOG_DETAIL("gap " << gap << " smaller than minimal gap " << minimalGapSize_ << " - using minimal gap");
      gap = minimalGapSize_;
   }

   const pe::Vec3 &posSphereI = sphereI->getPosition();
   const pe::Vec3 &posSphereJ = sphereJ->getPosition();
   pe::Vec3 fLub(0);

   // compute (global) coordinate between spheres' centers of gravity
   pe::Vec3 midPoint( (posSphereI + posSphereJ ) * real_c(0.5) );

   // Let process on which midPoint lies do the lubrication correction
   // or the local process of sphereI if sphereJ is global
   if ( blockAABB.contains(midPoint) || sphereJ->isGlobal() )
   {
      fLub = compLubricationSphrSphr(gap, sphereI, sphereJ);
      sphereI->addForce( fLub);
      sphereJ->addForce(-fLub);

      WALBERLA_LOG_DETAIL( "Lubrication force on sphere " << sphereI->getID() << " from sphere " << sphereJ->getID() << " is:" << fLub);
      WALBERLA_LOG_DETAIL( "Lubrication force on sphere " << sphereJ->getID() << " from sphere " << sphereI->getID() << " is:" << -fLub << "\n");
   }

}

void LubricationCorrection::treatLubricationSphrPlane( const pe::SphereID sphereI, const pe::ConstPlaneID planeJ )
{

   real_t gap = pe::getSurfaceDistance( sphereI, planeJ );

   if ( gap > cutOffDistance_ || gap < real_t(0) )
   {
      WALBERLA_LOG_DETAIL("gap " << gap << " larger than cutOff " << cutOffDistance_ << " - ignoring pair");
      return;
   }

   if ( gap < minimalGapSize_ )
   {
      WALBERLA_LOG_DETAIL("gap " << gap << " smaller than minimal gap " << minimalGapSize_ << " - using minimal gap");
      gap = minimalGapSize_;
   }

   pe::Vec3 fLub = compLubricationSphrPlane( gap, sphereI, planeJ);

   WALBERLA_LOG_DETAIL( "Lubrication force on sphere " << sphereI->getID() << " to plane with id " << planeJ->getID() << " is:" << fLub << std::endl );
   sphereI->addForce( fLub );

}



//*****************************************************************************************************************************************
/*! \brief Computes lubrication correction force between spheres.
 * \ingroup pe_coupling
 *
 * Lubrication correction according to Ladd and Verberg, 2001
 * ("Lattice-Boltzmann Simulations of Particle-Fluid Suspensions")
 *
 * Note: Verified quantitatively by computation in spreadsheet
 * and qualitatively by considering direction of force for example setup.
 */
//*****************************************************************************************************************************************
pe::Vec3 LubricationCorrection::compLubricationSphrSphr( real_t gap, const pe::SphereID sphereI, const pe::SphereID sphereJ) const
{
   const pe::Vec3 &posSphereI = sphereI->getPosition();
   const pe::Vec3 &posSphereJ = sphereJ->getPosition();

   const pe::Vec3 &tmpVec = posSphereJ - posSphereI;
   const pe::Vec3 &rIJ    = tmpVec.getNormalized();

   real_t radiusSphereI = sphereI->getRadius();
   real_t radiusSphereJ = sphereJ->getRadius();

   pe::Vec3 velDiff(sphereI->getLinearVel() - sphereJ->getLinearVel());

   real_t length = velDiff * rIJ;

   real_t radiiSQR    = ( radiusSphereI * radiusSphereJ ) * ( radiusSphereI * radiusSphereJ );
   real_t radiiSumSQR = ( radiusSphereI + radiusSphereJ ) * ( radiusSphereI + radiusSphereJ );

   pe::Vec3 fLub = ( -real_t(6) * dynamicViscosity_ * walberla::math::pi * radiiSQR / radiiSumSQR * ( real_t(1) / gap - real_t(1) / cutOffDistance_) * length * rIJ);

   WALBERLA_LOG_DETAIL_SECTION()
   {
      std::stringstream ss;
      ss << "Sphere I: \n uid:" << sphereI->getID() << "\n";
      ss << "vel: "  << sphereI->getLinearVel() << "\n";
      ss << "rad: "  << radiusSphereI << "\n";
      ss << "pos: "  << posSphereI << "\n\n";

      ss << "Sphere J: \n uid:" << sphereJ->getID() << "\n";
      ss << "vel: "  << sphereJ->getLinearVel() << "\n";
      ss << "rad: "  << radiusSphereJ << "\n";
      ss << "pos: "  << posSphereJ << "\n\n";

      real_t distance = gap + radiusSphereI + radiusSphereJ;
      ss << "distance: "  << distance << "\n";
      ss << "viscosity: " << dynamicViscosity_ << "\n";

      ss << "gap: "     << gap << "\n";
      ss << "cutOff: "  << cutOffDistance_ << "\n";
      ss << "velDiff "  << velDiff << "\n";
      ss << "rIJ "      << rIJ << "\n\n";

      ss << "Resulting lubrication force: " << fLub << "\n";

      WALBERLA_LOG_DETAIL( ss.str() );
   }

   return fLub;
}
//*****************************************************************************************************************************************



//*****************************************************************************************************************************************
/*! \brief Computes lubrication correction force between sphere and wall.
 * \ingroup pe_coupling
 *
 * Lubrication correction according to Ladd and Verberg, 2001
 * ("Lattice-Boltzmann Simulations of Particle-Fluid Suspensions" )
 *
 * Note: Verified quantitatively by computation in spreadsheet
 * and qualitatively by considering direction of force for example setup.
 */
//*****************************************************************************************************************************************
pe::Vec3 LubricationCorrection::compLubricationSphrPlane( real_t gap, const pe::SphereID sphereI, const pe::ConstPlaneID planeJ) const
{
   const pe::Vec3 &posSphereI( sphereI->getPosition() );
   real_t radiusSphereI = sphereI->getRadius();

   const pe::Vec3 &planeNormal( planeJ->getNormal() ); // took negative of normal from sphere to plane (plane's normal) - sign cancels out anyway
   pe::Vec3 rIJ( planeNormal.getNormalized() );        // for safety reasons, normalize normal

   real_t length = sphereI->getLinearVel() * rIJ;

   real_t radiiSQR = radiusSphereI * radiusSphereI;

   pe::Vec3 fLub( -real_t(6) * dynamicViscosity_ * walberla::math::pi * radiiSQR * (real_t(1) / gap - real_t(1) / cutOffDistance_) * length * rIJ);

   WALBERLA_LOG_DETAIL_SECTION() {
      std::stringstream ss;
      ss << "Sphere I: \n uid:" << sphereI->getID() << "\n";
      ss << "vel: "  << sphereI->getLinearVel() << "\n";
      ss << "rad: "  << radiusSphereI << "\n";
      ss << "pos: "  << posSphereI << "\n\n";

      real_t distance = gap + radiusSphereI;
      ss << "distance: "  << distance << "\n";
      ss << "viscosity: " << dynamicViscosity_ << "\n";

      ss << "gap: "     << gap << "\n";
      ss << "cutOff: "  << cutOffDistance_ << "\n";
      ss << "velDiff "  << sphereI->getLinearVel() << "\n";
      ss << "rIJ "      << -rIJ << "\n\n";

      ss << "Resulting lubrication force: " << fLub << "\n";

      WALBERLA_LOG_DETAIL( ss.str() );
   }

   return fLub;
}


} // pe_coupling
} // walberla
