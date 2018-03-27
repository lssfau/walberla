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
//! \file Integrators.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Integrators for the DEM collision resolver
//
//======================================================================================================================

#pragma once

#include "pe/Types.h"
#include "ICR.h"

namespace walberla {
namespace pe {
namespace cr {

//*************************************************************************************************
/*!\brief Integrate the trajectory of one body using implicit Euler.
*
* \param id Body ID.
* \param dt Time step size.
* \param solver Solver.
* \return void
*
* The implicit Euler algorithm, also known as backward Euler, is used. It is a first-order
* integrator that does conserves energy (i.e. it is symplectic.)
*/
class IntegrateImplicitEuler {
public:
   void operator()( BodyID id, real_t dt, ICR & solver ) const
   {
      // Calculating the linear acceleration by the equation
      //   force * m^(-1) + gravity
      const Vec3 vdot( id->getForce() * id->getInvMass() + solver.getGlobalLinearAcceleration() );

      // Calculating the angular acceleration by the equation
      //   R * Iinv * R^T * torque
      const Vec3 wdot( id->getInvInertia() * id->getTorque() );

      // Updating the linear velocity
      id->setLinearVel( id->getLinearVel() + vdot * dt );

      // Updating the angular velocity
      id->setAngularVel( id->getAngularVel() + wdot * dt );

      // Calculating the translational displacement
      id->setPosition( id->getPosition() + id->getLinearVel() * dt );

      // Calculating the rotation angle
      const Vec3 phi( id->getAngularVel() * dt );

      // Calculating the new orientation
      if (!floatIsEqual(phi.length(), 0))
         id->rotate( Quat( phi, phi.length() ) );
      WALBERLA_ASSERT_FLOAT_EQUAL( id->getRotation().getDeterminant(), real_c(1 ), "Corrupted rotation matrix determinant" );

      // Setting the axis-aligned bounding box
      id->calcBoundingBox();

      // Calculating the current motion of the body
      id->calcMotion();
   }
};

//*************************************************************************************************
/*!\brief Integrate the trajectory of one body using explicit Euler.
*
* \param id Body ID.
* \param dt Time step size.
* \param solver Solver.
* \return void
*
* The explicit Euler algorithm, also known as forward Euler, is used. It is a first-order
* integrator that does not conserve energy (i.e. it is not symplectic.)
*/
class IntegrateExplicitEuler {
public:
   void operator()( BodyID id, real_t dt, ICR & solver ) const
   {
      // Calculating the linear acceleration by the equation
      //   force * m^(-1) + gravity
      const Vec3 vdot( id->getForce() * id->getInvMass() + solver.getGlobalLinearAcceleration() );

      // Calculating the angular acceleration by the equation
      //   R * Iinv * R^T * torque
      const Vec3 wdot( id->getInvInertia() * id->getTorque() );

      // Calculating the translational displacement
      id->setPosition( id->getPosition() + id->getLinearVel() * dt + 0.5 * vdot * dt * dt );

      // Calculating the rotation angle
      const Vec3 phi( id->getAngularVel() * dt + 0.5 * wdot * dt * dt);

      // Calculating the new orientation
      if (!floatIsEqual(phi.length(), 0))
         id->rotate( Quat( phi, phi.length() ) );
      WALBERLA_ASSERT_FLOAT_EQUAL( id->getRotation().getDeterminant(), real_c(1 ), "Corrupted rotation matrix determinant" );

      // Updating the linear velocity
      id->setLinearVel( id->getLinearVel() + vdot * dt );

      // Updating the angular velocity
      id->setAngularVel( id->getAngularVel() + wdot * dt );

      // Setting the axis-aligned bounding box
      id->calcBoundingBox();

      // Calculating the current motion of the body
      id->calcMotion();
   }
};

}  // namespace cr
} // namespace pe
}  // namespace walberla
