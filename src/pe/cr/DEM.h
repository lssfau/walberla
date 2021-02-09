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
//! \file DEM.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file for the DEM solver
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "ICR.h"
#include "Integrators.h"
#include "ContactResolvers.h"
#include "pe/Types.h"

#include "domain_decomposition/BlockStorage.h"

namespace walberla {
namespace pe {
namespace cr {

/**
 * \ingroup pe
 */
template< typename Integrator, typename ContactResolver >
class DEMSolver : public ICR
{
public:
   DEMSolver( const Integrator & integrate, const ContactResolver & resolveContact
            , const shared_ptr<BodyStorage>&    globalBodyStorage
            , const shared_ptr<BlockStorage>&   blockStorage
            , domain_decomposition::BlockDataID storageID
            , domain_decomposition::BlockDataID ccdID
            , domain_decomposition::BlockDataID fcdID
            , WcTimingTree*                     tt = nullptr);

   /// forwards to timestep
   /// Convenience operator to make class a functor.
   void operator()(const real_t dt) { timestep(dt); }
   /// Advances the simulation dt seconds.
   void timestep( const real_t dt ) override;

   inline Integrator                getIntegrator()                const { return integrate_; }
   inline ContactResolver           getContactResolver()           const { return resolveContact_; }
   inline real_t            getMaximumPenetration()        const override { return maxPenetration_; }
   inline size_t            getNumberOfContacts()          const override { return numberOfContacts_; }
   inline size_t            getNumberOfContactsTreated()   const override { return numberOfContactsTreated_; }
private:
   Integrator                        integrate_;
   ContactResolver                   resolveContact_;
   shared_ptr<BodyStorage>           globalBodyStorage_;
   shared_ptr<BlockStorage>          blockStorage_;
   domain_decomposition::BlockDataID storageID_;
   domain_decomposition::BlockDataID ccdID_;
   domain_decomposition::BlockDataID fcdID_;
   WcTimingTree*                     tt_;

   real_t                            maxPenetration_;
   size_t                            numberOfContacts_;
   size_t                            numberOfContactsTreated_;
};

class DEM : public DEMSolver<IntegrateImplicitEuler, ResolveContactSpringDashpotHaffWerner>
{
public:
   DEM(  const shared_ptr<BodyStorage>&    globalBodyStorage
       , const shared_ptr<BlockStorage>&   blockStorage
       , domain_decomposition::BlockDataID storageID
       , domain_decomposition::BlockDataID ccdID
       , domain_decomposition::BlockDataID fcdID
       , WcTimingTree*                     tt = nullptr)
   : DEMSolver<IntegrateImplicitEuler, ResolveContactSpringDashpotHaffWerner>(
              IntegrateImplicitEuler(), ResolveContactSpringDashpotHaffWerner(),
              globalBodyStorage, blockStorage, storageID, ccdID, fcdID, tt )
   {
   }
};

}  // namespace cr
}  // namespace pe
}  // namespace walberla

#include "DEM.impl.h"
