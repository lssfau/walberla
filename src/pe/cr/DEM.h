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
//! \ingroup pe
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
#include "pe/Types.h"

#include "domain_decomposition/BlockStorage.h"

namespace walberla {
namespace pe {
namespace cr {

class DEM : public ICR
{
public:
   DEM(  const shared_ptr<BodyStorage>&    globalBodyStorage
       , const shared_ptr<BlockStorage>&   blockStorage
       , domain_decomposition::BlockDataID storageID
       , domain_decomposition::BlockDataID ccdID
       , domain_decomposition::BlockDataID fcdID
       , WcTimingTree*                     tt = NULL);

   /// forwards to timestep
   /// Convenience operator to make class a functor.
   void operator()(const real_t dt) { timestep(dt); }
   /// Advances the simulation dt seconds.
   void timestep( const real_t dt );

   virtual inline real_t            getMaximumPenetration()        const WALBERLA_OVERRIDE { return maxPenetration_; }
   virtual inline size_t            getNumberOfContacts()          const WALBERLA_OVERRIDE { return numberOfContacts_; }
   virtual inline size_t            getNumberOfContactsTreated()   const WALBERLA_OVERRIDE { return numberOfContactsTreated_; }
private:
   void resolveContact( ContactID c ) const;
   void move( BodyID id, real_t dt );

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

}  // namespace cr
} // namespace pe
}  // namespace walberla
