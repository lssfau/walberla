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
//! \file TimeStep.h
//! \ingroup pe_coupling
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "core/timing/TimingTree.h"

#include "domain_decomposition/BlockStorage.h"

#include "pe/cr/ICR.h"
#include "pe/rigidbody/BodyIterators.h"
#include "pe/synchronization/SyncForces.h"

#include <boost/function.hpp>

#include <map>


namespace walberla {
namespace pe_coupling {

class TimeStep
{
public:

    explicit TimeStep( const shared_ptr<StructuredBlockStorage> & blockStorage,
                       const BlockDataID & bodyStorageID,
                       pe::cr::ICR & collisionResponse,
                       boost::function<void (void)> synchronizeFunc,
                       const real_t timestep = real_t(1), const uint_t intermediateSteps = uint_c(1) )
        : timestep_( timestep )
        , intermediateSteps_( ( intermediateSteps == 0 ) ? uint_c(1) : intermediateSteps )
        , blockStorage_( blockStorage )
        , bodyStorageID_(bodyStorageID)
        , collisionResponse_(collisionResponse)
        , synchronizeFunc_(synchronizeFunc)
    {}

    void operator()()
    {
        if( intermediateSteps_ == 1 )
        {
            collisionResponse_.timestep( timestep_ );
            synchronizeFunc_( );
        }
        else
        {
            // sum up all forces from shadow copies on local body (owner)
            pe::reduceForces( blockStorage_->getBlockStorage(), bodyStorageID_ );
            // send total forces to shadow owners
            pe::distributeForces( blockStorage_->getBlockStorage(), bodyStorageID_ );

            // generate map from all known bodies (process local) to total forces
            std::map< walberla::id_t, std::vector< real_t > > forceMap;
            for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
            {
                for( auto bodyIt = pe::BodyIterator::begin(*blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
                {
                    auto & f = forceMap[ bodyIt->getSystemID() ];

                    const auto & force = bodyIt->getForce();
                    f.push_back( force[0] );
                    f.push_back( force[1] );
                    f.push_back( force[2] );

                    const auto & torque = bodyIt->getTorque();
                    f.push_back( torque[0] );
                    f.push_back( torque[1] );
                    f.push_back( torque[2] );

                    if ( bodyIt->isRemote() )
                    {
                        bodyIt->resetForceAndTorque();
                    }
                }
            }

            // perform pe time steps
            const real_t subTimestep = timestep_ / real_c( intermediateSteps_ );
            for( uint_t i = 0; i != intermediateSteps_; ++i )
            {
                // in the first set forces on local bodies are already set by force synchronization
                if( i != 0 ) {
                    for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
                    {
                        for( auto body = pe::LocalBodyIterator::begin(*blockIt, bodyStorageID_); body != pe::LocalBodyIterator::end(); ++body )
                        {
                            const auto & f = forceMap[ body->getSystemID() ];
                            body->addForce ( f[0], f[1], f[2] );
                            body->addTorque( f[3], f[4], f[5] );
                        }
                    }
                }

                collisionResponse_.timestep( subTimestep );
                synchronizeFunc_( );
            }
        }
    }

protected:

    const real_t timestep_;

    const uint_t intermediateSteps_;

    shared_ptr<StructuredBlockStorage> blockStorage_;
    const BlockDataID &  bodyStorageID_;

    pe::cr::ICR & collisionResponse_;
    boost::function<void (void)> synchronizeFunc_;

}; // class TimeStep



} // namespace pe_coupling
} // namespace walberla
