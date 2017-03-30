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
//! \file SweepWrappers.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once


namespace walberla {
namespace lbm {



template< typename Kernel >
class StreamSweep
{
public:

   StreamSweep( const shared_ptr< Kernel > & kernel ) : kernel_( kernel ) {}

   void operator()( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) )
   {
      kernel_->stream( block, numberOfGhostLayersToInclude );
   }

private:

   shared_ptr< Kernel > kernel_;
};


template< typename Kernel >
inline StreamSweep< Kernel > makeStreamSweep( const shared_ptr< Kernel > & kernel ) { return StreamSweep<Kernel>( kernel ); }


template< typename Kernel >
class CollideSweep
{
public:

   CollideSweep( const shared_ptr< Kernel > & kernel ) : kernel_( kernel ) {}

   void operator()( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) )
   {
      kernel_->collide( block, numberOfGhostLayersToInclude );
   }

private:

   shared_ptr< Kernel > kernel_;
};

template< typename Kernel >
inline CollideSweep< Kernel > makeCollideSweep( const shared_ptr< Kernel > & kernel ) { return CollideSweep<Kernel>( kernel ); }



} // namespace lbm
} // namespace walberla
