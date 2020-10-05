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
//! \file SimpleFCDDataHandling.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "SimpleFCD.h"

#include "blockforest/BlockDataHandling.h"

namespace walberla{
namespace pe{
namespace fcd {

template <typename BodyTypeTuple>
class SimpleFCDDataHandling : public blockforest::AlwaysInitializeBlockDataHandling<SimpleFCD<BodyTypeTuple> >{
public:
    SimpleFCD<BodyTypeTuple> * initialize( IBlock * const /*block*/ ) {return new SimpleFCD<BodyTypeTuple>();}
};

template <typename BodyTypeTuple>
[[deprecated("Use createGenericFCDDataHandling<BodyTypeTuple, AnalyticCollideFunctor>() instead")]]
shared_ptr<SimpleFCDDataHandling<BodyTypeTuple> > createSimpleFCDDataHandling();

template <typename BodyTypeTuple>
shared_ptr<SimpleFCDDataHandling<BodyTypeTuple> > createSimpleFCDDataHandling()
{
   return make_shared<SimpleFCDDataHandling<BodyTypeTuple> >( );
}

}
}
}
