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
//! \file basic.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Header file which includes common pe headers!
//
//======================================================================================================================

#pragma once

//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "pe/Materials.h"
#include "pe/Types.h"

#include "pe/rigidbody/SetBodyTypeIDs.h"

#include "pe/rigidbody/StorageDataHandling.h"
#include "pe/ccd/HashGridsDataHandling.h"
#include "pe/fcd/SimpleFCDDataHandling.h"
#include "pe/bg/SimpleBGDataHandling.h"
#include "pe/cr/DEM.h"
#include "pe/cr/HCSITS.h"

#include "pe/rigidbody/BodyIterators.h"
#include "pe/rigidbody/BoxFactory.h"
#include "pe/rigidbody/CapsuleFactory.h"
#include "pe/rigidbody/CylindricalBoundaryFactory.h"
#include "pe/rigidbody/PlaneFactory.h"
#include "pe/rigidbody/SphereFactory.h"
#include "pe/rigidbody/UnionFactory.h"
#include "pe/rigidbody/EllipsoidFactory.h"

#include "pe/synchronization/SyncNextNeighbors.h"
#include "pe/synchronization/SyncShadowOwners.h"

#include "pe/utility/GetBody.h"
