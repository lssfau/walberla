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
//! \\file D3Q27BoundaryCollection.h
//! \\author lbmpy
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "domain_decomposition/IBlock.h"

#include "OutflowD3Q27.h"
#include "FixedDensityD3Q27.h"
#include "FreeSlipD3Q27.h"
#include "NoSlipD3Q27.h"
#include "UBBD3Q27.h"



namespace walberla{
namespace lbm {

template <typename FlagField_T>
class D3Q27BoundaryCollection
{
 public:
   enum Type { ALL = 0, INNER = 1, OUTER = 2 };


   D3Q27BoundaryCollection(const shared_ptr<StructuredBlockForest> & blocks, BlockDataID flagID_, BlockDataID pdfsID_, FlagUID domainUID_, double density, double u_x, double u_y, double u_z)
      : blocks_(blocks), flagID(flagID_), pdfsID(pdfsID_), domainUID(domainUID_)
   {
      OutflowD3Q27Object = std::make_shared< lbm::OutflowD3Q27 >(blocks, pdfsID);
      FixedDensityD3Q27Object = std::make_shared< lbm::FixedDensityD3Q27 >(blocks, pdfsID, density);
      FreeSlipD3Q27Object = std::make_shared< lbm::FreeSlipD3Q27 >(blocks, pdfsID);
      NoSlipD3Q27Object = std::make_shared< lbm::NoSlipD3Q27 >(blocks, pdfsID);
      UBBD3Q27Object = std::make_shared< lbm::UBBD3Q27 >(blocks, pdfsID, u_x, u_y, u_z);
      

      OutflowD3Q27Object->fillFromFlagField<FlagField_T>(blocks, flagID, walberla::FlagUID("Outflow"), domainUID);
      FixedDensityD3Q27Object->fillFromFlagField<FlagField_T>(blocks, flagID, walberla::FlagUID("FixedDensity"), domainUID);
      FreeSlipD3Q27Object->fillFromFlagField<FlagField_T>(blocks, flagID, walberla::FlagUID("FreeSlip"), domainUID);
      NoSlipD3Q27Object->fillFromFlagField<FlagField_T>(blocks, flagID, walberla::FlagUID("NoSlip"), domainUID);
      UBBD3Q27Object->fillFromFlagField<FlagField_T>(blocks, flagID, walberla::FlagUID("UBB"), domainUID);
      
   }

   void run (IBlock * block)
   {
      OutflowD3Q27Object->run(block);
      FixedDensityD3Q27Object->run(block);
      FreeSlipD3Q27Object->run(block);
      NoSlipD3Q27Object->run(block);
      UBBD3Q27Object->run(block);
      
   }

   void inner (IBlock * block)
   {
      OutflowD3Q27Object->inner(block);
      FixedDensityD3Q27Object->inner(block);
      FreeSlipD3Q27Object->inner(block);
      NoSlipD3Q27Object->inner(block);
      UBBD3Q27Object->inner(block);
      
   }

   void outer (IBlock * block)
   {
      OutflowD3Q27Object->outer(block);
      FixedDensityD3Q27Object->outer(block);
      FreeSlipD3Q27Object->outer(block);
      NoSlipD3Q27Object->outer(block);
      UBBD3Q27Object->outer(block);
      
   }

   void operator() (IBlock * block)
   {
      run(block);
   }

   std::function<void (IBlock *)> getSweep(Type type = Type::ALL)
   {
      switch (type)
      {
      case Type::INNER:
         return [this](IBlock* block) { this->inner(block); };
      case Type::OUTER:
         return [this](IBlock* block) { this->outer(block); };
      default:
         return [this](IBlock* block) { this->run(block); };
      }
   }

   weak_ptr< StructuredBlockStorage > blocks_;
   BlockDataID flagID;
   BlockDataID pdfsID;
   walberla::FlagUID domainUID;

   shared_ptr<lbm::OutflowD3Q27> OutflowD3Q27Object;
   shared_ptr<lbm::FixedDensityD3Q27> FixedDensityD3Q27Object;
   shared_ptr<lbm::FreeSlipD3Q27> FreeSlipD3Q27Object;
   shared_ptr<lbm::NoSlipD3Q27> NoSlipD3Q27Object;
   shared_ptr<lbm::UBBD3Q27> UBBD3Q27Object;
   
};

}
}
