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
//! \file contact.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "blockforest/StructuredBlockForest.h"

#include "core/DataTypes.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"

#include "field/FlagField.h"
#include "field/GhostLayerField.h"

#include <set>
#include <vector>

namespace walberla
{
namespace lbm
{
class contact
{
 public:
   struct IndexInfo
   {
      int32_t x1;
      int32_t y1;
      int32_t z1;
      int32_t x2;
      int32_t y2;
      int32_t z2;
      IndexInfo(int32_t x1_, int32_t y1_, int32_t z1_, int32_t x2_, int32_t y2_, int32_t z2_)
         : x1(x1_), y1(y1_), z1(z1_), x2(x2_), y2(y2_), z2(z2_)
      {}
      bool operator==(const IndexInfo& o) const
      {
         return x1 == o.x1 && y1 == o.y1 && z1 == o.z1 && x2 == o.x2 && y2 == o.y2 && z2 == o.z2;
      }
   };

   class IndexVectors
   {
    public:
      using CpuIndexVector = std::vector< IndexInfo >;

      enum Type { ALL = 0, INNER = 1, OUTER = 2, NUM_TYPES = 3 };

      IndexVectors() : cpuVectors_(NUM_TYPES) {}
      bool operator==(IndexVectors& other) { return other.cpuVectors_ == cpuVectors_; }

      CpuIndexVector& indexVector(Type t) { return cpuVectors_[t]; }
      IndexInfo* pointerCpu(Type t) { return &(cpuVectors_[t][0]); }

    private:
      std::vector< CpuIndexVector > cpuVectors_;
   };

   contact(const shared_ptr< StructuredBlockForest >& blocks, BlockDataID phaseFieldID_, double alpha)
      : phaseFieldID(phaseFieldID_), alpha_(alpha)
   {
      auto createIdxVector = [](IBlock* const, StructuredBlockStorage* const) { return new IndexVectors(); };
      indexVectorID = blocks->addStructuredBlockData< IndexVectors >(createIdxVector, "IndexField_contact_angle");
   };

   void operator()(IBlock* block);
   void inner(IBlock* block);
   void outer(IBlock* block);

   template< typename NormalField_T >
   void fillFromNormalField(const shared_ptr< StructuredBlockForest >& blocks, ConstBlockDataID normalFieldID)
   {
      for (auto& block : *blocks)
      {
         auto* indexVectors     = block.getData< IndexVectors >(indexVectorID);
         auto& indexVectorAll   = indexVectors->indexVector(IndexVectors::ALL);
         auto& indexVectorInner = indexVectors->indexVector(IndexVectors::INNER);
         auto& indexVectorOuter = indexVectors->indexVector(IndexVectors::OUTER);

         auto* normalField = block.getData< NormalField_T >(normalFieldID);

         auto inner = normalField->xyzSize();
         inner.expand(cell_idx_t(-1));

         indexVectorAll.clear();
         indexVectorInner.clear();
         indexVectorOuter.clear();
         // clang-format off
         WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(normalField, Cell globalCell;
            blocks->transformBlockLocalToGlobalCell(globalCell, block, Cell(x, y, z));
            if(normalField->get(x, y, z, 0) != 0 || normalField->get(x, y, z, 1) != 0 || normalField->get(x, y, z, 2) != 0)
               {
                  auto element = IndexInfo(x, y,  z, normalField->get(x, y, z, 0), normalField->get(x, y, z, 1), normalField->get(x, y, z, 2));
                  indexVectorAll.push_back( element );
                  if( inner.contains( x, y, z ) )
                      indexVectorInner.push_back( element );
                  else
                     indexVectorOuter.push_back( element );
               }
         )
         // clang-format on
      }
   }

 private:
   void run(IBlock* block, IndexVectors::Type type);

   BlockDataID indexVectorID;

 public:
   BlockDataID phaseFieldID;
   double alpha_;
};

} // namespace lbm
} // namespace walberla