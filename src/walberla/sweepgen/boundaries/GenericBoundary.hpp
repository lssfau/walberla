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
//! \file GenericBoundary.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//! \author Philipp Suffa <philipp.suffa@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "core/DataTypes.h"
#include "core/cell/Cell.h"
#include "core/logging/Logging.h"
#include "domain_decomposition/BlockDataID.h"
#include "field/FlagField.h"
#include "geometry/bodies/BodyOverlapFunctions.h"
#include "stencil/Directions.h"

#include <memory>
#include <tuple>

#include "walberla/experimental/sweep/SparseIndexList.hpp"
#include "walberla/experimental/sweep/Sweeper.hpp"

namespace walberla::sweepgen
{

template< typename Impl >
struct GenericBoundaryFactoryImplTraits
{
   using Stencil_T         = typename Impl::Stencil;
   using IdxStruct_T       = typename Impl::IdxStruct;
   using DataStruct_T      = typename Impl::DataStruct;
   using IndexListMemTag_T = typename Impl::memtag_t;
   using IndexList_T       = walberla::experimental::sweep::SparseIndexList< IdxStruct_T, IndexListMemTag_T >;
};

/**
 * @brief Represents a potential boundary link.
 */
struct PotentialBoundaryLink
{
   const Block & block;
   Cell fluidCell;
   stencil::Direction dir;
   Cell wallCell;
};


template< typename F, typename LinkData >
concept LinkPredicate = requires(F f, PotentialBoundaryLink link) {
        { f(link) } -> std::same_as< std::conditional_t< std::same_as< LinkData, void >, bool, std::optional< LinkData > > >;
};

template< typename F, typename LinkData >
concept ReturnsLinkData = requires(F f, PotentialBoundaryLink link) {
   { f(link) } -> std::same_as< LinkData >;
};

template< typename F >
concept ReturnsRealT = requires(F f, Vector3<real_t> v) {
   { f(v) } -> std::same_as< real_t >;
};

template< typename F >
concept ContainsCallable = requires(F f, Vector3<real_t> v) {
   { walberla::geometry::contains(f,v) } -> std::same_as< bool >;
};


template< typename Impl >
class GenericBoundaryFactory
{
   using ImplTraits = GenericBoundaryFactoryImplTraits< Impl >;

 public:
   GenericBoundaryFactory(const shared_ptr< StructuredBlockForest > blocks)
      : blocks_{ blocks }
   {}

   template< typename FlagField_T >
   auto fromFlagField(BlockDataID flagFieldID, field::FlagUID boundaryFlagUID, field::FlagUID fluidFlagUID)
   {
      using flag_t      = typename FlagField_T::flag_t;
      using IndexList_T = typename ImplTraits::IndexList_T;
      using Stencil_T     = typename ImplTraits::Stencil_T;

      IndexList_T indexList{ *blocks_ };

      for (auto& block : *blocks_)
      {
         FlagField_T* flagField = block.template getData< FlagField_T >(flagFieldID);
         flag_t boundaryFlag{ flagField->getFlag(boundaryFlagUID) };
         flag_t fluidFlag{ flagField->getFlag(fluidFlagUID) };

         if (!(flagField->flagExists(boundaryFlagUID) && flagField->flagExists(fluidFlagUID))) continue;

         auto& idxVector = indexList.getVector(block);

         for (auto it = flagField->beginXYZ(); it != flagField->end(); ++it)
         {
            Cell fluidCell{ it.cell() };
            if (!isFlagSet(it, fluidFlag)) { continue; }

            for (auto dIt = Stencil_T::beginNoCenter(); dIt != Stencil_T::end(); ++dIt)
            {
               stencil::Direction dir{ *dIt };
               Cell neighbor{ fluidCell + dir };
               if (flagField->isFlagSet(neighbor, boundaryFlag)) { idxVector.emplace_back(fluidCell, dir); }
            }
         }
      }

      return impl().irregularFromIndexVector(indexList);
   }

   /**
    * @brief Construct a boundary condition sweep by selecting boundary links via a predicate
    * 
    * Construct a boundary condition sweep by iterating over all lattice links
    * and setting the boundary condition on all links where the given predicate `pred`
    * evaluates to true.
    * 
    * @param pred Predicate to select boundary links.
    */
   template< LinkPredicate<typename ImplTraits::DataStruct_T> P >
   auto fromLinks(P pred)
   {
      using IndexList_T = typename ImplTraits::IndexList_T;
      using Stencil_T     = typename ImplTraits::Stencil_T;

      IndexList_T indexList{ *blocks_ };
      walberla::experimental::sweep::SerialSweeper sweeper{ blocks_ };

      for (auto& ib : *blocks_)
      {
         Block & block = dynamic_cast< Block & >(ib);
         auto& idxVector = indexList.getVector(block);

         sweeper.forAllCells([&](Cell fluidCell) {
            for (auto dIt = Stencil_T::beginNoCenter(); dIt != Stencil_T::end(); ++dIt)
            {
               stencil::Direction dir{ *dIt };
               Cell wallCell{ fluidCell + dir };
               PotentialBoundaryLink link { .block=block, .fluidCell=fluidCell, .dir=dir, .wallCell=wallCell };

               if constexpr ( std::same_as< typename ImplTraits::DataStruct_T, void > ) {
                  if (pred(link)) {
                     idxVector.emplace_back(fluidCell.x(), fluidCell.y(), fluidCell.z(), dir);
                  }
               } else {
                  auto result = pred(link);
                  if(result.has_value()) {
                     idxVector.emplace_back(fluidCell.x(), fluidCell.y(), fluidCell.z(), dir, *result);
                  }
               }
            }
         });
      }
      return impl().irregularFromIndexVector(indexList);
   }


   template <ContainsCallable Body>
   requires(std::same_as<typename ImplTraits::DataStruct_T, void >)
   auto fromBody(Body body)
   {
      using IndexList_T = typename ImplTraits::IndexList_T;
      using Stencil_T     = typename ImplTraits::Stencil_T;

      IndexList_T indexList{ *blocks_ };
      walberla::experimental::sweep::SerialSweeper sweeper{ blocks_ };

      for (auto& ib : *blocks_)
      {
         Block & block = dynamic_cast< Block & >(ib);
         auto& idxVector = indexList.getVector(block);

         sweeper.forAllCells([&](Cell fluidCell) {
            auto fluidCellCenter = blocks_->getGlobalCellCenterFromBlockLocalCell(fluidCell, ib);
            if (!walberla::geometry::contains(body, fluidCellCenter)) {
               for (auto dIt = Stencil_T::beginNoCenter(); dIt != Stencil_T::end(); ++dIt)
               {
                  stencil::Direction dir{ *dIt };
                  Cell wallCell{ fluidCell + dir };
                  auto wallCellCenter = blocks_->getGlobalCellCenterFromBlockLocalCell(wallCell, ib);
                  if (walberla::geometry::contains(body, wallCellCenter)) {
                     idxVector.emplace_back(fluidCell.x(), fluidCell.y(), fluidCell.z(), dir);
                  }
               }
            }
         });
      }
      return impl().irregularFromIndexVector(indexList);
   }


   template < ContainsCallable Body, ReturnsLinkData<typename ImplTraits::DataStruct_T> DataStruct>
   requires(!std::same_as< typename ImplTraits::DataStruct_T, void >)
   auto fromBody(Body body, DataStruct data)
   {
      using IndexList_T = typename ImplTraits::IndexList_T;
      using Stencil_T     = typename ImplTraits::Stencil_T;

      IndexList_T indexList{ *blocks_ };
      walberla::experimental::sweep::SerialSweeper sweeper{ blocks_ };

      for (auto& ib : *blocks_)
      {
         Block & block = dynamic_cast< Block & >(ib);
         auto& idxVector = indexList.getVector(block);

         sweeper.forAllCells([&](Cell fluidCell) {
            auto fluidCellCenter = blocks_->getGlobalCellCenterFromBlockLocalCell(fluidCell, ib);
            if (!walberla::geometry::contains(body, fluidCellCenter)) {
               for (auto dIt = Stencil_T::beginNoCenter(); dIt != Stencil_T::end(); ++dIt)
               {
                  stencil::Direction dir{ *dIt };
                  Cell wallCell{ fluidCell + dir };
                  auto wallCellCenter = blocks_->getGlobalCellCenterFromBlockLocalCell(wallCell, ib);
                  if (walberla::geometry::contains(body, wallCellCenter)) {
                     PotentialBoundaryLink link { .block=block, .fluidCell=fluidCell, .dir=dir, .wallCell=wallCell };
                     idxVector.emplace_back(fluidCell.x(), fluidCell.y(), fluidCell.z(), dir, data(link));
                  }
               }
            }
         });
      }
      return impl().irregularFromIndexVector(indexList);
   }


   template <ReturnsRealT DistFunc>
   requires(std::same_as<typename ImplTraits::DataStruct_T, void >)
   auto fromDistanceFunction(DistFunc dist)
   {
      WALBERLA_LOG_WARNING_ON_ROOT("Creating boundary links from a distance function might be very slow for large meshes. You may use a voxelisation onto a flag field and create the boundary links from flag field instead.")
      WcTimer simTimer;
      simTimer.start();

      using IndexList_T = typename ImplTraits::IndexList_T;
      using Stencil_T     = typename ImplTraits::Stencil_T;

      IndexList_T indexList{ *blocks_ };

      walberla::experimental::sweep::SerialSweeper sweeper{ blocks_ };

      for (auto& ib : *blocks_)
      {
         Block & block = dynamic_cast< Block & >(ib);
         auto& idxVector = indexList.getVector(block);

         sweeper.forAllCells([&](Cell fluidCell) {
            auto fluidCellCenter = blocks_->getGlobalCellCenterFromBlockLocalCell(fluidCell, ib);
            if(dist(fluidCellCenter) > 0.0) {
               for (auto dIt = Stencil_T::beginNoCenter(); dIt != Stencil_T::end(); ++dIt)
               {
                  stencil::Direction dir{ *dIt };
                  Cell wallCell{ fluidCell + dir };
                  auto wallCellCenter = blocks_->getGlobalCellCenterFromBlockLocalCell(wallCell, ib);
                  if (dist(wallCellCenter) <= 0.0) {
                     idxVector.emplace_back(fluidCell.x(), fluidCell.y(), fluidCell.z(), dir);
                  }
               }
            }
         });
      }
      simTimer.end();
      double time = simTimer.max();
      WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(time, walberla::mpi::MAX); }
      WALBERLA_LOG_INFO_ON_ROOT("Created boundary links from distance function in " << time << "s")

      return impl().irregularFromIndexVector(indexList);
   }


   template < ReturnsRealT DistFunc, ReturnsLinkData<typename ImplTraits::DataStruct_T> DataStruct>
   requires(!std::same_as<typename ImplTraits::DataStruct_T, void >)
   auto fromDistanceFunction(DistFunc dist, DataStruct data)
   {
      WALBERLA_LOG_WARNING_ON_ROOT("Creating boundary links from a distance function might be very slow for large meshes. You may use a voxelisation onto a flag field and create the boundary links from flag field instead.")
      WcTimer simTimer;
      simTimer.start();

      using IndexList_T = typename ImplTraits::IndexList_T;
      using Stencil_T     = typename ImplTraits::Stencil_T;

      IndexList_T indexList{ *blocks_ };

      walberla::experimental::sweep::SerialSweeper sweeper{ blocks_ };

      for (auto& ib : *blocks_)
      {
         Block & block = dynamic_cast< Block & >(ib);
         auto& idxVector = indexList.getVector(block);

         sweeper.forAllCells([&](Cell fluidCell) {
            auto fluidCellCenter = blocks_->getGlobalCellCenterFromBlockLocalCell(fluidCell, ib);
            if(dist(fluidCellCenter) > 0.0) {
               for (auto dIt = Stencil_T::beginNoCenter(); dIt != Stencil_T::end(); ++dIt)
               {
                  stencil::Direction dir{ *dIt };
                  Cell wallCell{ fluidCell + dir };
                  auto wallCellCenter = blocks_->getGlobalCellCenterFromBlockLocalCell(wallCell, ib);
                  if (dist(wallCellCenter) <= 0.0) {
                     PotentialBoundaryLink link { .block=block, .fluidCell=fluidCell, .dir=dir, .wallCell=wallCell };
                     idxVector.emplace_back(fluidCell.x(), fluidCell.y(), fluidCell.z(), dir, data(link));
                  }
               }
            }
         });
      }
      simTimer.end();
      double time = simTimer.max();
      WALBERLA_MPI_SECTION() { walberla::mpi::reduceInplace(time, walberla::mpi::MAX); }
      WALBERLA_LOG_INFO_ON_ROOT("Created boundary links from distance function in " << time << "s")

      return impl().irregularFromIndexVector(indexList);
   }

 protected:
   shared_ptr< StructuredBlockForest > blocks_;

 private:
   Impl& impl() { return *static_cast< Impl* >(this); }

   const Impl& impl() const { return *static_cast< const Impl* >(this); }
};

} // namespace walberla::sweepgen
