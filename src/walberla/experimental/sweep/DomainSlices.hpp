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
//! \file DomainSlices.hpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "core/cell/CellInterval.h"

#include "domain_decomposition/IBlock.h"

#include "stencil/Directions.h"

#include <memory>

namespace walberla::experimental::sweep
{

template< typename T >
concept Sweep = requires(T obj, IBlock* block) {
   { obj(block) } -> std::same_as< void >;
};

template< typename T >
concept CellIntervalSweep = requires(T obj, IBlock* block, const CellInterval& ci) {
   { obj.runOnCellInterval(block, ci) } -> std::same_as< void >;
};

namespace detail
{
template< stencil::Direction BorderDir >
struct BorderSweepSlice
{
   shared_ptr< StructuredBlockForest > blocks;
   cell_idx_t offset;
   cell_idx_t thickness;
   cell_idx_t expand;

   CellInterval ci();
   bool isBlockAtBorder(IBlock& b);
};

template<>
CellInterval BorderSweepSlice< stencil::Direction::W >::ci()
{
   return { { offset, -expand, -expand },
            { offset + thickness - 1,                                       //
              cell_idx_c(blocks->getNumberOfYCellsPerBlock()) - 1 + expand, //
              cell_idx_c(blocks->getNumberOfZCellsPerBlock()) - 1 + expand } };
};

template<>
CellInterval BorderSweepSlice< stencil::Direction::E >::ci()
{
   return { { cell_idx_c(blocks->getNumberOfXCellsPerBlock()) - offset - thickness, -expand, -expand },
            { cell_idx_c(blocks->getNumberOfXCellsPerBlock()) - (offset + 1), //
              cell_idx_c(blocks->getNumberOfYCellsPerBlock()) - 1 + expand,   //
              cell_idx_c(blocks->getNumberOfZCellsPerBlock()) - 1 + expand } };
};

template<>
CellInterval BorderSweepSlice< stencil::Direction::S >::ci()
{
   return { { -expand, offset, -expand },
            { cell_idx_c(blocks->getNumberOfXCellsPerBlock()) - 1 + expand, //
              offset + thickness - 1,                                       //
              cell_idx_c(blocks->getNumberOfZCellsPerBlock()) - 1 + expand } };
};

template<>
CellInterval BorderSweepSlice< stencil::Direction::N >::ci()
{
   return { { -expand, cell_idx_c(blocks->getNumberOfYCellsPerBlock()) - offset - thickness, -expand },
            { cell_idx_c(blocks->getNumberOfXCellsPerBlock()) - 1 + expand,   //
              cell_idx_c(blocks->getNumberOfYCellsPerBlock()) - (offset + 1), //
              cell_idx_c(blocks->getNumberOfZCellsPerBlock()) - 1 + expand } };
};

template<>
CellInterval BorderSweepSlice< stencil::Direction::B >::ci()
{
   return { { -expand, -expand, offset },
            { cell_idx_c(blocks->getNumberOfXCellsPerBlock()) - 1 + expand, //
              cell_idx_c(blocks->getNumberOfYCellsPerBlock()) - 1 + expand, //
              offset + thickness - 1 } };
};

template<>
CellInterval BorderSweepSlice< stencil::Direction::T >::ci()
{
   return { { -expand, -expand, cell_idx_c(blocks->getNumberOfZCellsPerBlock()) - offset - thickness },
            { cell_idx_c(blocks->getNumberOfXCellsPerBlock()) - 1 + expand, //
              cell_idx_c(blocks->getNumberOfYCellsPerBlock()) - 1 + expand, //
              cell_idx_c(blocks->getNumberOfZCellsPerBlock()) - (offset + 1) } };
};

template<>
bool BorderSweepSlice< stencil::Direction::W >::isBlockAtBorder(IBlock& b)
{
   return blocks->atDomainXMinBorder(b);
}

template<>
bool BorderSweepSlice< stencil::Direction::E >::isBlockAtBorder(IBlock& b)
{
   return blocks->atDomainXMaxBorder(b);
}

template<>
bool BorderSweepSlice< stencil::Direction::S >::isBlockAtBorder(IBlock& b)
{
   return blocks->atDomainYMinBorder(b);
}

template<>
bool BorderSweepSlice< stencil::Direction::N >::isBlockAtBorder(IBlock& b)
{
   return blocks->atDomainYMaxBorder(b);
}

template<>
bool BorderSweepSlice< stencil::Direction::B >::isBlockAtBorder(IBlock& b)
{
   return blocks->atDomainZMinBorder(b);
}

template<>
bool BorderSweepSlice< stencil::Direction::T >::isBlockAtBorder(IBlock& b)
{
   return blocks->atDomainZMaxBorder(b);
}

} // namespace detail

/**
 * @brief Sweep restriction to a domain border.
 * 
 * Wrapper sweep that restricts a sweep to a slice of cells at a domain border.
 * To create border sweeps, use `SweepFactory::atDomainBorder`.
 */
template< stencil::Direction BorderDir, CellIntervalSweep Sweep_T >
class BorderSweep
{
 private:
   Sweep_T sweep_;
   detail::BorderSweepSlice< BorderDir > slice_;
   CellInterval ci_;

 public:
   BorderSweep(const shared_ptr< StructuredBlockForest >& blocks, const Sweep_T& sweep, cell_idx_t offset = 0,
               uint_t thickness = 1, cell_idx_t expand = 0)
      : sweep_{ sweep }, slice_{ blocks, offset, cell_idx_c(thickness), expand }, ci_{ slice_.ci() }
   {}

   void operator()(IBlock* block)
   {
      if (slice_.isBlockAtBorder(*block)) { sweep_.runOnCellInterval(block, ci_); }
   }
};

/**
 * @brief Fusion of two sweeps.
 * 
 * Wrapper sweep that fuses two sub-sweeps, executing them on each block in order.
 * To create a fused sweep, use `SweepFactory::fuse`.
 */
template< Sweep Sw1, Sweep Sw2 >
class FusedSweep
{
 private:
   Sw1 sw1_;
   Sw2 sw2_;

 public:
   FusedSweep(const Sw1& sw1, const Sw2& sw2) : sw1_{ sw1 }, sw2_{ sw2 } {}

   void operator()(IBlock* b)
   {
      sw1_(b);
      sw2_(b);
   }

   template< typename = void >
      requires(CellIntervalSweep< Sw1 > && CellIntervalSweep< Sw2 >)
   void runOnCellInterval(IBlock* b, const CellInterval& ci)
   {
      sw1_.runOnCellInterval(b, ci);
      sw2_.runOnCellInterval(b, ci);
   }
};

namespace detail
{
auto doFuse(const Sweep auto& sw1, const Sweep auto& sw2) { return FusedSweep{ sw1, sw2 }; }

auto doFuse(const Sweep auto& first, const Sweep auto&... others)
{
   return FusedSweep{ first, doFuse(std::forward< decltype(others) >(others)...) };
}
} // namespace detail

/**
 * @brief Factory for derived and wrapper sweeps.
 * 
 * Objects of this class can be used to construct derived sweeps
 * from existing sweep objects.
 * This inclused fusing sweeps,
 * restricting them to certain cell intervals (TODO) or domain borders,
 * or selectively activating them according to certain conditions (TODO).
 */
class SweepFactory
{
 private:
   std::shared_ptr< StructuredBlockForest > blocks_;

 public:
   SweepFactory(const shared_ptr< StructuredBlockForest >& blocks) : blocks_{ blocks } {}

   /**
    * @brief Restrict a sweep to a slice at the simulation domain border.
    * 
    * Return a wrapper sweep which executes the given `sweep`
    * only on cells at the simulation domain border specified by the `Dir` template argument.
    * 
    * @param sweep The original sweep object
    * @param offset Offset of the outermost cell layer of the border slice from the actual border.
    *    `0` means the outermost layer of valid cells,
    *    negative numbers shift the slice into the ghost layer,
    *    while positive numbers shift it into the domain.
    * @param thickness Thickness of the iteration slice; must be a positive number.
    * @param expand Number of cells by which to expand the slice into ghost layers along the border.
    * 
    */
   template< stencil::Direction Dir, CellIntervalSweep Sweep_T >
   BorderSweep< Dir, Sweep_T > atDomainBorder( //
      const Sweep_T& sweep, //
      cell_idx_t offset = 0, //
      uint_t thickness = 1, //
      cell_idx_t expand = 0 //
   )
   {
      WALBERLA_ASSERT_GREATER_EQUAL(thickness, 1);
      return { blocks_, sweep, offset, thickness, expand };
   }

   /**
    * @brief Fuse a sequence of sweeps.
    * 
    * Return a new wrapper sweep which, on each block, executes all of the given sweeps in
    * the specified order.
    */
   auto fuse(const Sweep auto&... sweeps) { return detail::doFuse(std::forward< decltype(sweeps) >(sweeps)...); }
};

} // namespace walberla::experimental::sweep
