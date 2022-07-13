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
//! \file SimplePressureWithFreeSurface.h
//! \ingroup free_surface
//! \author Martin Bauer
//! \author Christian Godenschwager
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief SimplePressure boundary condition for the free surface LBM.
//
//======================================================================================================================

#pragma once

#include "boundary/Boundary.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/math/Vector3.h"

#include "field/FlagField.h"

#include "lbm/field/PdfField.h"
#include "lbm/free_surface/bubble_model/BubbleModel.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"
#include "lbm/lattice_model/ForceModel.h"

#include "stencil/Directions.h"

#include <vector>

namespace walberla
{
namespace free_surface
{
/***********************************************************************************************************************
 * SimplePressure boundary condition for the free surface LBM. The implementation is almost identical to the general
 * lbm/boundary/SimplePressure.h boundary condition, however, the boundary pressure (density) is also set in the bubble
 * model.
 **********************************************************************************************************************/
template< typename LatticeModel_T, typename FlagField_T >
class SimplePressureWithFreeSurface : public boundary::Boundary< typename FlagField_T::flag_t >
{
   using PdfField_T = lbm::PdfField< LatticeModel_T >;
   using Stencil_T  = typename LatticeModel_T::Stencil;
   using flag_t     = typename FlagField_T::flag_t;

 public:
   static const bool threadsafe = true;

   static std::shared_ptr< BoundaryConfiguration > createConfiguration(const Config::BlockHandle& /*config*/)
   {
      return std::make_shared< BoundaryConfiguration >();
   }

   SimplePressureWithFreeSurface(const BoundaryUID& boundaryUID, const FlagUID& uid, IBlock* block,
                                 PdfField_T* const pdfField, FlagField_T* flagField, flag_t interfaceFlag,
                                 const real_t latticeDensity)
      : Boundary< flag_t >(boundaryUID), uid_(uid), block_(block), pdfs_(pdfField), flagField_(flagField),
        interfaceFlag_(interfaceFlag), bubbleModel_(nullptr), latticeDensity_(latticeDensity)
   {
      WALBERLA_ASSERT_NOT_NULLPTR(pdfs_);
      WALBERLA_ASSERT_NOT_NULLPTR(flagField);
   }

   void pushFlags(std::vector< FlagUID >& uids) const { uids.push_back(uid_); }

   void beforeBoundaryTreatment() const {}
   void afterBoundaryTreatment() const {}

   template< typename Buffer_T >
   void packCell(Buffer_T&, const cell_idx_t, const cell_idx_t, const cell_idx_t) const
   {}

   template< typename Buffer_T >
   void registerCell(Buffer_T&, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t)
   {}

   void registerCell(const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t, const BoundaryConfiguration&)
   {}
   void registerCells(const flag_t, const CellInterval&, const BoundaryConfiguration&) const {}
   template< typename CellIterator >
   void registerCells(const flag_t, const CellIterator&, const CellIterator&, const BoundaryConfiguration&) const
   {}

   void unregisterCell(const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t) const {}

   void setLatticeDensity(real_t newLatticeDensity) { latticeDensity_ = newLatticeDensity; }

   void setBubbleModel(BubbleModelBase* bubbleModel) { bubbleModel_ = bubbleModel; }

#ifndef NDEBUG
   inline void treatDirection(const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const stencil::Direction dir,
                              const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask)
#else
   inline void treatDirection(const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const stencil::Direction dir,
                              const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/)
#endif
   {
      WALBERLA_ASSERT_EQUAL(nx, x + cell_idx_c(stencil::cx[dir]));
      WALBERLA_ASSERT_EQUAL(ny, y + cell_idx_c(stencil::cy[dir]));
      WALBERLA_ASSERT_EQUAL(nz, z + cell_idx_c(stencil::cz[dir]));

      WALBERLA_ASSERT_UNEQUAL((mask & this->mask_), numeric_cast< flag_t >(0));
      WALBERLA_ASSERT_EQUAL((mask & this->mask_),
                            this->mask_); // only true if "this->mask_" only contains one single flag, which is the case
                                          // for the current implementation of this boundary condition (SimplePressure)
      Vector3< real_t > u = pdfs_->getVelocity(x, y, z);

      // set density in bubble model according to pressure boundary condition
      if (bubbleModel_ && flagField_->isFlagSet(x, y, z, interfaceFlag_))
      {
         bubbleModel_->setDensity(block_, Cell(x, y, z), latticeDensity_);
      }

      // result will be streamed to (x,y,z, stencil::inverseDir[d]) during sweep
      pdfs_->get(nx, ny, nz, Stencil_T::invDirIdx(dir)) =
         -pdfs_->get(x, y, z, Stencil_T::idx[dir]) // anti-bounce-back
         + real_c(2) * lbm::EquilibriumDistribution< LatticeModel_T >::getSymmetricPart(
                          dir, u, latticeDensity_); // pressure term
   }

 protected:
   FlagUID uid_;

   IBlock* block_;
   PdfField_T* pdfs_;
   FlagField_T* flagField_;
   flag_t interfaceFlag_;
   BubbleModelBase* bubbleModel_;

   real_t latticeDensity_;
};

} // namespace free_surface
} // namespace walberla
