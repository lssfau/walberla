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
//! \file Reconstructor.h
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================


#pragma once

#include "domain_decomposition/StructuredBlockStorage.h"

#include "mesa_pd/data/IAccessor.h"
#include "mesa_pd/common/ParticleFunctions.h"

#include "lbm/lattice_model/EquilibriumDistribution.h"

#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"

namespace walberla {
namespace lbm_mesapd_coupling {


//**************************************************************************************************************************************
/*!
 *   \brief Classes to be used together with the PDFReconstruction class to reconstruct PDFs
 *
 *   Each reconstructor must exactly implement the member function
 *    template< typename PdfField_T, typename ParticleAccessor_T  >
 *    void operator()( IBlock * const block, const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, PdfField_T * const pdfField,
 *                     const size_t particleIdx, const ParticleAccessor_T & ac);
 *   that reconstructs all PDFs in a specific cell with cell indices x,y,z on a given block.
 *   Additionally, a pointer to the pdfField and information about the formerly present particle is provided (via idx and accessor).
 *
 *   Different variants are available:
 *
 *   - EquilibriumReconstructor:
 *     Determines a local average density and sets the PDFs to their equilibrium values based on this density and the particle's velocity
 *
 *   - EquilibriumAndNonEquilibriumReconstructor:
 *     First reconstructs the equilibrium values with the help of the EquilibriumReconstructor.
 *     Then extrapolates the non-equilibrium part from neighboring cells and adds it to the reconstructed PDFs for better accuracy.
 *     Defaults to EquilibriumReconstructor if no suitable cell for extrapolation is available.
 *
 *   - ExtrapolationReconstructor:
 *     Extrapolates the PDFs of three or two neighboring cells that lie in extrapolation direction.
 *     Optionally, a special treatment after the extrapolation can be used that sets certain moments directly (only D3Q19).
 *     Defaults to EquilibriumAndNonEquilibriumReconstructor if not enough cells in this direction are available.
 *
 *   - GradsMomentReconstructor:
 *     Uses Grad's moment approximation, i.e. equilibrium distribution function extended by pressure gradient information.
 *
 *   EquilibriumAndNonEquilibriumReconstructor and ExtrapolationReconstructor need an extrapolation direction
 *   which is provided by the ExtrapolationDirectionFinder that is to be chosen (see ExtrapolationDirectionFinder.h).
 *
 *  For convenient construction, use the provided makeXReconstructor functions.
 *
 */
//**************************************************************************************************************************************


/*
 * function that determines the number of valid (=fluid) cells for extrapolation in the given extrapolation direction
 * If useDataFromGhostLayers = true, also cells inside the ghostlayer are considered valid so they have to contain valid data by communicating the PDF values beforehand
 * Depending on the scheme, this can be an optimized PDF communication or has to be a full PDF communication.
 */
template< typename BoundaryHandling_T>
uint_t getNumberOfExtrapolationCells( IBlock * const block, const BlockDataID & boundaryHandlingID,
                                      const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z,
                                      const Vector3<cell_idx_t> & extrapolationDirection, const uint_t maximumNumberOfNeededExtrapolationCells,
                                      bool useDataFromGhostLayers )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingID );
   WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );

   if( extrapolationDirection == cell_idx_t(0) ) return uint_t(0);

   CellInterval localDomain = (useDataFromGhostLayers) ? boundaryHandling->getFlagField()->xyzSizeWithGhostLayer() : boundaryHandling->getFlagField()->xyzSize();

   for( auto numCells = uint_t(1); numCells <= maximumNumberOfNeededExtrapolationCells; ++numCells )
   {
      Cell checkCell( x + cell_idx_c(numCells) * extrapolationDirection[0], y + cell_idx_c(numCells) * extrapolationDirection[1], z + cell_idx_c(numCells) * extrapolationDirection[2] );

      // check if cell is inside domain & fluid
      if( !localDomain.contains( checkCell ) )
         return numCells - uint_t(1);

      if( !boundaryHandling->isDomain( checkCell ) )
         return numCells - uint_t(1);
   }
   return maximumNumberOfNeededExtrapolationCells;
}


template< typename BoundaryHandling_T >
class EquilibriumReconstructor
{
public:

   EquilibriumReconstructor( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & boundaryHandlingID,
                             bool useDataFromGhostLayers = false)
   : blockStorage_( blockStorage ), boundaryHandlingID_( boundaryHandlingID ), useDataFromGhostLayers_( useDataFromGhostLayers )
   {}


   template< typename PdfField_T, typename ParticleAccessor_T  >
   void operator()( IBlock * const block, const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, PdfField_T * const pdfField,
                    const size_t particleIdx, const ParticleAccessor_T & ac)
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      WALBERLA_ASSERT_NOT_NULLPTR( block );
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
      WALBERLA_ASSERT_UNEQUAL(particleIdx, ac.getInvalidIdx());

      setLocalEquilibrium(block, x, y, z, pdfField, particleIdx, ac);
   }

private:

   template< typename PdfField_T, typename ParticleAccessor_T  >
   void setLocalEquilibrium( IBlock * const block, const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, PdfField_T * const pdfField,
                             const size_t particleIdx, const ParticleAccessor_T & ac )
   {
      const real_t averageDensity = getLocalAverageDensity( block, x, y, z, pdfField );

      real_t cx, cy, cz;
      blockStorage_->getBlockLocalCellCenter( *block, Cell(x,y,z), cx, cy, cz );

      const auto velocity = mesa_pd::getVelocityAtWFPoint(particleIdx, ac, Vector3<real_t>(cx,cy,cz));

      pdfField->setToEquilibrium( x, y, z, velocity, averageDensity );
   }

   template< typename PdfField_T>
   real_t getLocalAverageDensity( IBlock * const block, const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, PdfField_T * const pdfField ) const
   {
      BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );

      WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );

      CellInterval localDomain = useDataFromGhostLayers_ ? pdfField->xyzSizeWithGhostLayer() : pdfField->xyzSize();

      auto nAverage = uint_t(0);
      auto averageDensity = real_t(0);
      for( auto neighborDir = stencil::D3Q27::beginNoCenter(); neighborDir != stencil::D3Q27::end(); ++neighborDir )
      {
         Cell neighbor( x + neighborDir.cx(), y + neighborDir.cy(), z + neighborDir.cz() );

         // check if neighbor cell is inside domain
         if( !( localDomain.contains( neighbor ) ) )
            continue;

         // check if neighbor cell is a fluid cell
         if( boundaryHandling->isDomain( neighbor ) )
         {
            // obtain density and add it up
            averageDensity += pdfField->getDensity( neighbor );
            ++nAverage;
         }
      }

      return ( nAverage > uint_t( 0 ) ) ? averageDensity / real_c( nAverage ) : real_t(1.0);
   }


   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID boundaryHandlingID_;
   const bool useDataFromGhostLayers_;

};


template< typename BoundaryHandling_T, typename ExtrapolationDirectionFinder_T >
class EquilibriumAndNonEquilibriumReconstructor
{
public:

   EquilibriumAndNonEquilibriumReconstructor( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & boundaryHandlingID,
                                              const shared_ptr<ExtrapolationDirectionFinder_T> & extrapolationDirectionFinder,
                                              uint_t maximumNumberOfExtrapolationCells = uint_t(3),
                                              bool useDataFromGhostLayers = false)
   : blockStorage_( blockStorage ), boundaryHandlingID_( boundaryHandlingID ),
     extrapolationDirectionFinder_( extrapolationDirectionFinder ),
     maximumNumberOfExtrapolationCells_( maximumNumberOfExtrapolationCells ), useDataFromGhostLayers_( useDataFromGhostLayers ),
     equilibriumReconstructor_( EquilibriumReconstructor< BoundaryHandling_T >( blockStorage, boundaryHandlingID, useDataFromGhostLayers ) )
   {
      WALBERLA_ASSERT_LESS_EQUAL(maximumNumberOfExtrapolationCells, uint_t(3), "Only supports up to quadratic extrapolation!");
   }

   template< typename PdfField_T, typename ParticleAccessor_T >
   void operator()( IBlock * const block, const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, PdfField_T * const pdfField,
                    const size_t particleIdx, const ParticleAccessor_T & ac)
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      WALBERLA_ASSERT_NOT_NULLPTR( block );
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
      WALBERLA_ASSERT_UNEQUAL(particleIdx, ac.getInvalidIdx());

      Vector3<cell_idx_t> extrapolationDirection = (*extrapolationDirectionFinder_)( block, x, y, z, particleIdx, ac );

      (*this)(block, x, y, z, pdfField, particleIdx, ac, extrapolationDirection);
   }

   template< typename PdfField_T, typename ParticleAccessor_T >
   void operator()( IBlock * const block, const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, PdfField_T * const pdfField,
                    const size_t particleIdx, const ParticleAccessor_T & ac, const Vector3<cell_idx_t> & extrapolationDirection )
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      WALBERLA_ASSERT_NOT_NULLPTR( block );
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
      WALBERLA_ASSERT_UNEQUAL(particleIdx, ac.getInvalidIdx());

      equilibriumReconstructor_( block, x, y, z, pdfField, particleIdx, ac );

      uint_t numberOfCellsForExtrapolation = getNumberOfExtrapolationCells<BoundaryHandling_T>( block, boundaryHandlingID_, x, y, z, extrapolationDirection, maximumNumberOfExtrapolationCells_, useDataFromGhostLayers_ );

      if( enoughCellsForExtrapolation( numberOfCellsForExtrapolation ) )
      {
         extrapolateNonEquilibrium( x, y, z, pdfField, extrapolationDirection, numberOfCellsForExtrapolation );
      }

   }
private:

   template< typename PdfField_T >
   void extrapolateNonEquilibrium( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, PdfField_T * const pdfField,
                                   const Vector3<cell_idx_t> & extrapolationDirection, const uint_t & numberOfCellsForExtrapolation)
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField );

      WALBERLA_ASSERT_LESS_EQUAL(numberOfCellsForExtrapolation, maximumNumberOfExtrapolationCells_);
      WALBERLA_ASSERT_GREATER(numberOfCellsForExtrapolation, uint_t(0), "Requires at least one point for extrapolation!");

      if( enoughCellsForQuadraticExtrapolation( numberOfCellsForExtrapolation ) )
      {
         applyQuadraticExtrapolation(x, y, z, pdfField, extrapolationDirection);
      }
      else if ( enoughCellsForLinearExtrapolation( numberOfCellsForExtrapolation ) )
      {
         applyLinearExtrapolation(x, y, z, pdfField, extrapolationDirection);
      }
      else
      {
         applyConstantExtrapolation(x, y, z, pdfField, extrapolationDirection);
      }
   }

   bool enoughCellsForQuadraticExtrapolation(uint_t numberOfCellsForExtrapolation )
   {
      return numberOfCellsForExtrapolation >= uint_t(3);
   }

   bool enoughCellsForLinearExtrapolation(uint_t numberOfCellsForExtrapolation )
   {
      return numberOfCellsForExtrapolation >= uint_t(2);
   }

   bool enoughCellsForExtrapolation(uint_t numberOfCellsForExtrapolation )
   {
      return numberOfCellsForExtrapolation >= uint_t(1);
   }

   template< typename PdfField_T >
   void applyQuadraticExtrapolation( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, PdfField_T * const pdfField,
                                     const Vector3<cell_idx_t> & extrapolationDirection)
   {
      WALBERLA_ASSERT_NOT_NULLPTR(pdfField);

      using LatticeModel_T = typename PdfField_T::LatticeModel;

      auto pdfsXf   = getNonEquilibriumPdfsInCell(x +   extrapolationDirection[0], y +   extrapolationDirection[1], z +   extrapolationDirection[2], pdfField);
      auto pdfsXff  = getNonEquilibriumPdfsInCell(x + 2*extrapolationDirection[0], y + 2*extrapolationDirection[1], z + 2*extrapolationDirection[2], pdfField);
      auto pdfsXfff = getNonEquilibriumPdfsInCell(x + 3*extrapolationDirection[0], y + 3*extrapolationDirection[1], z + 3*extrapolationDirection[2], pdfField);

      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         pdfField->get( x, y, z, d.toIdx() ) += real_t(3) * pdfsXf[d.toIdx()] - real_t(3) * pdfsXff[d.toIdx()] + pdfsXfff[d.toIdx()];
      }
   }

   template< typename PdfField_T >
   void applyLinearExtrapolation( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, PdfField_T * const pdfField,
                                  const Vector3<cell_idx_t> & extrapolationDirection)
   {
      WALBERLA_ASSERT_NOT_NULLPTR(pdfField);

      using LatticeModel_T = typename PdfField_T::LatticeModel;

      auto pdfsXf  = getNonEquilibriumPdfsInCell(x +   extrapolationDirection[0], y +   extrapolationDirection[1], z +   extrapolationDirection[2], pdfField);
      auto pdfsXff = getNonEquilibriumPdfsInCell(x + 2*extrapolationDirection[0], y + 2*extrapolationDirection[1], z + 2*extrapolationDirection[2], pdfField);

      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         pdfField->get( x, y, z, d.toIdx() ) += real_t(2) * pdfsXf[d.toIdx()] - pdfsXff[d.toIdx()];
      }
   }

   template< typename PdfField_T >
   void applyConstantExtrapolation( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, PdfField_T * const pdfField,
                                    const Vector3<cell_idx_t> & extrapolationDirection)
   {
      WALBERLA_ASSERT_NOT_NULLPTR(pdfField);

      using LatticeModel_T = typename PdfField_T::LatticeModel;

      // = copy from neighbor cell
      auto pdfsXf = getNonEquilibriumPdfsInCell(x + extrapolationDirection[0], y + extrapolationDirection[1], z + extrapolationDirection[2], pdfField);

      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         pdfField->get( x, y, z, d.toIdx() ) += pdfsXf[d.toIdx()];
      }
   }

   template< typename PdfField_T >
   std::vector<real_t> getNonEquilibriumPdfsInCell( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z,  PdfField_T * const pdfField ) const
   {
      using LatticeModel_T = typename PdfField_T::LatticeModel;

      std::vector< real_t > nonEquilibriumPartOfPdfs(LatticeModel_T::Stencil::Size);

      Vector3< real_t > velocity;
      real_t density = pdfField->getDensityAndVelocity( velocity, x,y,z );

      for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
      {
         nonEquilibriumPartOfPdfs[d.toIdx()] = pdfField->get( x, y, z, d.toIdx() ) - lbm::EquilibriumDistribution< LatticeModel_T >::get( *d, velocity, density );
      }
      return nonEquilibriumPartOfPdfs;
   }

   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID boundaryHandlingID_;
   shared_ptr<ExtrapolationDirectionFinder_T> extrapolationDirectionFinder_;
   const uint_t maximumNumberOfExtrapolationCells_;
   const bool useDataFromGhostLayers_;
   EquilibriumReconstructor< BoundaryHandling_T > equilibriumReconstructor_;

};

template< typename BoundaryHandling_T, typename ExtrapolationDirectionFinder_T, bool EnforceNoSlipConstraintAfterExtrapolation = false >
class ExtrapolationReconstructor
{
public:

   ExtrapolationReconstructor( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & boundaryHandlingID,
                               const shared_ptr<ExtrapolationDirectionFinder_T> & extrapolationDirectionFinder,
                               uint_t maximumNumberOfExtrapolationCells = uint_t(3),
                               bool useDataFromGhostLayers = false)
   : blockStorage_( blockStorage ), boundaryHandlingID_( boundaryHandlingID ),
     extrapolationDirectionFinder_( extrapolationDirectionFinder ),
     maximumNumberOfExtrapolationCells_(maximumNumberOfExtrapolationCells), useDataFromGhostLayers_(useDataFromGhostLayers),
     alternativeReconstructor_( EquilibriumAndNonEquilibriumReconstructor< BoundaryHandling_T, ExtrapolationDirectionFinder_T >
      ( blockStorage, boundaryHandlingID, extrapolationDirectionFinder, uint_t(1), useDataFromGhostLayers ) )
   {
      WALBERLA_ASSERT_LESS_EQUAL(maximumNumberOfExtrapolationCells, uint_t(3), "Only supports up to quadratic extrapolation!");
   }

   template< typename PdfField_T, typename ParticleAccessor_T >
   void operator()( IBlock * const block, const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, PdfField_T * const pdfField,
                    const size_t particleIdx, const ParticleAccessor_T & ac)
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      WALBERLA_ASSERT_NOT_NULLPTR(block);
      WALBERLA_ASSERT_NOT_NULLPTR(pdfField);
      WALBERLA_ASSERT_UNEQUAL(particleIdx, ac.getInvalidIdx());

      Vector3<cell_idx_t> extrapolationDirection = (*extrapolationDirectionFinder_)(block, x, y, z, particleIdx, ac);

      uint_t numberOfCellsForExtrapolation = getNumberOfExtrapolationCells<BoundaryHandling_T>( block, boundaryHandlingID_, x, y, z, extrapolationDirection,
                                                                                                maximumNumberOfExtrapolationCells_, useDataFromGhostLayers_ );

      if( enoughCellsForExtrapolation(numberOfCellsForExtrapolation) )
      {
         extrapolatePDFs( x, y, z, pdfField, extrapolationDirection, numberOfCellsForExtrapolation );
         if( EnforceNoSlipConstraintAfterExtrapolation )
         {
            real_t cx, cy, cz;
            blockStorage_->getBlockLocalCellCenter( *block, Cell(x,y,z), cx, cy, cz );
            const auto localParticleVelocity = mesa_pd::getVelocityAtWFPoint(particleIdx, ac, Vector3<real_t>(cx,cy,cz));

            enforceNoSlipConstraint( x, y, z, pdfField, localParticleVelocity );
         }
      }
      else
      {
         alternativeReconstructor_( block, x, y, z, pdfField, particleIdx, ac, extrapolationDirection );
      }

   }

private:

   bool enoughCellsForExtrapolation( uint_t numberOfCellsForExtrapolation )
   {
      return numberOfCellsForExtrapolation >= uint_t(2);
   }

   template< typename PdfField_T >
   void extrapolatePDFs( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, PdfField_T * const pdfField,
                         const Vector3<cell_idx_t> & extrapolationDirection, const uint_t & numberOfCellsForExtrapolation)
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField );

      using LatticeModel_T = typename PdfField_T::LatticeModel;

      if( numberOfCellsForExtrapolation == uint_t(3) )
      {
         // quadratic normal extrapolation
         for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
         {
            pdfField->get( x, y, z, d.toIdx() ) =   real_t(3) * pdfField->get( x +   extrapolationDirection[0], y +   extrapolationDirection[1], z +   extrapolationDirection[2], d.toIdx() )
                                                  - real_t(3) * pdfField->get( x + 2*extrapolationDirection[0], y + 2*extrapolationDirection[1], z + 2*extrapolationDirection[2], d.toIdx() )
                                                  +             pdfField->get( x + 3*extrapolationDirection[0], y + 3*extrapolationDirection[1], z + 3*extrapolationDirection[2], d.toIdx() );
         }
      } else { // numberOfCellsForExtrapolation == 2
         // linear normal extrapolation
         for( auto d = LatticeModel_T::Stencil::begin(); d != LatticeModel_T::Stencil::end(); ++d )
         {
            pdfField->get( x, y, z, d.toIdx() ) =   real_t(2) * pdfField->get( x +   extrapolationDirection[0], y +   extrapolationDirection[1], z +   extrapolationDirection[2], d.toIdx() )
                                                  -             pdfField->get( x + 2*extrapolationDirection[0], y + 2*extrapolationDirection[1], z + 2*extrapolationDirection[2], d.toIdx() );
         }
      }
   }

   template< typename PdfField_T >
   void enforceNoSlipConstraint( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, PdfField_T * const pdfField, const Vector3<real_t> & localParticleVelocity)
   {
      using LatticeModel_T = typename PdfField_T::LatticeModel;

      static_assert(std::is_same<typename LatticeModel_T::Stencil, stencil::D3Q19>::value || !EnforceNoSlipConstraintAfterExtrapolation, "Enforcing no-slip constraint after extrapolation currently only works with D3Q19 stencil!");

      WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
      WALBERLA_ASSERT( !math::isnan(localParticleVelocity) );

      // transforms to moment space (see MRT collision model) to set the particle's velocity in cell without affecting other moments
      const real_t _1_2  = real_t(1) / real_t(2);
      const real_t _1_3  = real_t(1) / real_t(3);
      const real_t _1_4  = real_t(1) / real_t(4);
      const real_t _1_6  = real_t(1) / real_t(6);
      const real_t _1_8  = real_t(1) / real_t(8);
      const real_t _1_12 = real_t(1) / real_t(12);
      const real_t _1_16 = real_t(1) / real_t(16);
      const real_t _1_18 = real_t(1) / real_t(18);
      const real_t _1_24 = real_t(1) / real_t(24);
      const real_t _1_36 = real_t(1) / real_t(36);
      const real_t _1_48 = real_t(1) / real_t(48);
      const real_t _1_72 = real_t(1) / real_t(72);

      // restriction to D3Q19!
      const real_t vC  = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::C]  );
      const real_t vN  = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::N]  );
      const real_t vS  = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::S]  );
      const real_t vW  = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::W]  );
      const real_t vE  = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::E]  );
      const real_t vT  = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::T]  );
      const real_t vB  = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::B]  );
      const real_t vNW = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::NW] );
      const real_t vNE = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::NE] );
      const real_t vSW = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::SW] );
      const real_t vSE = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::SE] );
      const real_t vTN = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TN] );
      const real_t vTS = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TS] );
      const real_t vTW = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TW] );
      const real_t vTE = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TE] );
      const real_t vBN = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BN] );
      const real_t vBS = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BS] );
      const real_t vBW = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BW] );
      const real_t vBE = pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BE] );

      // transform to moment space and change momentum to particle's velocity ( * rho_0 )
      const real_t m0  = vC + vN + vS + vW + vE + vT + vB + vNW + vNE + vSW + vSE + vTN + vTS + vTW + vTE + vBN + vBS + vBW + vBE;
      const real_t m1  = -vC  + vNW + vNE + vSW + vSE + vTN + vTS + vTW + vTE + vBN + vBS + vBW + vBE;
      const real_t m2  = vC - real_t(2) * ( vN + vS + vW + vE + vT + vB ) + vNW + vNE + vSW + vSE + vTN + vTS + vTW + vTE + vBN + vBS + vBW + vBE;
      const real_t m3  = localParticleVelocity[0];
      const real_t m4  = real_t(2) * vW - real_t(2) * vE - vNW + vNE - vSW + vSE - vTW + vTE - vBW + vBE;
      const real_t m5  = localParticleVelocity[1];
      const real_t m6  = real_t(-2) * vN + real_t(2) * vS + vNW + vNE - vSW - vSE + vTN - vTS + vBN - vBS;
      const real_t m7  = localParticleVelocity[2];
      const real_t m8  = real_t(-2) * vT + real_t(2) * vB + vTN + vTS + vTW + vTE - vBN - vBS - vBW - vBE;
      const real_t m9  = -vN - vS + real_t(2) * vW + real_t(2) * vE - vT - vB + vNW + vNE + vSW + vSE - real_t(2) * vTN
                         - real_t(2) * vTS + vTW + vTE - real_t(2) * vBN - real_t(2) * vBS + vBW + vBE;
      const real_t m10 = vN + vS - real_t(2) * vW - real_t(2) * vE + vT + vB + vNW + vNE + vSW + vSE - real_t(2) * vTN
                         - real_t(2) * vTS + vTW + vTE - real_t(2) * vBN - real_t(2) * vBS + vBW + vBE;
      const real_t m11 = vN  + vS  - vT  - vB  + vNW + vNE + vSW + vSE - vTW - vTE - vBW - vBE;
      const real_t m12 = -vN - vS  + vT  + vB  + vNW + vNE + vSW + vSE - vTW - vTE - vBW - vBE;
      const real_t m13 = -vNW + vNE + vSW - vSE;
      const real_t m14 =  vTN - vTS - vBN + vBS;
      const real_t m15 = -vTW + vTE + vBW - vBE;
      const real_t m16 = -vNW + vNE - vSW + vSE + vTW - vTE + vBW - vBE;
      const real_t m17 = -vNW - vNE + vSW + vSE + vTN - vTS + vBN - vBS;
      const real_t m18 = -vTN - vTS + vTW + vTE + vBN + vBS - vBW - vBE;

      // transform back to PDF space
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::C] )  = _1_3  * m0  - _1_2  * m1  + _1_6  * m2;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::N] )  = _1_18 * m0  - _1_18 * m2  + _1_6  * m5  - _1_6  * m6  - _1_24 * m9  +
                                                                            _1_24 * m10 + _1_8  * m11 - _1_8  * m12;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::S] )  = _1_18 * m0  - _1_18 * m2  - _1_6  * m5  + _1_6  * m6  - _1_24 * m9  +
                                                                            _1_24 * m10 + _1_8  * m11 - _1_8  * m12;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::W] )  = _1_18 * m0  - _1_18 * m2  - _1_6  * m3  + _1_6  * m4  + _1_12 * m9  - _1_12 * m10;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::E] )  = _1_18 * m0  - _1_18 * m2  + _1_6  * m3  - _1_6  * m4  + _1_12 * m9  - _1_12 * m10;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::T] )  = _1_18 * m0  - _1_18 * m2  + _1_6  * m7  - _1_6  * m8  - _1_24 * m9  +
                                                                            _1_24 * m10 - _1_8  * m11 + _1_8  * m12;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::B] )  = _1_18 * m0  - _1_18 * m2  - _1_6  * m7  + _1_6  * m8  - _1_24 * m9  +
                                                                            _1_24 * m10 - _1_8  * m11 + _1_8  * m12;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::NW] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  - _1_12 * m3  - _1_24 * m4  +
                                                                            _1_12 * m5  + _1_24 * m6  + _1_48 * m9  + _1_48 * m10 + _1_16 * m11 +
                                                                            _1_16 * m12 - _1_4  * m13 - _1_8  * m16 - _1_8  * m17;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::NE] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  + _1_12 * m3  + _1_24 * m4  +
                                                                            _1_12 * m5  + _1_24 * m6  + _1_48 * m9  + _1_48 * m10 + _1_16 * m11 +
                                                                            _1_16 * m12 + _1_4  * m13 + _1_8  * m16 - _1_8  * m17;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::SW] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  - _1_12 * m3  - _1_24 * m4  -
                                                                            _1_12 * m5  - _1_24 * m6  + _1_48 * m9  + _1_48 * m10 + _1_16 * m11 +
                                                                            _1_16 * m12 + _1_4  * m13 - _1_8  * m16 + _1_8  * m17;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::SE] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  + _1_12 * m3  + _1_24 * m4  -
                                                                            _1_12 * m5  - _1_24 * m6  + _1_48 * m9  + _1_48 * m10 + _1_16 * m11 +
                                                                            _1_16 * m12 - _1_4  * m13 + _1_8  * m16 + _1_8  * m17;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TN] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  + _1_12 * m5  + _1_24 * m6  +
                                                                            _1_12 * m7  + _1_24 * m8  - _1_24 * m9  - _1_24 * m10 + _1_4  * m14 +
                                                                            _1_8  * m17 - _1_8  * m18;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TS] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  - _1_12 * m5  - _1_24 * m6  +
                                                                            _1_12 * m7  + _1_24 * m8  - _1_24 * m9  - _1_24 * m10 - _1_4  * m14 -
                                                                            _1_8  * m17 - _1_8  * m18;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TW] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  - _1_12 * m3  - _1_24 * m4  +
                                                                            _1_12 * m7  + _1_24 * m8  + _1_48 * m9  + _1_48 * m10 - _1_16 * m11 -
                                                                            _1_16 * m12 - _1_4  * m15 + _1_8  * m16 + _1_8  * m18;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::TE] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  + _1_12 * m3  + _1_24 * m4  +
                                                                            _1_12 * m7  + _1_24 * m8  + _1_48 * m9  + _1_48 * m10 - _1_16 * m11 -
                                                                            _1_16 * m12 + _1_4  * m15 - _1_8  * m16 + _1_8  * m18;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BN] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  + _1_12 * m5  + _1_24 * m6  -
                                                                            _1_12 * m7  - _1_24 * m8  - _1_24 * m9  - _1_24 * m10 - _1_4  * m14 +
                                                                            _1_8  * m17 + _1_8  * m18;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BS] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  - _1_12 * m5  - _1_24 * m6  -
                                                                            _1_12 * m7  - _1_24 * m8  - _1_24 * m9  - _1_24 * m10 + _1_4  * m14 -
                                                                            _1_8  * m17 + _1_8  * m18;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BW] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  - _1_12 * m3  - _1_24 * m4  -
                                                                            _1_12 * m7  - _1_24 * m8  + _1_48 * m9  + _1_48 * m10 - _1_16 * m11 -
                                                                            _1_16 * m12 + _1_4  * m15 + _1_8  * m16 - _1_8  * m18;
      pdfField->get( x, y, z, LatticeModel_T::Stencil::idx[stencil::BE] ) = _1_36 * m0  + _1_24 * m1  + _1_72 * m2  + _1_12 * m3  + _1_24 * m4  -
                                                                            _1_12 * m7  - _1_24 * m8  + _1_48 * m9  + _1_48 * m10 - _1_16 * m11 -
                                                                            _1_16 * m12 - _1_4  * m15 - _1_8  * m16 - _1_8  * m18;
   }

   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID boundaryHandlingID_;
   shared_ptr<ExtrapolationDirectionFinder_T> extrapolationDirectionFinder_;
   const uint_t maximumNumberOfExtrapolationCells_;
   const bool useDataFromGhostLayers_;

   EquilibriumAndNonEquilibriumReconstructor< BoundaryHandling_T, ExtrapolationDirectionFinder_T > alternativeReconstructor_;
};


/*
 * similar to Chikatamarla, Ansumali and Karlin - Grad's approximation for missing data in lattice Boltzmann simulations, EPL (2006)
 * also in Dorschner, Chikatamarla, Boesch, Karlin - Grad's approximation for moving and stationary walls in entropic lattice Boltzmann simulations, Journal of Computational Physics, 2015
 * omegaShear: relaxation RATE that determines the kinematic viscosity
 *
 * To obtain the pressure gradient information, finite differences with central differences (useCentralDifferences) are used if enough information available, else upwinding differences are applied.
 * When useDataFromGhostLayers = true, a full ghost layer sync is required just before the reconstruction step. This is required to avoid inconsistencies when using parallelization as else the behavior close to block boarders is altered.
 */
template< typename BoundaryHandling_T >
class GradsMomentApproximationReconstructor
{
public:

   GradsMomentApproximationReconstructor( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & boundaryHandlingID, real_t omegaShear,
                                          bool recomputeTargetDensity = false, bool useCentralDifferences = true, bool useDataFromGhostLayers = false )
         : blockStorage_( blockStorage ), boundaryHandlingID_( boundaryHandlingID ), omegaShear_(omegaShear),
           recomputeTargetDensity_(recomputeTargetDensity), useCentralDifferences_(useCentralDifferences), useDataFromGhostLayers_(useDataFromGhostLayers)
   {}


   template< typename PdfField_T, typename ParticleAccessor_T  >
   void operator()( IBlock * const block, const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, PdfField_T * const pdfField,
                    const size_t particleIdx, const ParticleAccessor_T & ac)
   {
      static_assert(std::is_base_of<mesa_pd::data::IAccessor, ParticleAccessor_T>::value, "Provide a valid accessor as template");

      WALBERLA_ASSERT_NOT_NULLPTR( block );
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
      WALBERLA_ASSERT_UNEQUAL(particleIdx, ac.getInvalidIdx());

      using LatticeModel = typename PdfField_T::LatticeModel;
      using Stencil = typename LatticeModel::Stencil;

      BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );
      WALBERLA_ASSERT_NOT_NULLPTR( boundaryHandling );

      CellInterval localDomain = (useDataFromGhostLayers_) ? pdfField->xyzSizeWithGhostLayer() : pdfField->xyzSize();

      std::vector<bool> availableStencilIndices(Stencil::Size, false);

      auto nAverage = uint_t(0);
      auto averageDensity = real_t(0);

      // density and velocity used in the reconstruction
      auto targetDensity = real_t(0);
      Vector3<real_t> targetVelocity;

      real_t cx, cy, cz;
      blockStorage_->getBlockLocalCellCenter( *block, Cell(x,y,z), cx, cy, cz );

      // 1. evaluate local average density and find available (=fluid) stencil directions
      for( auto neighborDir = Stencil::beginNoCenter(); neighborDir != Stencil::end(); ++neighborDir )
      {
         Cell neighbor( x + neighborDir.cx(), y + neighborDir.cy(), z + neighborDir.cz() );

         // check if neighbor cell is inside domain
         if( !( localDomain.contains( neighbor ) ) )
            continue;

         // check if neighbor cell is a fluid cell
         if( boundaryHandling->isDomain( neighbor ) )
         {
            // obtain density and add it up
            averageDensity += pdfField->getDensity( neighbor );
            ++nAverage;
            availableStencilIndices[neighborDir.toIdx()] = true;
         }
      }
      averageDensity = ( nAverage > uint_t( 0 ) ) ? averageDensity / real_c( nAverage ) : real_t(1.0);

      // 2. evaluate target velocity
      // we simply use the body velocity in the cell center here since the cell center is probably not far from the particle surface so this is a valid approximation
      // alternatively, some interpolated velocity might be used here as well (see Chikatamarla et al)
      targetVelocity = mesa_pd::getVelocityAtWFPoint(particleIdx, ac, Vector3<real_t>(cx,cy,cz));

      // 3. compute target density (and adapt to compressible/incompressible)

      // note: for compressible lattice models, the targetDensity is centered around 1
      // for incompressible, the targetDensity is just the deviation from the average (=1)

      if(recomputeTargetDensity_)
      {
         // in this variant the target density gets explicitly computed by considering the surrounding PDFs and using "stream" and "velocity bounce back" considerations like in Chikatamarla et al
         for( auto q = Stencil::begin(); q != Stencil::end(); ++q )
         {
            if( availableStencilIndices[q.toInvIdx()])
            {
               // "stream"
               auto invDir = q.inverseDir();
               Cell neighbor( x + stencil::cx[invDir], y + stencil::cy[invDir], z + stencil::cz[invDir] );
               targetDensity += pdfField->get(neighbor, *q);
            }
            else if(availableStencilIndices[q.toIdx()])
            {
               // "velocity bounce back"
               Cell neighbor( x + stencil::cx[*q], y + stencil::cy[*q], z + stencil::cz[*q] );
               targetDensity += pdfField->get(neighbor, q.inverseDir()); // "bounce back part"
               if(LatticeModel::compressible)
               {
                  targetDensity += real_t(6) * averageDensity * LatticeModel::w[ Stencil::idx[*q] ] * ( real_c( stencil::cx[ *q ] ) * targetVelocity[0] +
                                                                                                        real_c( stencil::cy[ *q ] ) * targetVelocity[1] +
                                                                                                        real_c( stencil::cz[ *q ] ) * targetVelocity[2] ); //TODO use wall velocity here?
               } else {
                  targetDensity += real_t(6) * LatticeModel::w[ Stencil::idx[*q] ] * ( real_c( stencil::cx[ *q ] ) * targetVelocity[0] +
                                                                                       real_c( stencil::cy[ *q ] ) * targetVelocity[1] +
                                                                                       real_c( stencil::cz[ *q ] ) * targetVelocity[2] ); //TODO use wall velocity here?
               }

            } else{
               // no neighboring information available (e.g. for center), so we use feq based on target velocity and average density
               targetDensity += lbm::EquilibriumDistribution< LatticeModel >::get( *q, targetVelocity, averageDensity );
            }
         }
      } else {
         // alternatively, we can simply use the average density from the surrounding fluid cells
         // only minor differences have been seen in comparison to recomputing
         targetDensity = (LatticeModel::compressible) ? averageDensity : averageDensity - real_t(1);
      }

      // 4. compute pressure tensor from eq and neq parts
      // 4.1 evaluate Peq -> will be done directly lateron

      // 4.2 evaluate Pneq -> needs velocity gradient tensor
      // grad(u) =
      // | du1/dx1 du2/dx1 du3/dx1 |   | 0 1 2 |   | 0,0  0,1  0,2 |
      // | du1/dx2 du2/dx2 du3/dx2 | = | 3 4 5 | = | 1,0  1,1  1,2 |
      // | du1/dx3 du2/dx3 du3/dx3 |   | 6 7 8 |   | 2,0  2,1  2,2 |
      // evaluation of velocity gradients done via central finite differences if two neighbors available
      // else first-order FD upwinding is carried out if upwinding neighbor available
      // else first-order finite differences if available if other neighbor is available
      // else we assume gradient of 0 (if no fluid neighbors available in respective direction)

      Matrix3<real_t> velocityGradient(real_t(0));
      // check if both neighbors are available for central differences
      if(useCentralDifferences_ && availableStencilIndices[Stencil::idx[stencil::E]] && availableStencilIndices[Stencil::idx[stencil::W]])
      {
         auto neighborVelocity1 = pdfField->getVelocity(Cell(x,y,z)+stencil::E);
         auto neighborVelocity2 = pdfField->getVelocity(Cell(x,y,z)+stencil::W);
         velocityGradient[0] = real_t(0.5) * ( neighborVelocity1[0] - neighborVelocity2[0]); // assuming dx = 1
         velocityGradient[1] = real_t(0.5) * ( neighborVelocity1[1] - neighborVelocity2[1]); // assuming dx = 1
         velocityGradient[2] = real_t(0.5) * ( neighborVelocity1[2] - neighborVelocity2[2]); // assuming dx = 1

      } else {
         //upwinding
         stencil::Direction upwindingXDirection = (targetVelocity[0] > real_t(0) ) ? stencil::W : stencil::E;
         stencil::Direction sourceXDirection = (availableStencilIndices[Stencil::idx[upwindingXDirection]]) ? upwindingXDirection
                                                                                                            : ((availableStencilIndices[Stencil::idx[stencil::inverseDir[upwindingXDirection]]]) ? stencil::inverseDir[upwindingXDirection]
                                                                                                                                                                                                 : stencil::C );
         if(sourceXDirection == stencil::E)
         {
            auto neighborVelocity = pdfField->getVelocity(Cell(x,y,z)+sourceXDirection);
            velocityGradient[0] = neighborVelocity[0] - targetVelocity[0]; // assuming dx = 1
            velocityGradient[1] = neighborVelocity[1] - targetVelocity[1]; // assuming dx = 1
            velocityGradient[2] = neighborVelocity[2] - targetVelocity[2]; // assuming dx = 1
         }
         if(sourceXDirection == stencil::W)
         {
            auto neighborVelocity = pdfField->getVelocity(Cell(x,y,z)+sourceXDirection);
            velocityGradient[0] = targetVelocity[0] - neighborVelocity[0]; // assuming dx = 1
            velocityGradient[1] = targetVelocity[1] - neighborVelocity[1]; // assuming dx = 1
            velocityGradient[2] = targetVelocity[2] - neighborVelocity[2]; // assuming dx = 1
         }
         // else: 0
      }

      // check if both neighbors are available for central differences
      if(useCentralDifferences_ && availableStencilIndices[Stencil::idx[stencil::N]] && availableStencilIndices[Stencil::idx[stencil::S]]) {
         auto neighborVelocity1 = pdfField->getVelocity(Cell(x, y, z) + stencil::N);
         auto neighborVelocity2 = pdfField->getVelocity(Cell(x, y, z) + stencil::S);
         velocityGradient[3] = real_t(0.5) * (neighborVelocity1[0] - neighborVelocity2[0]); // assuming dx = 1
         velocityGradient[4] = real_t(0.5) * (neighborVelocity1[1] - neighborVelocity2[1]); // assuming dx = 1
         velocityGradient[5] = real_t(0.5) * (neighborVelocity1[2] - neighborVelocity2[2]); // assuming dx = 1
      } else {
         //upwinding
         stencil::Direction upwindingYDirection = (targetVelocity[1] > real_t(0) ) ? stencil::S : stencil::N;
         stencil::Direction sourceYDirection = (availableStencilIndices[Stencil::idx[upwindingYDirection]]) ? upwindingYDirection
                                                                                                            : ((availableStencilIndices[Stencil::idx[stencil::inverseDir[upwindingYDirection]]]) ? stencil::inverseDir[upwindingYDirection]
                                                                                                                                                                                                 : stencil::C );
         if(sourceYDirection == stencil::N)
         {
            auto neighborVelocity = pdfField->getVelocity(Cell(x,y,z)+sourceYDirection);
            velocityGradient[3] = neighborVelocity[0] - targetVelocity[0]; // assuming dx = 1
            velocityGradient[4] = neighborVelocity[1] - targetVelocity[1]; // assuming dx = 1
            velocityGradient[5] = neighborVelocity[2] - targetVelocity[2]; // assuming dx = 1
         }
         if(sourceYDirection == stencil::S)
         {
            auto neighborVelocity = pdfField->getVelocity(Cell(x,y,z)+sourceYDirection);
            velocityGradient[3] = targetVelocity[0] - neighborVelocity[0]; // assuming dx = 1
            velocityGradient[4] = targetVelocity[1] - neighborVelocity[1]; // assuming dx = 1
            velocityGradient[5] = targetVelocity[2] - neighborVelocity[2]; // assuming dx = 1
         }
         // else: 0
      }

      if(Stencil::D == 3)
      {
         // only in 3D
         // check if both neighbors are available for central differences
         if(useCentralDifferences_ && availableStencilIndices[Stencil::idx[stencil::T]] && availableStencilIndices[Stencil::idx[stencil::B]]) {
            auto neighborVelocity1 = pdfField->getVelocity(Cell(x, y, z) + stencil::T);
            auto neighborVelocity2 = pdfField->getVelocity(Cell(x, y, z) + stencil::B);
            velocityGradient[6] = real_t(0.5) * (neighborVelocity1[0] - neighborVelocity2[0]); // assuming dx = 1
            velocityGradient[7] = real_t(0.5) * (neighborVelocity1[1] - neighborVelocity2[1]); // assuming dx = 1
            velocityGradient[8] = real_t(0.5) * (neighborVelocity1[2] - neighborVelocity2[2]); // assuming dx = 1
         } else {
            //upwinding
            stencil::Direction upwindingZDirection = (targetVelocity[2] > real_t(0)) ? stencil::B : stencil::T;
            stencil::Direction sourceZDirection = (availableStencilIndices[Stencil::idx[upwindingZDirection]]) ? upwindingZDirection
                                                                                                               : ((availableStencilIndices[Stencil::idx[stencil::inverseDir[upwindingZDirection]]]) ? stencil::inverseDir[upwindingZDirection]
                                                                                                                                                                                                    : stencil::C);
            if (sourceZDirection == stencil::T) {
               auto neighborVelocity = pdfField->getVelocity(Cell(x, y, z) + sourceZDirection);
               velocityGradient[6] = neighborVelocity[0] - targetVelocity[0]; // assuming dx = 1
               velocityGradient[7] = neighborVelocity[1] - targetVelocity[1]; // assuming dx = 1
               velocityGradient[8] = neighborVelocity[2] - targetVelocity[2]; // assuming dx = 1
            }
            if (sourceZDirection == stencil::B) {
               auto neighborVelocity = pdfField->getVelocity(Cell(x, y, z) + sourceZDirection);
               velocityGradient[6] = targetVelocity[0] - neighborVelocity[0]; // assuming dx = 1
               velocityGradient[7] = targetVelocity[1] - neighborVelocity[1]; // assuming dx = 1
               velocityGradient[8] = targetVelocity[2] - neighborVelocity[2]; // assuming dx = 1
            }
            // else: 0
         }
      }


      Matrix3<real_t> pressureTensorNeq(real_t(0)); // without prefactor of rho, added later
      const real_t preFac = - real_t(1)  / (real_t(3) * omegaShear_); // 2 * beta (in Chikatamarla et al) = omega related to kinematic viscosity
      for(auto j = uint_t(0); j <= uint_t(2); ++j)
      {
         for(auto i = uint_t(0); i <= uint_t(2); ++i)
         {
            pressureTensorNeq(i,j) += preFac * (velocityGradient(i,j)+velocityGradient(j,i));
         }
      }

      // 5. set PDFs to approximated equilibrium by Grad
      // this is just the regular feq with additional Pneq part
      if(LatticeModel::compressible)
      {
         for( auto q = Stencil::begin(); q != Stencil::end(); ++q )
         {
            const real_t velci = lbm::internal::multiplyVelocityDirection( *q, targetVelocity );

            auto contributionFromPneq = real_t(0);
            for(auto j = uint_t(0); j <= uint_t(2); ++j) {
               for (auto i = uint_t(0); i <= uint_t(2); ++i) {
                  //Pneq_a,b * c_q,a * c_q,b
                  contributionFromPneq += pressureTensorNeq(i,j) * real_c(stencil::c[i][*q]) * real_c(stencil::c[j][*q]);
               }
            }
            //- Pneq_a,b * cs**2 * delta_a,b
            contributionFromPneq -= (pressureTensorNeq(0,0) + pressureTensorNeq(1,1) + pressureTensorNeq(2,2)) / real_t(3);

            // all terms are multiplied with the density
            real_t fGrad = LatticeModel::w[ q.toIdx() ] * targetDensity * (real_t(1) + real_t(3) * velci - real_t(1.5) * targetVelocity.sqrLength() + real_t(4.5) * velci * velci + real_t(4.5) * contributionFromPneq); // standard comp. feq + comp. Pneq
            pdfField->get(x,y,z,*q) = fGrad;
         }
      }else{
         // assume rho = 1 almost everywhere, targetDensity is just the deviation from 1
         for( auto q = Stencil::begin(); q != Stencil::end(); ++q )
         {
            const real_t velci = lbm::internal::multiplyVelocityDirection( *q, targetVelocity );

            auto contributionFromPneq = real_t(0);
            for(auto j = uint_t(0); j <= uint_t(2); ++j) {
               for (auto i = uint_t(0); i <= uint_t(2); ++i) {
                  //Pneq_a,b * c_q,a * c_q,b
                  contributionFromPneq += pressureTensorNeq(i,j) * real_c(stencil::c[i][*q]) * real_c(stencil::c[j][*q]);
               }
            }
            //- Pneq_a,b * cs**2 * delta_a,b
            contributionFromPneq -= (pressureTensorNeq(0,0) + pressureTensorNeq(1,1) + pressureTensorNeq(2,2)) / real_t(3);

            // density deviation just appears as the leading order term
            real_t fGrad = LatticeModel::w[ q.toIdx() ] * ( targetDensity + real_t(3) * velci - real_t(1.5) * targetVelocity.sqrLength() + real_t(4.5) * velci * velci + real_t(4.5) * contributionFromPneq); // standard incomp. feq + incomp Pneq
            pdfField->get(x,y,z,*q) = fGrad;
         }
      }
   }

   void setOmegaShear(real_t omegaShear)
   {
      omegaShear_ = omegaShear;
   }

   real_t getOmegaShear()
   {
      return omegaShear_;
   }

private:

   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID boundaryHandlingID_;
   real_t omegaShear_;
   const bool recomputeTargetDensity_;
   const bool useCentralDifferences_;
   const bool useDataFromGhostLayers_;
};



// make functionality

template< typename BoundaryHandling_T>
shared_ptr<EquilibriumReconstructor<BoundaryHandling_T> >
makeEquilibriumReconstructor(const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & boundaryHandlingID, bool useDataFromGhostLayers = false)
{
   using Rec_T = EquilibriumReconstructor<BoundaryHandling_T>;
   return make_shared<Rec_T>(blockStorage, boundaryHandlingID, useDataFromGhostLayers);
}

template< typename BoundaryHandling_T, typename ExtrapolationDirectionFinder_T >
shared_ptr<EquilibriumAndNonEquilibriumReconstructor<BoundaryHandling_T,ExtrapolationDirectionFinder_T> >
makeEquilibriumAndNonEquilibriumReconstructor(const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & boundaryHandlingID, const shared_ptr<ExtrapolationDirectionFinder_T> & extrapolationDirectionFinder,
                                              uint_t maximumNumberOfExtrapolationCells = uint_t(3), bool useDataFromGhostLayers = false)
{
   using Rec_T = EquilibriumAndNonEquilibriumReconstructor<BoundaryHandling_T,ExtrapolationDirectionFinder_T>;
   return make_shared<Rec_T>(blockStorage, boundaryHandlingID, extrapolationDirectionFinder, maximumNumberOfExtrapolationCells, useDataFromGhostLayers);
}

template< typename BoundaryHandling_T, typename ExtrapolationDirectionFinder_T, bool EnforceNoSlipConstraintAfterExtrapolation = false >
shared_ptr<ExtrapolationReconstructor<BoundaryHandling_T,ExtrapolationDirectionFinder_T, EnforceNoSlipConstraintAfterExtrapolation> >
makeExtrapolationReconstructor(const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & boundaryHandlingID, const shared_ptr<ExtrapolationDirectionFinder_T> & extrapolationDirectionFinder,
                               uint_t maximumNumberOfExtrapolationCells = uint_t(3), bool useDataFromGhostLayers = false)
{
   using Rec_T = ExtrapolationReconstructor<BoundaryHandling_T,ExtrapolationDirectionFinder_T, EnforceNoSlipConstraintAfterExtrapolation>;
   return make_shared<Rec_T>(blockStorage, boundaryHandlingID, extrapolationDirectionFinder, maximumNumberOfExtrapolationCells, useDataFromGhostLayers);
}

template< typename BoundaryHandling_T>
shared_ptr<GradsMomentApproximationReconstructor<BoundaryHandling_T> >
makeGradsMomentApproximationReconstructor(const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & boundaryHandlingID, real_t omegaShear,
                                          bool recomputeTargetDensity = false, bool useCentralDifferences = true, bool useDataFromGhostLayers = false)
{
   using Rec_T = GradsMomentApproximationReconstructor<BoundaryHandling_T>;
   return make_shared<Rec_T>(blockStorage, boundaryHandlingID, omegaShear, recomputeTargetDensity, useCentralDifferences, useDataFromGhostLayers);
}


} // namespace lbm_mesapd_coupling
} // namespace walberla
