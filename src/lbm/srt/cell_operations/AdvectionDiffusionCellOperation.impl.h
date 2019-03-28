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
//! \file AdvectionDiffusionCellOperation.impl.h
//! \ingroup lbm
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#include "lbm/field/MacroscopicValueCalculation.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"
#include "lbm/lattice_model/LatticeModelBase.h"

#include <type_traits>


namespace walberla {
namespace lbm {


   ////////////////////////////////////////////////////////////////////////
// Available SRT implementations:                                     //
//                                                                    //
// Generic (D*Q*) version:                                            //
//                                      incompressible | compressible //
//                           no forces:       x               x       //
//                                                                    //
// Optimized D3Q19 implementation:                                    //
//                                      incompressible | compressible //
//                           no forces:       x               x       //
////////////////////////////////////////////////////////////////////////


///////////////////////////////
// Specialization for:       //
// - no additional forces    //
///////////////////////////////

template< typename LM_AdvDiff, typename LM_Hydro >
class AdvectionDiffusionCellOperation< LM_AdvDiff, LM_Hydro, typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                                                      LM_AdvDiff::CollisionModel::constant &&
                                                                                      LM_AdvDiff::compressible &&
                                                                                      std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value
                                                                                      >::type >
{
public:

   static_assert( (std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( LM_AdvDiff::compressible,                                                                      "Only works with compressible models!" );
   static_assert( (std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::None_tag >::value),        "Only works without additional forces!" );

   typedef PdfField< LM_AdvDiff >        AdvDiffPdfField_T;
   typedef PdfField< LM_Hydro   >        HydroPdfField_T;
   typedef typename LM_AdvDiff::Stencil  Stencil;

   AdvectionDiffusionCellOperation() : omega_( real_t(0) ), advDiffLatticeModel_( NULL ), hydroLatticeModel_(NULL) {}

   void configure( const LM_AdvDiff & advDiffLatticeModel, const LM_Hydro & hydroLatticeModel )
   {
      omega_ = advDiffLatticeModel.collisionModel().omega();
      advDiffLatticeModel_ = &advDiffLatticeModel;
      hydroLatticeModel_   = &hydroLatticeModel;
   }

   void operator()( AdvDiffPdfField_T * src, AdvDiffPdfField_T * dst, HydroPdfField_T * hydro, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const
   {
      // stream pull
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         dst->get( x,y,z,d.toIdx() ) = src->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );

      Vector3<real_t> velocity = hydro->getVelocity( x, y, z );
      real_t scalar = dst->getDensity( x, y, z );

      // collide
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         dst->get( x, y, z, d.toIdx() ) = ( real_t(1.0) - omega_ ) * dst->get( x, y, z, d.toIdx() ) +
            omega_ * EquilibriumDistribution< LM_AdvDiff >::get( *d, velocity, scalar );
      }
   }

   template< typename AdvDiffFieldPtrOrIterator, typename HydroFieldPtrOrIterator >
   void operator()( AdvDiffFieldPtrOrIterator & src, AdvDiffFieldPtrOrIterator & dst, HydroFieldPtrOrIterator & hydro ) const
   {
      // stream pull
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         dst[ d.toIdx() ] = src.neighbor( d.inverseDir(), d.toIdx() );

      Vector3<real_t> velocity = getVelocity( *hydroLatticeModel_, hydro );
      real_t scalar = getDensity( *advDiffLatticeModel_, dst );

      // collide
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         dst[ d.toIdx() ] = ( real_t(1.0) - omega_ ) * dst[ d.toIdx() ] +
            omega_   * EquilibriumDistribution< LM_AdvDiff >::get( *d, velocity, scalar );
      }
   }

private:

   real_t omega_;
   const LM_AdvDiff * advDiffLatticeModel_;
   const LM_Hydro   * hydroLatticeModel_;
};



///////////////////////////////
// Specialization for:       //
// - correction force        //
///////////////////////////////

template< typename LM_AdvDiff, typename LM_Hydro >
class AdvectionDiffusionCellOperation< LM_AdvDiff, LM_Hydro, typename std::enable_if< std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value &&
                                                                                      LM_AdvDiff::CollisionModel::constant &&
                                                                                      LM_AdvDiff::compressible &&
                                                                                      std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::Correction_tag >::value
                                                                                      >::type >
{
public:

   static_assert( (std::is_same< typename LM_AdvDiff::CollisionModel::tag, collision_model::SRT_tag >::value), "Only works with SRT!" );
   static_assert( LM_AdvDiff::compressible,                                                                      "Only works with compressible models!" );
   static_assert( (std::is_same< typename LM_AdvDiff::ForceModel::tag, force_model::Correction_tag >::value),  "Only works with correction force!" );

   typedef PdfField< LM_AdvDiff >        AdvDiffPdfField_T;
   typedef PdfField< LM_Hydro   >        HydroPdfField_T;
   typedef typename LM_AdvDiff::Stencil  Stencil;

   AdvectionDiffusionCellOperation() : omega_( real_t(0) ), advDiffLatticeModel_( NULL ), hydroLatticeModel_(NULL) {}

   void configure( const LM_AdvDiff & advDiffLatticeModel, const LM_Hydro & hydroLatticeModel )
   {
      omega_ = advDiffLatticeModel.collisionModel().omega();
      advDiffLatticeModel_ = &advDiffLatticeModel;
      hydroLatticeModel_   = &hydroLatticeModel;
   }

   void operator()( AdvDiffPdfField_T * src, AdvDiffPdfField_T * dst, HydroPdfField_T * hydro, Vector3<real_t> dtsv, cell_idx_t x, cell_idx_t y, cell_idx_t z ) const
   {
      // stream pull
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         dst->get( x,y,z,d.toIdx() ) = src->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );

      Vector3<real_t> velocity = hydro->getVelocity( x, y, z );
      real_t scalar = dst->getDensity( x, y, z );
      dtsv *= real_t(3) - real_c(1.5)*omega_;

      // collide
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         dst->get( x, y, z, d.toIdx() ) = ( real_t(1.0) - omega_ ) * dst->get( x, y, z, d.toIdx() ) +
                                                          omega_   * EquilibriumDistribution< LM_AdvDiff >::get( *d, velocity, scalar ) +
                                                                     LM_AdvDiff::w[ d.toIdx() ] * ( d.cx()*dtsv[0] + d.cy()*dtsv[1] + d.cz()*dtsv[2] );
      }
   }

   template< typename AdvDiffFieldPtrOrIterator, typename HydroFieldPtrOrIterator >
   void operator()( AdvDiffFieldPtrOrIterator & src, AdvDiffFieldPtrOrIterator & dst, HydroFieldPtrOrIterator & hydro, Vector3<real_t> dtsv ) const
   {
      // stream pull
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         dst[ d.toIdx() ] = src.neighbor( d.inverseDir(), d.toIdx() );

      Vector3<real_t> velocity = getVelocity( *hydroLatticeModel_, hydro );
      real_t scalar = getDensity( *advDiffLatticeModel_, dst );
      dtsv *= real_t(3) - real_c(1.5)*omega_;

      // collide
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         dst[ d.toIdx() ] = ( real_t(1.0) - omega_ ) * dst[ d.toIdx() ] +
                                            omega_   * EquilibriumDistribution< LM_AdvDiff >::get( *d, velocity, scalar ) +
                                                       LM_AdvDiff::w[ d.toIdx() ] * ( d.cx()*dtsv[0] + d.cy()*dtsv[1] + d.cz()*dtsv[2] );
      }
   }

private:

   real_t omega_;
   const LM_AdvDiff * advDiffLatticeModel_;
   const LM_Hydro   * hydroLatticeModel_;
};




} // namespace lbm
} // namespace walberla
