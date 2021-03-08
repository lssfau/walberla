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
//! \file DynamicUBB.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/DensityAndVelocity.h"
#include "lbm/field/PdfField.h"
#include "lbm/refinement/TimeTracker.h"

#include "boundary/Boundary.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/math/Vector3.h"

#include "field/FlagUID.h"

#include "stencil/Directions.h"

#include <vector>



namespace walberla {
namespace lbm {



// VelocityFunctor_T: functor that requires to implement two member functions:
//   1. A member function "void operator()( const real_t t )" that is called once before the boundary treatement with the current time
//   2. A member function "Vector3< real_t > operator()( const Vector3< real_t > & x, const real_t t )" that is called for every
//      boundary link treated by "treatDirection". The arguments are the position 'x' of the boudnary cell in the simulation space and the current time 't'.
//      The functon is supposed to return the velocity used by the boundary treatment.
template< typename LatticeModel_T, typename flag_t, typename VelocityFunctor_T, bool AdaptVelocityToExternalForce = false, bool StoreForce = false >
class DynamicUBB : public Boundary<flag_t>
{
   using PDFField = lbm::PdfField<LatticeModel_T>;
   using Stencil = typename LatticeModel_T::Stencil;

   using ForceField = GhostLayerField<Vector3<real_t>, 1>;

public:

   static const bool threadsafe = true;

   static shared_ptr<BoundaryConfiguration> createConfiguration( const Config::BlockHandle & /*config*/ )
      { return make_shared<BoundaryConfiguration>(); }



   DynamicUBB( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField * const pdfField,
               const shared_ptr< TimeTracker > & timeTracker, const uint_t level, const VelocityFunctor_T & velocity, const AABB & aabb ) :
      Boundary<flag_t>( boundaryUID ), uid_( uid ), pdfField_( pdfField ),
      timeTracker_( timeTracker ), time_( real_t(0) ), level_( level ), velocity_( velocity )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ );
      dx_[0] = aabb.xSize() / real_c( pdfField_->xSize() );
      dx_[1] = aabb.ySize() / real_c( pdfField_->ySize() );
      dx_[2] = aabb.zSize() / real_c( pdfField_->zSize() );
      origin_[0] = aabb.xMin() + real_c(0.5) * dx_[0];
      origin_[1] = aabb.yMin() + real_c(0.5) * dx_[1];
      origin_[2] = aabb.zMin() + real_c(0.5) * dx_[2];

      if( !timeTracker_ )
      {
         velocity_( time_ );
      }

      if (StoreForce)
         force_ = make_shared<ForceField>( pdfField_->xSize(), pdfField_->ySize(), pdfField_->zSize(), pdfField_->nrOfGhostLayers(), field::zyxf );
   }
   DynamicUBB( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField * const pdfField,
               const uint_t level, const VelocityFunctor_T & velocity, const AABB & aabb ) :
      DynamicUBB( boundaryUID, uid, pdfField, nullptr, level, velocity, aabb )
   {}

   shared_ptr< TimeTracker > getTimeTracker() { return timeTracker_; }

   void pushFlags( std::vector< FlagUID > & uids ) const { uids.push_back( uid_ ); }

   void beforeBoundaryTreatment()
   {
     if( timeTracker_ )
     {
        time_ = timeTracker_->getTime( level_ );
        velocity_( time_ );
     }

     if (StoreForce)
        force_->setWithGhostLayer( Vector3<real_t>() );
   }
   void  afterBoundaryTreatment() const {}

   template< typename Buffer_T >
   void packCell( Buffer_T &, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

   template< typename Buffer_T >
   void registerCell( Buffer_T &, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) {}

   void registerCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t, const BoundaryConfiguration & ) {}
   void registerCells( const flag_t, const CellInterval &, const BoundaryConfiguration & ) const {}
   template< typename CellIterator >
   void registerCells( const flag_t, const CellIterator &, const CellIterator &, const BoundaryConfiguration & ) const {}

   void unregisterCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

#ifndef NDEBUG
   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ )
#endif
   {
      WALBERLA_ASSERT_EQUAL( nx, x + cell_idx_c( stencil::cx[ dir ] ) );
      WALBERLA_ASSERT_EQUAL( ny, y + cell_idx_c( stencil::cy[ dir ] ) );
      WALBERLA_ASSERT_EQUAL( nz, z + cell_idx_c( stencil::cz[ dir ] ) );
      WALBERLA_ASSERT_UNEQUAL( mask & this->mask_, numeric_cast<flag_t>(0) );
      WALBERLA_ASSERT_EQUAL( mask & this->mask_, this->mask_ ); // only true if "this->mask_" only contains one single flag, which is the case for the
                                                                // current implementation of this boundary condition (DynamicUBB)

      const real_t pdf_old = pdfField_->get( x, y, z, Stencil::idx[dir] );

      const Vector3< real_t > pos( origin_[0] + real_c(nx) * dx_[0],
                                   origin_[1] + real_c(ny) * dx_[1],
                                   origin_[2] + real_c(nz) * dx_[2] );

      const auto velocity = velocity_( pos, time_ );

      if( LatticeModel_T::compressible )
      {
         const auto density  = pdfField_->getDensity(x,y,z);
         const auto vel = AdaptVelocityToExternalForce ? lbm::internal::AdaptVelocityToForce<LatticeModel_T>::get( x, y, z, pdfField_->latticeModel(), velocity, density ) :
                                                         velocity;

         pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) ) = pdfField_->get( x, y, z, Stencil::idx[dir] ) -
                                                                 ( real_c(6) * density * real_c(LatticeModel_T::w[ Stencil::idx[dir] ]) *
                                                                    ( real_c(stencil::cx[ dir ]) * vel[0] +
                                                                      real_c(stencil::cy[ dir ]) * vel[1] +
                                                                      real_c(stencil::cz[ dir ]) * vel[2] ) );
      }
      else
      {
         const auto vel = AdaptVelocityToExternalForce ? lbm::internal::AdaptVelocityToForce<LatticeModel_T>::get( x, y, z, pdfField_->latticeModel(), velocity, real_t(1) ) :
                                                         velocity;

         pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) ) = pdfField_->get( x, y, z, Stencil::idx[dir] ) -
                                                                 ( real_c(6) * real_c(LatticeModel_T::w[ Stencil::idx[dir] ]) *
                                                                    ( real_c(stencil::cx[ dir ]) * vel[0] +
                                                                      real_c(stencil::cy[ dir ]) * vel[1] +
                                                                      real_c(stencil::cz[ dir ]) * vel[2] ) );
      }

      if (StoreForce && pdfField_->isInInnerPart( Cell(x,y,z) ))
      {
         const real_t forceMEM = pdf_old + pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) );
         Vector3<real_t> force( real_c( stencil::cx[dir] ) * forceMEM,
                                real_c( stencil::cy[dir] ) * forceMEM,
                                real_c( stencil::cz[dir] ) * forceMEM );
         force_->get( nx, ny, nz ) += force;
      }
   }

   const typename ForceField::value_type & getForce( const cell_idx_t x, cell_idx_t y, cell_idx_t z ) const
   {
      static_assert(StoreForce, "this member function is only available if the fourth template argument on the class is true");
      return force_->get(x,y,z);
   }

private:

   const FlagUID uid_;

   PDFField * const pdfField_;

   Vector3< real_t > origin_;
   Vector3< real_t > dx_;

   // required to keep track of the simulation time
   shared_ptr< TimeTracker > timeTracker_; // -> addPostBoundaryHandlingVoidFunction (when used with refinement time step)
   real_t time_;
   uint_t level_;

   VelocityFunctor_T velocity_;
   shared_ptr<ForceField> force_;

}; // class DynamicUBB



} // namespace lbm
} // namespace walberla
