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
//! \file ParserUBB.h
//! \ingroup lbm
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
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
#include "core/math/ParserOMP.h"

#include "field/FlagUID.h"

#include "stencil/Directions.h"

#include <vector>



namespace walberla {
namespace lbm {



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce = false >
class ParserUBB : public Boundary<flag_t>
{
   typedef lbm::PdfField< LatticeModel_T >   PDFField;
   typedef typename LatticeModel_T::Stencil  Stencil;

public:

   static const bool threadsafe = true;

   class Parser : public BoundaryConfiguration
   {
   public:
      inline Parser( const Config::BlockHandle & config );
      Vector3< real_t > operator()( const Vector3< real_t > & x, const real_t t ) const;
      Vector3< real_t > operator()( const Vector3< real_t > & x ) const;
      bool isTimeDependent() const { return timeDependent_; }

   private:
      std::array< math::FunctionParserOMP, 3 > parsers_;
      bool timeDependent_;
   }; // class Parser

   typedef GhostLayerField< shared_ptr<Parser>, 1 > ParserField;
   typedef GhostLayerField< Vector3<real_t>, 1 >    VelocityField;

   static shared_ptr<Parser> createConfiguration( const Config::BlockHandle & config )
      { return make_shared<Parser>( config ); }



   ParserUBB( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField * const pdfField,
              FlagField<flag_t> * const flagField, const shared_ptr< TimeTracker > & timeTracker,
              const uint_t level, const AABB & aabb );
   ParserUBB( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField * const pdfField,
              FlagField<flag_t> * const flagField,  const uint_t level, const AABB & aabb );

   shared_ptr< TimeTracker > getTimeTracker() { return timeTracker_; }

   void pushFlags( std::vector< FlagUID > & uids ) const { uids.push_back( uid_ ); }

   void beforeBoundaryTreatment()
   {
      if( timeTracker_ )
         time_ = timeTracker_->getTime( level_ );
   }
   void  afterBoundaryTreatment() const {}

   template< typename Buffer_T >
   void packCell( Buffer_T &, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

   template< typename Buffer_T >
   void registerCell( Buffer_T &, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) {}

   inline void registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const BoundaryConfiguration & parser );
   inline void registerCells( const flag_t, const CellInterval & cells, const BoundaryConfiguration & parser );
   template< typename CellIterator >
   inline void registerCells( const flag_t, const CellIterator & begin, const CellIterator & end, const BoundaryConfiguration & parser );

   void unregisterCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

#ifndef NDEBUG
   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask );
#else
   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ );
#endif

private:

   const FlagUID uid_;

   PDFField * const pdfField_;

   Vector3< real_t > origin_;
   Vector3< real_t > dx_;

   // required to keep track of the simulation time
   shared_ptr< TimeTracker > timeTracker_; // -> addPostBoundaryHandlingVoidFunction (when used with refinement time step)
   real_t time_;
   uint_t level_;

   shared_ptr<ParserField> parserField_;
   shared_ptr<VelocityField> velocityField_;

}; // class ParserUBB

template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce>
inline ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce>::Parser::Parser( const Config::BlockHandle & config )
: parsers_(), timeDependent_( false )
{
   if( !config )
      return;

   if( config.isDefined( "x" ) )
   {
      parsers_[0].parse( config.getParameter<std::string>( "x" ) );
      if( parsers_[0].symbolExists( "t" ) )
         timeDependent_ = true;
   }
   if( config.isDefined( "y" ) )
   {
      parsers_[1].parse( config.getParameter<std::string>( "y" ) );
      if( parsers_[1].symbolExists( "t" ) )
         timeDependent_ = true;
   }
   if( config.isDefined( "z" ) )
   {
      parsers_[2].parse( config.getParameter<std::string>( "z" ) );
      if( parsers_[2].symbolExists( "t" ) )
         timeDependent_ = true;
   }
}

template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce>
Vector3< real_t > ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce>::Parser::operator()( const Vector3< real_t > & x, const real_t t ) const
{
   std::map< std::string, double > symbols;
   symbols["x"] = x[0];
   symbols["y"] = x[1];
   symbols["z"] = x[2];
   symbols["t"] = t;

   Vector3< real_t > v;
   v[0] = parsers_[0].evaluate( symbols );
   v[1] = parsers_[1].evaluate( symbols );
   v[2] = parsers_[2].evaluate( symbols );
   return v;
}

template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce>
Vector3< real_t > ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce>::Parser::operator()( const Vector3< real_t > & x ) const
{
   WALBERLA_ASSERT( !timeDependent_ );

   std::map< std::string, double > symbols;
   symbols["x"] = x[0];
   symbols["y"] = x[1];
   symbols["z"] = x[2];

   Vector3< real_t > v;
   v[0] = parsers_[0].evaluate( symbols );
   v[1] = parsers_[1].evaluate( symbols );
   v[2] = parsers_[2].evaluate( symbols );
   return v;
}



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce>
inline ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce>::ParserUBB::ParserUBB( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField * const pdfField,
                                                                                           FlagField<flag_t> * const flagField, const shared_ptr< TimeTracker > & timeTracker,
                                                                                           const uint_t level, const AABB & aabb )
   : Boundary<flag_t>( boundaryUID ), uid_( uid ), pdfField_( pdfField ), timeTracker_( timeTracker ), time_( real_t(0) ), level_( level )
{
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField_ );
   dx_[0] = aabb.xSize() / real_c( pdfField_->xSize() );
   dx_[1] = aabb.ySize() / real_c( pdfField_->ySize() );
   dx_[2] = aabb.zSize() / real_c( pdfField_->zSize() );
   origin_[0] = aabb.xMin() + real_c(0.5) * dx_[0];
   origin_[1] = aabb.yMin() + real_c(0.5) * dx_[1];
   origin_[2] = aabb.zMin() + real_c(0.5) * dx_[2];

   if(flagField != NULL)
   {
      parserField_   = make_shared<ParserField>  ( pdfField->xSize(), pdfField->ySize(), pdfField->zSize(), flagField->nrOfGhostLayers(), field::zyxf );
      velocityField_ = make_shared<VelocityField>( pdfField->xSize(), pdfField->ySize(), pdfField->zSize(), flagField->nrOfGhostLayers(), field::zyxf );
   }
   else
   {
      parserField_   = make_shared<ParserField>  ( pdfField->xSize(), pdfField->ySize(), pdfField->zSize(), pdfField->nrOfGhostLayers(),  field::zyxf );
      velocityField_ = make_shared<VelocityField>( pdfField->xSize(), pdfField->ySize(), pdfField->zSize(), pdfField->nrOfGhostLayers(),  field::zyxf );
   }
}



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce>
inline ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce>::ParserUBB::ParserUBB( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField * const pdfField,
                                                                                           FlagField<flag_t> * const flagField, const uint_t level, const AABB & aabb )
   : ParserUBB( boundaryUID, uid, pdfField, flagField, nullptr, level, aabb )
{}



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce >
inline void ParserUBB< LatticeModel_T, flag_t, AdaptVelocityToExternalForce >::registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                                             const BoundaryConfiguration & parser )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Parser * >( &parser ), &parser );
   WALBERLA_ASSERT_NOT_NULLPTR( parserField_ );

   auto & p = dynamic_cast< const Parser & >( parser );

   if( p.isTimeDependent() )
   {
      parserField_->get( x, y, z ) = make_shared<Parser>(p);
   }
   else
   {
      const Vector3< real_t > pos( origin_[0] + real_c(x) * dx_[0],
                                   origin_[1] + real_c(y) * dx_[1],
                                   origin_[2] + real_c(z) * dx_[2] );
      velocityField_->get( x, y, z ) = p(pos);
      parserField_->get( x, y, z ) = nullptr;
   }
}



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce >
inline void ParserUBB< LatticeModel_T, flag_t, AdaptVelocityToExternalForce >::registerCells( const flag_t, const CellInterval & cells, const BoundaryConfiguration & parser )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Parser * >( &parser ), &parser );
   WALBERLA_ASSERT_NOT_NULLPTR( parserField_ );

   auto & p = dynamic_cast< const Parser & >( parser );

   if( p.isTimeDependent() )
   {
      auto shared_p = make_shared<Parser>(p);
      for( auto cell = parserField_->beginSliceXYZ( cells ); cell != parserField_->end(); ++cell )
         *cell = shared_p;
   }
   else
   {
      for( auto cell = parserField_->beginSliceXYZ( cells ); cell != parserField_->end(); ++cell )
      {
         const Vector3< real_t > pos( origin_[0] + real_c(cell.x()) * dx_[0],
                                      origin_[1] + real_c(cell.y()) * dx_[1],
                                      origin_[2] + real_c(cell.z()) * dx_[2] );
         velocityField_->get( cell.x(), cell.y(), cell.z() ) = p(pos);
         *cell = nullptr;
      }
   }
}



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce >
template< typename CellIterator >
inline void ParserUBB< LatticeModel_T, flag_t, AdaptVelocityToExternalForce >::registerCells( const flag_t, const CellIterator & begin, const CellIterator & end,
                                                                                              const BoundaryConfiguration & parser )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Parser * >( &parser ), &parser );
   WALBERLA_ASSERT_NOT_NULLPTR( parserField_ );

   auto & p = dynamic_cast< const Parser & >( parser );

   if( p.isTimeDependent() )
   {
      auto shared_p = make_shared<Parser>(p);
      for( auto cell = begin; cell != end; ++cell )
      {
         parserField_->get( cell->x(), cell->y(), cell->z() ) = shared_p;
      }
   }
   else
   {
      for( auto cell = begin; cell != end; ++cell )
      {
         const Vector3< real_t > pos( origin_[0] + real_c(cell->x()) * dx_[0],
                                      origin_[1] + real_c(cell->y()) * dx_[1],
                                      origin_[2] + real_c(cell->z()) * dx_[2] );
         velocityField_->get( cell->x(), cell->y(), cell->z() ) = p(pos);
         parserField_->get( cell->x(), cell->y(), cell->z() ) = nullptr;
      }
   }
}



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce>
#ifndef NDEBUG
   inline void ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce>::ParserUBB::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                                           const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
   inline void ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce>::ParserUBB::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                                           const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ )
#endif
   {
      WALBERLA_ASSERT_EQUAL( nx, x + cell_idx_c( stencil::cx[ dir ] ) );
      WALBERLA_ASSERT_EQUAL( ny, y + cell_idx_c( stencil::cy[ dir ] ) );
      WALBERLA_ASSERT_EQUAL( nz, z + cell_idx_c( stencil::cz[ dir ] ) );
      WALBERLA_ASSERT_UNEQUAL( mask & this->mask_, numeric_cast<flag_t>(0) );
      WALBERLA_ASSERT_EQUAL( mask & this->mask_, this->mask_ ); // only true if "this->mask_" only contains one single flag, which is the case for the
                                                                // current implementation of this boundary condition (ParserUBB)

      const Vector3< real_t > pos( origin_[0] + real_c(nx) * dx_[0],
                                   origin_[1] + real_c(ny) * dx_[1],
                                   origin_[2] + real_c(nz) * dx_[2] );

      Vector3<real_t> velocity;
      if( parserField_->get(nx,ny,nz) )
      {
         WALBERLA_ASSERT_NOT_NULLPTR( getTimeTracker() );
         velocity = (*parserField_->get(nx,ny,nz))( pos, time_ );
      }
      else
      {
         velocity = velocityField_->get(nx,ny,nz);
      }

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
   }


} // namespace lbm
} // namespace walberla
