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



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce = false, bool StoreForce = false >
class ParserUBB : public Boundary<flag_t>
{
   typedef lbm::PdfField< LatticeModel_T >   PDFField;
   typedef typename LatticeModel_T::Stencil  Stencil;

   typedef GhostLayerField< Vector3<real_t>, 1 > ForceField;

public:

   static const bool threadsafe = true;

   class Parser : public BoundaryConfiguration
   {
   public:
      inline Parser( const Config::BlockHandle & config );
      inline Parser( const std::array< std::string, 3 > & equations );
      Vector3< real_t > operator()( const Vector3< real_t > & x, const real_t t ) const;
      Vector3< real_t > operator()( const Vector3< real_t > & x ) const;
      bool isTimeDependent() const { return timeDependent_; }
      const std::array< std::string, 3 > & equations() const { return equations_; }

   private:
      std::array< math::FunctionParserOMP, 3 > parsers_;
      std::array< std::string, 3 > equations_;
      bool timeDependent_;
   }; // class Parser

   // constant velocity class is for API compatibility with UBB
   class Velocity : public BoundaryConfiguration {
   public:
             Velocity( const Vector3< real_t > & _velocity ) : velocity_( _velocity ) {}
             Velocity( const real_t _x, const real_t _y, const real_t _z ) : velocity_(_x,_y,_z) {}
      inline Velocity( const Config::BlockHandle & config );

      const Vector3< real_t > & velocity() const { return velocity_; }

      const real_t & x() const { return velocity_[0]; }
      const real_t & y() const { return velocity_[1]; }
      const real_t & z() const { return velocity_[2]; }

      Vector3< real_t > & velocity() { return velocity_; }

      real_t & x() { return velocity_[0]; }
      real_t & y() { return velocity_[1]; }
      real_t & z() { return velocity_[2]; }

   private:
      Vector3< real_t > velocity_;
   }; // class Velocity

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

      if (StoreForce)
         force_->setWithGhostLayer( Vector3<real_t>() );
   }
   void  afterBoundaryTreatment() const {}

   template< typename Buffer_T >
   inline void packCell( Buffer_T &, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const;

   template< typename Buffer_T >
   inline void registerCell( Buffer_T &, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t );

   inline void registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const BoundaryConfiguration & parser );
   inline void registerCells( const flag_t, const CellInterval & cells, const BoundaryConfiguration & parser );
   template< typename CellIterator >
   inline void registerCells( const flag_t, const CellIterator & begin, const CellIterator & end, const BoundaryConfiguration & parser );

   inline void unregisterCell( const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const;

#ifndef NDEBUG
   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask );
#else
   inline void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                               const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ );
#endif

   inline Vector3<real_t> getValue( const cell_idx_t x, cell_idx_t y, cell_idx_t z ) const;
   inline Vector3<real_t> getValue( const cell_idx_t x, cell_idx_t y, cell_idx_t z, real_t t ) const;

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

   shared_ptr<ParserField> parserField_;
   shared_ptr<VelocityField> velocityField_;
   shared_ptr<ForceField> force_;

}; // class ParserUBB

template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce>
inline ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce>::Parser::Parser( const Config::BlockHandle & config )
: parsers_(), equations_(), timeDependent_( false )
{
   if( !config )
      return;

   if( config.isDefined( "x" ) )
   {
      equations_[0] = config.getParameter<std::string>( "x" );
      parsers_[0].parse( equations_[0] );
      if( parsers_[0].symbolExists( "t" ) )
         timeDependent_ = true;
   }
   if( config.isDefined( "y" ) )
   {
      equations_[1] = config.getParameter<std::string>( "y" );
      parsers_[1].parse( equations_[1] );
      if( parsers_[1].symbolExists( "t" ) )
         timeDependent_ = true;
   }
   if( config.isDefined( "z" ) )
   {
      equations_[2] = config.getParameter<std::string>( "z" );
      parsers_[2].parse( equations_[2] );
      if( parsers_[2].symbolExists( "t" ) )
         timeDependent_ = true;
   }
}

template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce>
inline ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce>::Parser::Parser( const std::array< std::string, 3 > & equations )
: parsers_(), equations_( equations ), timeDependent_( false )
{
   if( equations_[0].length() > 0 )
   {
      parsers_[0].parse( equations_[0] );
      if( parsers_[0].symbolExists( "t" ) )
         timeDependent_ = true;
   }
   if( equations_[1].length() > 0 )
   {
      parsers_[1].parse( equations_[1] );
      if( parsers_[1].symbolExists( "t" ) )
         timeDependent_ = true;
   }
   if( equations_[2].length() > 0 )
   {
      parsers_[2].parse( equations_[2] );
      if( parsers_[2].symbolExists( "t" ) )
         timeDependent_ = true;
   }
}

template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce >
inline ParserUBB< LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce >::Velocity::Velocity( const Config::BlockHandle & config  )
{
   velocity_[0] = ( config && config.isDefined( "x" ) ) ? config.getParameter<real_t>( "x" ) : real_c(0.0);
   velocity_[1] = ( config && config.isDefined( "y" ) ) ? config.getParameter<real_t>( "y" ) : real_c(0.0);
   velocity_[2] = ( config && config.isDefined( "z" ) ) ? config.getParameter<real_t>( "z" ) : real_c(0.0);
}

template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce>
Vector3< real_t > ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce>::Parser::operator()( const Vector3< real_t > & x, const real_t t ) const
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

template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce>
Vector3< real_t > ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce>::Parser::operator()( const Vector3< real_t > & x ) const
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



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce>
inline ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce>::ParserUBB( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField * const pdfField,
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

   if (StoreForce)
      force_ = make_shared<ForceField>( pdfField_->xSize(), pdfField_->ySize(), pdfField_->zSize(), pdfField_->nrOfGhostLayers(), field::zyxf );
}



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce >
template< typename Buffer_T >
inline void ParserUBB< LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce >::packCell( Buffer_T & buffer, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   if( parserField_->get( x, y, z ) )
   {
      auto & eqs = parserField_->get( x, y, z )->equations();
      buffer << true << eqs[0] << eqs[1] << eqs[2];
   }
   else
   {
      buffer << false << velocityField_->get( x, y, z );
   }
}



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce >
template< typename Buffer_T >
inline void ParserUBB< LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce >::registerCell( Buffer_T & buffer, const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   bool isparser;
   buffer >> isparser;
   if( isparser )
   {
      std::array< std::string, 3> eqs;
      buffer >> eqs[0] >> eqs[1] >> eqs[2];

      auto p = make_shared<Parser>(eqs);
      parserField_->get( x, y, z ) = p;
   }
   else
   {
      buffer >> velocityField_->get( x, y, z );
      parserField_->get( x, y, z ) = nullptr;
   }
}



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce>
inline ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce>::ParserUBB( const BoundaryUID & boundaryUID, const FlagUID & uid, PDFField * const pdfField,
                                                                                   FlagField<flag_t> * const flagField, const uint_t level, const AABB & aabb )
   : ParserUBB( boundaryUID, uid, pdfField, flagField, nullptr, level, aabb )
{}



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce >
inline void ParserUBB< LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce >::registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                                             const BoundaryConfiguration & parser )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Parser * >( &parser ), &parser );
   WALBERLA_ASSERT_NOT_NULLPTR( parserField_ );

   if( auto v = dynamic_cast< const Velocity * >( &parser ) )
   {
      velocityField_->get( x, y, z ) = v->velocity();
      parserField_->get( x, y, z ) = nullptr;
      return;
   }

   auto & p = dynamic_cast< const Parser & >( parser );

   if( p.isTimeDependent() )
   {
      parserField_->get( x, y, z ) = make_shared<Parser>( p.equations() );
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



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce >
inline void ParserUBB< LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce >::registerCells( const flag_t, const CellInterval & cells, const BoundaryConfiguration & parser )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Parser * >( &parser ), &parser );
   WALBERLA_ASSERT_NOT_NULLPTR( parserField_ );

   if( auto v = dynamic_cast< const Velocity * >( &parser ) )
   {
      for( auto cell = parserField_->beginSliceXYZ( cells ); cell != parserField_->end(); ++cell )
      {
         velocityField_->get( cell.x(), cell.y(), cell.z() ) = v->velocity();
         *cell = nullptr;
      }
      return;
   }

   auto & p = dynamic_cast< const Parser & >( parser );

   if( p.isTimeDependent() )
   {
      auto shared_p = make_shared<Parser>( p.equations() );
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



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce >
template< typename CellIterator >
inline void ParserUBB< LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce >::registerCells( const flag_t, const CellIterator & begin, const CellIterator & end,
                                                                                              const BoundaryConfiguration & parser )
{
   WALBERLA_ASSERT_EQUAL( dynamic_cast< const Parser * >( &parser ), &parser );
   WALBERLA_ASSERT_NOT_NULLPTR( parserField_ );

   if( auto v = dynamic_cast< const Velocity * >( &parser ) )
   {
      for( auto cell = begin; cell != end; ++cell )
      {
         velocityField_->get( cell.x(), cell.y(), cell.z() ) = v->velocity();
         parserField_->get( cell->x(), cell->y(), cell->z() ) = nullptr;
      }
      return;
   }

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



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce >
inline void ParserUBB< LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce >::unregisterCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   parserField_->get(x,y,z) = nullptr;
}



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce>
#ifndef NDEBUG
inline void ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce>::ParserUBB::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                                        const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )
#else
inline void ParserUBB<LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce>::ParserUBB::treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
                                                                                                        const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t /*mask*/ )
#endif
   {
      WALBERLA_ASSERT_EQUAL( nx, x + cell_idx_c( stencil::cx[ dir ] ) );
      WALBERLA_ASSERT_EQUAL( ny, y + cell_idx_c( stencil::cy[ dir ] ) );
      WALBERLA_ASSERT_EQUAL( nz, z + cell_idx_c( stencil::cz[ dir ] ) );
      WALBERLA_ASSERT_UNEQUAL( mask & this->mask_, numeric_cast<flag_t>(0) );
      WALBERLA_ASSERT_EQUAL( mask & this->mask_, this->mask_ ); // only true if "this->mask_" only contains one single flag, which is the case for the
                                                                // current implementation of this boundary condition (ParserUBB)

      const real_t pdf_old = pdfField_->get( x, y, z, Stencil::idx[dir] );

      Vector3<real_t> velocity;
      if( parserField_->get(nx,ny,nz) )
      {
         WALBERLA_ASSERT_NOT_NULLPTR( getTimeTracker(), "A TimeTracker is needed for time-dependent equations" );
         const Vector3< real_t > pos( origin_[0] + real_c(nx) * dx_[0],
                                      origin_[1] + real_c(ny) * dx_[1],
                                      origin_[2] + real_c(nz) * dx_[2] );
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

      if (StoreForce && pdfField_->isInInnerPart( Cell(x,y,z) ))
      {
         const real_t forceMEM = pdf_old + pdfField_->get( nx, ny, nz, Stencil::invDirIdx(dir) );
         Vector3<real_t> force( real_c( stencil::cx[dir] ) * forceMEM,
                                real_c( stencil::cy[dir] ) * forceMEM,
                                real_c( stencil::cz[dir] ) * forceMEM );
         force_->get( nx, ny, nz ) += force;
      }
   }



template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce >
inline Vector3<real_t> ParserUBB< LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce >::getValue( const cell_idx_t x, cell_idx_t y, cell_idx_t z ) const
{
   return getValue( x, y, z, time_ );
}

template< typename LatticeModel_T, typename flag_t, bool AdaptVelocityToExternalForce, bool StoreForce >
inline Vector3<real_t> ParserUBB< LatticeModel_T, flag_t, AdaptVelocityToExternalForce, StoreForce >::getValue( const cell_idx_t x, cell_idx_t y, cell_idx_t z, real_t t ) const
{
   Vector3<real_t> velocity;
   if( parserField_->get(x,y,z) )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( getTimeTracker(), "A TimeTracker is needed for time-dependent equations" );
      const Vector3< real_t > pos( origin_[0] + real_c(x) * dx_[0],
                                   origin_[1] + real_c(y) * dx_[1],
                                   origin_[2] + real_c(z) * dx_[2] );
      velocity = (*parserField_->get(x,y,z))( pos, t );
   }
   else
   {
      velocity = velocityField_->get(x,y,z);
   }
   return velocity;
}


} // namespace lbm
} // namespace walberla
