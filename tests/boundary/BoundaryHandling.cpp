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
//! \file BoundaryHandling.cpp
//! \ingroup boundary
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "boundary/BoundaryHandling.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/cell/CellSet.h"
#include "core/cell/CellVector.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/MPIManager.h"

#include "field/FlagField.h"
#include "field/GhostLayerField.h"

#include "stencil/D3Q27.h"


namespace walberla {



using flag_t = uint8_t;

using FlagField_T = FlagField<flag_t>;
typedef GhostLayerField< uint_t, stencil::D3Q27::Size > WorkField_T;
typedef GhostLayerField< uint_t, 1 >                    FactorField_T;



/////////////////////////////
// COPY BOUNDARY CONDITION //
/////////////////////////////

class CopyBoundary : public Boundary<flag_t> {

public:

   static shared_ptr<BoundaryConfiguration> createConfiguration( const Config::BlockHandle& /*config*/ )
      { return make_shared<BoundaryConfiguration>(); }

   CopyBoundary( const FlagUID& uid1, const FlagUID& uid2, WorkField_T* const field ) :
      Boundary<flag_t>("CopyBoundary"), uid1_( uid1 ), uid2_( uid2 ), field_( field ), beforeCounter_(0), afterCounter_(0) { WALBERLA_ASSERT_NOT_NULLPTR( field_ ); }

   void pushFlags( std::vector< FlagUID >& uids ) const { uids.push_back( uid1_ ); uids.push_back( uid2_ ); }

   void beforeBoundaryTreatment() { ++beforeCounter_; }
   void  afterBoundaryTreatment() { ++afterCounter_; }

   void registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const BoundaryConfiguration& ) {
      registeredCells_.push_back(x,y,z); }
   void registerCells( const flag_t, const CellInterval& interval, const BoundaryConfiguration& ) {
      for( auto cell = interval.begin(); cell != interval.end(); ++cell ) registeredCells_.push_back( cell->x(), cell->y(), cell->z() ); }
   template< typename CellIterator >
   void registerCells( const flag_t, const CellIterator& begin, const CellIterator& end, const BoundaryConfiguration& ) {
      for( auto cell = begin; cell != end; ++cell ) registeredCells_.push_back( cell->x(), cell->y(), cell->z() ); }

   void unregisterCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) { unregisteredCells_.push_back(x,y,z); }

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
      WALBERLA_ASSERT_EQUAL( mask & this->mask_, mask );

      field_->get( nx, ny, nz, stencil::D3Q27::invDirIdx(dir) ) = field_->get( x, y, z, stencil::D3Q27::idx[dir] );
   }

   // non concept functions

   uint_t getBeforeCounter() const { return beforeCounter_; }
   uint_t getAfterCounter()  const { return afterCounter_; }

   const CellVector& getRegisteredCells()   const { return registeredCells_; }
   const CellVector& getUnregisteredCells() const { return unregisteredCells_; }

private:

   const FlagUID uid1_;
   const FlagUID uid2_;

   WorkField_T* const field_;

   uint_t beforeCounter_;
   uint_t afterCounter_;

   CellVector registeredCells_;
   CellVector unregisteredCells_;

}; // class CopyBoundary



////////////////////////////
// ADD BOUNDARY CONDITION //
////////////////////////////

class AddBoundary : public Boundary<flag_t> {

public:

   class Configuration : public BoundaryConfiguration {
   public:
      Configuration( const uint_t _f ) : f_(_f) {}
      Configuration( const Config::BlockHandle& config ) { f_ = ( config && config.isDefined( "f" ) ) ? config.getParameter<uint_t>( "f" ) : uint_c(0); }
      const uint_t& f() const { return f_; }
            uint_t& f()       { return f_; }
   private:
      uint_t f_;
   };

   static shared_ptr<Configuration> createConfiguration( const Config::BlockHandle& config ) { return make_shared<Configuration>( config ); }



   AddBoundary( const FlagUID& uid, WorkField_T* const field ) : Boundary<flag_t>("AddBoundary"), uid_( uid ), field_( field ) {
      WALBERLA_ASSERT_NOT_NULLPTR( field_ );
      factorField_ = make_shared<FactorField_T>( field_->xSize(), field_->ySize(), field_->zSize(), field_->nrOfGhostLayers(), uint_c(0) );
   }

   void pushFlags( std::vector< FlagUID >& uids ) const { uids.push_back( uid_ ); }

   void beforeBoundaryTreatment() const {}
   void  afterBoundaryTreatment() const {}

   inline void registerCell( const flag_t, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const BoundaryConfiguration& factor ) {
      factorField_->get(x,y,z) = dynamic_cast< const Configuration& >( factor ).f();
   }
   inline void registerCells( const flag_t, const CellInterval& cells, const BoundaryConfiguration& factor ) {
      for( auto cell = cells.begin(); cell != cells.end(); ++cell )
         factorField_->get( cell->x(), cell->y(), cell->z() ) = dynamic_cast< const Configuration& >( factor ).f();
   }
   template< typename CellIterator >
   inline void registerCells( const flag_t, const CellIterator& begin, const CellIterator& end, const BoundaryConfiguration& factor ) {
      for( auto cell = begin; cell != end; ++cell )
         factorField_->get( cell->x(), cell->y(), cell->z() ) = dynamic_cast< const Configuration& >( factor ).f();
   }

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
      WALBERLA_ASSERT_EQUAL( mask & this->mask_, mask );

      field_->get( nx, ny, nz, stencil::D3Q27::invDirIdx(dir) ) = field_->get( x, y, z, stencil::D3Q27::idx[dir] ) + factorField_->get( nx, ny, nz );
   }

   // non concept functions

   const FactorField_T& getFactorField() const { return *factorField_; }

private:

   const FlagUID uid_;

   WorkField_T* const field_;

   shared_ptr<FactorField_T> factorField_;

}; // class AddBoundary



////////////////////////////
// TEST BOUNDARY HANDLING //
////////////////////////////

typedef BoundaryHandling< FlagField_T, stencil::D3Q27, CopyBoundary, AddBoundary > TestBoundaryHandling;



///////////////////////////////
// HARD CODED BOUNDARY SWEEP //
///////////////////////////////

static void boundarySweep( FlagField_T& flagField, WorkField_T& workField, FactorField_T& factorField,
                           const FlagUID& nearBoundary, const FlagUID& copyBoundary1, const FlagUID& copyBoundary2, const FlagUID& addBoundary )
{
   flag_t near  = flagField.getFlag( nearBoundary );
   flag_t copy1 = flagField.getFlag( copyBoundary1 );
   flag_t copy2 = flagField.getFlag( copyBoundary2 );
   flag_t add   = flagField.getFlag(  addBoundary );
   flag_t boundary = numeric_cast< flag_t >( copy1 | copy2 | add );

   for( auto cell = flagField.begin(); cell != flagField.end(); ++cell )
   {
      const auto x = cell.x();
      const auto y = cell.y();
      const auto z = cell.z();

      if( isFlagSet( cell, near ) )
      {
         for( auto d = stencil::D3Q27::begin(); d != stencil::D3Q27::end(); ++d )
         {
            const auto nx = cell_idx_c( x + d.cx() );
            const auto ny = cell_idx_c( y + d.cy() );
            const auto nz = cell_idx_c( z + d.cz() );

            if( isPartOfMaskSet( cell.neighbor(*d), boundary ) )
            {
               if( flagField.isPartOfMaskSet( nx, ny, nz, numeric_cast<flag_t>( copy1 | copy2 ) ) ) // behavior identical to CopyBoundary
               {
                  workField.get( nx, ny, nz, d.toInvIdx() ) = workField.get( x, y, z, d.toIdx() );
               }
               else if( flagField.isFlagSet( nx, ny, nz, add ) ) // behavior identical to AddBoundary
               {
                  workField.get( nx, ny, nz, d.toInvIdx() ) = workField.get( x, y, z, d.toIdx() ) + factorField.get( nx, ny, nz );
               }
            }
         }
      }
   }
}



static void resetNearFlags( FlagField_T& flagField, const FlagUID& domainFlag1, const FlagUID& domainFlag2,
                            const FlagUID& nearBoundary, const FlagUID& copyBoundary1, const FlagUID& copyBoundary2, const FlagUID& addBoundary )
{
   flag_t domain1 = flagField.getFlag( domainFlag1 );
   flag_t domain2 = flagField.getFlag( domainFlag2 );
   flag_t near    = flagField.getFlag( nearBoundary );
   flag_t copy1   = flagField.getFlag( copyBoundary1 );
   flag_t copy2   = flagField.getFlag( copyBoundary2 );
   flag_t add     = flagField.getFlag(  addBoundary );
   flag_t boundary = numeric_cast< flag_t >( copy1 | copy2 | add );

   for( auto cell = flagField.beginWithGhostLayer(); cell != flagField.end(); ++cell )
      removeFlag( cell, near );

   for( auto cell = flagField.begin(); cell != flagField.end(); ++cell ) {
      if( isPartOfMaskSet( cell, numeric_cast< flag_t >( domain1 |  domain2 ) ) ) {
         for( auto d = stencil::D3Q27::beginNoCenter(); d != stencil::D3Q27::end(); ++d ) {
            if( isPartOfMaskSet( cell.neighbor(*d), boundary ) )
            {
               addFlag( cell, near );
               break;
            }
         }
      }
   }
}



static void clear( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, FlagField_T& flagField,
                   const FlagUID& copyBoundary1, const FlagUID& copyBoundary2, CellVector& unregisteredCells )
{
   flag_t copy1 = flagField.getFlag( copyBoundary1 );
   flag_t copy2 = flagField.getFlag( copyBoundary2 );
   flag_t boundary = numeric_cast< flag_t >( copy1 | copy2 );

   if( flagField.isPartOfMaskSet(x,y,z,boundary) )
      unregisteredCells.push_back(x,y,z);

   flagField(x,y,z) = 0;
}



static void remove( flag_t mask, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, FlagField_T& flagField,
                    const FlagUID& copyBoundary1, const FlagUID& copyBoundary2, const FlagUID& addBoundary, CellVector& unregisteredCells )
{
   flag_t copy1   = flagField.getFlag( copyBoundary1 );
   flag_t copy2   = flagField.getFlag( copyBoundary2 );
   flag_t add     = flagField.getFlag(  addBoundary );
   flag_t boundary = numeric_cast< flag_t >( copy1 | copy2 | add );

   if( flagField.isPartOfMaskSet(x,y,z,boundary & mask) )
      clear( x, y, z, flagField, copyBoundary1, copyBoundary2, unregisteredCells );
}



//////////
// MAIN //
//////////

static int main( int argc, char **argv )
{
   MPIManager::instance()->initializeMPI( &argc, &argv );

   // INITIALIZATION

   const uint_t xSize = 10;
   const uint_t ySize = 5;
   const uint_t zSize = 8;

   const uint_t gl = 1;

   FlagField_T flagField_Ref( xSize, ySize, zSize, gl );
   FlagField_T flagField_BH( xSize, ySize, zSize, gl );

   FlagUID domain1( "domain1" );
   FlagUID domain2( "domain2" );
   FlagUID add( "add" );
   FlagUID copy1( "copy1" );
   FlagUID copy2( "copy2" );

   flag_t domainFlag1 = flagField_Ref.registerFlag( domain1, uint_c(0) );
                        flagField_BH.registerFlag( domain1, uint_c(0) );
   flag_t domainFlag2 = flagField_Ref.registerFlag( domain2, uint_c(3) );
                        flagField_BH.registerFlag( domain2, uint_c(3) );

   flag_t addFlag = flagField_Ref.registerFlag( add, uint_c(2) );
                    flagField_BH.registerFlag( add, uint_c(2) );

   flag_t copyFlag1 = flagField_Ref.registerFlag( copy1, uint_c(1) );
                      flagField_BH.registerFlag( copy1, uint_c(1) );

   flag_t copyFlag2 = flagField_Ref.registerFlag( copy2, uint_c(4) );
                      flagField_BH.registerFlag( copy2, uint_c(4) );

   WorkField_T workField_Ref( xSize, ySize, zSize, gl );
   WorkField_T workField_BH( xSize, ySize, zSize, gl );

   uint_t value = 0;
   for( cell_idx_t z = -cell_idx_c(gl); z != cell_idx_c(zSize+gl); ++z ) {
      for( cell_idx_t y = -cell_idx_c(gl); y != cell_idx_c(ySize+gl); ++y ) {
         for( cell_idx_t x = -cell_idx_c(gl); x != cell_idx_c(xSize+gl); ++x ) {
            for( uint_t i = 0; i != stencil::D3Q27::Size; ++i, ++value )
            {
               workField_Ref(x,y,z,i) = value;
               workField_BH(x,y,z,i)  = value;
            }
         }
      }
   }

   FactorField_T factorField_Ref( xSize, ySize, zSize, gl, uint_c(0) );

   TestBoundaryHandling handling( "test boundary handling", &flagField_BH, numeric_cast<flag_t>( domainFlag1 | domainFlag2 ),
         CopyBoundary( copy1, copy2, &workField_BH ), AddBoundary( add, &workField_BH ) );

   // START TESTING / COMPARING

   FlagUID near( "near" );
   flag_t nearFlag = flagField_Ref.registerFlag( near );
   WALBERLA_CHECK_EQUAL( nearFlag, handling.getNearBoundaryFlag() ); // TODO: force near flag in flagField_Ref to be identical to near flag of boundary handler

   WALBERLA_CHECK_EQUAL( numeric_cast<flag_t>( 9 ), handling.getDomainMask() );
   WALBERLA_CHECK_EQUAL( numeric_cast<flag_t>( 22 ), handling.getBoundaryMask() );

   for( cell_idx_t z = -cell_idx_c(gl); z != cell_idx_c(zSize+gl); ++z ) {
      for( cell_idx_t y = -cell_idx_c(gl); y != cell_idx_c(ySize+gl); ++y ) {
         for( cell_idx_t x = -cell_idx_c(gl); x != cell_idx_c(xSize+gl); ++x )
         {
            WALBERLA_CHECK( !handling.isBoundary(x,y,z) );
            WALBERLA_CHECK( !handling.isDomain(x,y,z) );
            WALBERLA_CHECK( handling.isEmpty(x,y,z) );
            WALBERLA_CHECK( !handling.isNearBoundary(x,y,z) );
         }
      }
   }

   WALBERLA_CHECK( handling.containsBoundaryCondition( BoundaryUID("CopyBoundary") ) );
   WALBERLA_CHECK( handling.containsBoundaryCondition( BoundaryUID("AddBoundary") ) );
   WALBERLA_CHECK( !handling.containsBoundaryCondition( BoundaryUID("MultiplyBoundary") ) );

   WALBERLA_CHECK( handling.containsBoundaryCondition( copy1 ) );
   WALBERLA_CHECK( handling.containsBoundaryCondition( copy2 ) );
   WALBERLA_CHECK( handling.containsBoundaryCondition( add ) );
   WALBERLA_CHECK( !handling.containsBoundaryCondition( domain1 ) );
   WALBERLA_CHECK( !handling.containsBoundaryCondition( domain2 ) );

   WALBERLA_CHECK( handling.containsBoundaryCondition( copyFlag1 ) );
   WALBERLA_CHECK( handling.containsBoundaryCondition( copyFlag2 ) );
   WALBERLA_CHECK( handling.containsBoundaryCondition( addFlag ) );
   WALBERLA_CHECK( handling.containsBoundaryCondition( numeric_cast<flag_t>( copyFlag1 | copyFlag2 ) ) );
   WALBERLA_CHECK( handling.containsBoundaryCondition( numeric_cast<flag_t>( copyFlag1 | addFlag ) ) );
   WALBERLA_CHECK( handling.containsBoundaryCondition( numeric_cast<flag_t>( copyFlag1 | copyFlag2 | addFlag ) ) );
   WALBERLA_CHECK( !handling.containsBoundaryCondition( 1 ) );
   WALBERLA_CHECK( !handling.containsBoundaryCondition( 8 ) );
   WALBERLA_CHECK( !handling.containsBoundaryCondition( 9 ) );

   WALBERLA_CHECK_EQUAL( handling.numberOfMatchingBoundaryConditions(1), 0 );
   WALBERLA_CHECK_EQUAL( handling.numberOfMatchingBoundaryConditions(2), 1 );
   WALBERLA_CHECK_EQUAL( handling.numberOfMatchingBoundaryConditions(4), 1 );
   WALBERLA_CHECK_EQUAL( handling.numberOfMatchingBoundaryConditions(8), 0 );
   WALBERLA_CHECK_EQUAL( handling.numberOfMatchingBoundaryConditions(16), 1 );
   WALBERLA_CHECK_EQUAL( handling.numberOfMatchingBoundaryConditions(6), 2 );
   WALBERLA_CHECK_EQUAL( handling.numberOfMatchingBoundaryConditions(18), 1 );
   WALBERLA_CHECK_EQUAL( handling.numberOfMatchingBoundaryConditions(20), 2 );
   WALBERLA_CHECK_EQUAL( handling.numberOfMatchingBoundaryConditions(22), 2 );
   WALBERLA_CHECK_EQUAL( handling.numberOfMatchingBoundaryConditions(32), 0 );
   WALBERLA_CHECK_EQUAL( handling.numberOfMatchingBoundaryConditions(96), 0 );
   WALBERLA_CHECK_EQUAL( handling.numberOfMatchingBoundaryConditions(38), 2 );

   WALBERLA_CHECK_EQUAL( handling.getBoundaryMask( "CopyBoundary" ), numeric_cast<flag_t>( copyFlag1 | copyFlag2 ) );
   WALBERLA_CHECK_EQUAL( handling.getBoundaryMask( "AddBoundary" ), addFlag );

   // setDomain - x,y,z

   Cell cell(0,0,5);
   flagField_Ref.addMask( cell.x(), cell.y(), cell.z(), domainFlag1 );
   handling.setDomain( domainFlag1, cell.x(), cell.y(), cell.z() );
   cell = Cell(0,0,8);
   flagField_Ref.addFlag( cell.x(), cell.y(), cell.z(), domainFlag1 );
   handling.setDomain( domainFlag1, cell.x(), cell.y(), cell.z() );
   cell = Cell(-1,0,2);
   flagField_Ref.addFlag( cell.x(), cell.y(), cell.z(), domainFlag2 );
   handling.setDomain( domainFlag2, cell.x(), cell.y(), cell.z() );

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );

   // setDomain - iterators

   CellVector cellVec;
   for( uint_t i = 3; i != 8; ++i )
      cellVec.push_back( cell_idx_c(1), cell_idx_c(5), cell_idx_c(i) );
   for( auto c = cellVec.begin(); c != cellVec.end(); ++c )
      flagField_Ref.addMask( c->x(), c->y(), c->z(), domainFlag1 );
   handling.setDomain( domainFlag1, cellVec.begin(), cellVec.end() );

   cellVec.clear();
   for( uint_t i = 1; i != 5; ++i )
      cellVec.push_back( cell_idx_c(1), cell_idx_c(3), cell_idx_c(i) );
   for( auto c = cellVec.begin(); c != cellVec.end(); ++c )
      flagField_Ref.addFlag( c->x(), c->y(), c->z(), domainFlag1 );
   handling.setDomain( domainFlag1, cellVec.begin(), cellVec.end() );

   CellSet cellSet;
   for( cell_idx_t i = -1; i != 4; ++i )
      cellSet.insert( cell_idx_c(1), cell_idx_c(1), i );
   for( auto c = cellSet.begin(); c != cellSet.end(); ++c )
      flagField_Ref.addFlag( c->x(), c->y(), c->z(), domainFlag2 );
   handling.setDomain( domainFlag2, cellSet.begin(), cellSet.end() );

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );

   // setDomain - interval

   CellInterval cells( cell_idx_c(-1), cell_idx_c(2), cell_idx_c(2), cell_idx_c(-1), cell_idx_c(4), cell_idx_c(6) );
   for( auto c = cells.begin(); c != cells.end(); ++c )
      flagField_Ref.addMask( c->x(), c->y(), c->z(), domainFlag1 );
   handling.setDomain( domainFlag1, cells );

   cells = CellInterval( cell_idx_c(2), cell_idx_c(-1), cell_idx_c(-1), cell_idx_c(3), cell_idx_c(5), cell_idx_c(8) );
   for( auto c = cells.begin(); c != cells.end(); ++c )
      flagField_Ref.addFlag( c->x(), c->y(), c->z(), domainFlag2 );
   handling.setDomain( domainFlag2, cells );

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );

   // getBoundaryCondition

   CopyBoundary& copyBC = handling.getBoundaryCondition<CopyBoundary>( "CopyBoundary" );
   AddBoundary&  addBC  = handling.getBoundaryCondition<AddBoundary>( "AddBoundary" );

   const auto & copy1BCUid = handling.getBoundaryUID( copy1 );
   WALBERLA_CHECK_EQUAL( copy1BCUid.getIdentifier(), "CopyBoundary" );

   const auto & copy2BCUid = handling.getBoundaryUID( copy2 );
   WALBERLA_CHECK_EQUAL( copy2BCUid.getIdentifier(), "CopyBoundary" );

   const auto & addBCUid = handling.getBoundaryUID( add );
   WALBERLA_CHECK_EQUAL( addBCUid.getIdentifier(), "AddBoundary" );

   CellVector registeredCells;
   CellVector unregisteredCells;

   // setBoundary - x,y,z

   cell = Cell(4,0,0);
   flagField_Ref.addFlag( cell.x(), cell.y(), cell.z(), addFlag );
   factorField_Ref( cell.x(), cell.y(), cell.z() ) = uint_c(23);
   handling.setBoundary( add, cell.x(), cell.y(), cell.z(), AddBoundary::Configuration( uint_c(23) ) );

   cell = Cell(4,0,3);
   flagField_Ref.addFlag( cell.x(), cell.y(), cell.z(), copyFlag1 );
   registeredCells.push_back( cell.x(), cell.y(), cell.z() );
   handling.setBoundary( copy1, cell.x(), cell.y(), cell.z() );

   cell = Cell(4,0,5);
   flagField_Ref.addFlag( cell.x(), cell.y(), cell.z(), copyFlag2 );
   registeredCells.push_back( cell.x(), cell.y(), cell.z() );
   handling.setBoundary( copyFlag2, cell.x(), cell.y(), cell.z() );

   // setBoundary - interval

   cells = CellInterval( cell_idx_c(4), cell_idx_c(1), cell_idx_c(-1), cell_idx_c(4), cell_idx_c(2), cell_idx_c(8) );
   for( auto c = cells.begin(); c != cells.end(); ++c )
   {
      flagField_Ref.addFlag( c->x(), c->y(), c->z(), addFlag );
      factorField_Ref( c->x(), c->y(), c->z() ) = uint_c(42);
   }
   handling.setBoundary( add, cells, AddBoundary::Configuration( uint_c(42) ) );

   cells = CellInterval( cell_idx_c(4), cell_idx_c(3), cell_idx_c(-1), cell_idx_c(4), cell_idx_c(5), cell_idx_c(8) );
   for( auto c = cells.begin(); c != cells.end(); ++c )
   {
      flagField_Ref.addMask( c->x(), c->y(), c->z(), copyFlag1 );
      registeredCells.push_back( c->x(), c->y(), c->z() );
   }
   handling.setBoundary( copyFlag1, cells );

   // setBoundary - iterators

   cellVec.clear();
   for( uint_t i = 1; i != 4; ++i )
      cellVec.push_back( cell_idx_c(6), cell_idx_c(2), cell_idx_c(i) );
   for( auto c = cellVec.begin(); c != cellVec.end(); ++c )
   {
      flagField_Ref.addFlag( c->x(), c->y(), c->z(), addFlag );
      factorField_Ref( c->x(), c->y(), c->z() ) = uint_c(5);
   }
   handling.setBoundary( add, cellVec.begin(), cellVec.end(), AddBoundary::Configuration( uint_c(5) ) );

   cellSet.clear();
   for( uint_t i = 1; i != 4; ++i )
      cellSet.insert( cell_idx_c(6), cell_idx_c(3), cell_idx_c(i) );
   for( auto c = cellSet.begin(); c != cellSet.end(); ++c )
   {
      flagField_Ref.addFlag( c->x(), c->y(), c->z(), copyFlag2 );
      registeredCells.push_back( c->x(), c->y(), c->z() );
   }
   handling.setBoundary( copyFlag2, cellSet.begin(), cellSet.end() );

   resetNearFlags( flagField_Ref, domain1, domain2, near, copy1, copy2, add );

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );
   WALBERLA_CHECK( CellSet( registeredCells ) == CellSet( copyBC.getRegisteredCells() ) );
   WALBERLA_CHECK_EQUAL( factorField_Ref, addBC.getFactorField() );

   // fillWithDomain

   cells = CellInterval( cell_idx_c(-1), cell_idx_c(-1), cell_idx_c(-1), cell_idx_c(5), cell_idx_c(5), cell_idx_c(8) );
   for( cell_idx_t z = cells.zMin(); z <= cells.zMax(); ++z ) {
      for( cell_idx_t y = cells.yMin(); y <= cells.yMax(); ++y ) {
         for( cell_idx_t x = cells.xMin(); x <= cells.xMax(); ++x )
         {
            if( flagField_Ref(x,y,z) == 0 )
               flagField_Ref.addMask( x, y, z, domainFlag2 );
         }
      }
   }
   handling.fillWithDomain( domainFlag2, cells );

   resetNearFlags( flagField_Ref, domain1, domain2, near, copy1, copy2, add );

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );

   // perform boundary handling

   boundarySweep( flagField_Ref, workField_Ref, factorField_Ref, near, copy1, copy2, add );
   cells = CellInterval( cell_idx_c(0), cell_idx_c(0), cell_idx_c(0), cell_idx_c(9), cell_idx_c(4), cell_idx_c(7) );
   handling( cells );

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );
   WALBERLA_CHECK_EQUAL( workField_Ref, workField_BH );
   WALBERLA_CHECK_EQUAL( factorField_Ref, addBC.getFactorField() );

   boundarySweep( flagField_Ref, workField_Ref, factorField_Ref, near, copy1, copy2, add );
   handling();

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );
   WALBERLA_CHECK_EQUAL( workField_Ref, workField_BH );
   WALBERLA_CHECK_EQUAL( factorField_Ref, addBC.getFactorField() );

   WALBERLA_CHECK_EQUAL( copyBC.getBeforeCounter(), 1 );
   WALBERLA_CHECK_EQUAL( copyBC.getAfterCounter(), 1 );

   // fillWithDomain

   cells = CellInterval( cell_idx_c(5), cell_idx_c(-1), cell_idx_c(-1), cell_idx_c(7), cell_idx_c(5), cell_idx_c(8) );
   for( cell_idx_t z = cells.zMin(); z <= cells.zMax(); ++z ) {
      for( cell_idx_t y = cells.yMin(); y <= cells.yMax(); ++y ) {
         for( cell_idx_t x = cells.xMin(); x <= cells.xMax(); ++x )
         {
            if( flagField_Ref(x,y,z) == 0 )
               flagField_Ref.addMask( x, y, z, domainFlag1 );
         }
      }
   }
   handling.fillWithDomain( domainFlag1, cells.begin(), cells.end() );

   resetNearFlags( flagField_Ref, domain1, domain2, near, copy1, copy2, add );

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );

   // perform boundary handling

   boundarySweep( flagField_Ref, workField_Ref, factorField_Ref, near, copy1, copy2, add );
   handling();

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );
   WALBERLA_CHECK_EQUAL( workField_Ref, workField_BH );
   WALBERLA_CHECK_EQUAL( factorField_Ref, addBC.getFactorField() );

   WALBERLA_CHECK_EQUAL( copyBC.getBeforeCounter(), 2 );
   WALBERLA_CHECK_EQUAL( copyBC.getAfterCounter(), 2 );

   // fillWithDomain

   cells = CellInterval( cell_idx_c(0), cell_idx_c(0), cell_idx_c(0), cell_idx_c(9), cell_idx_c(4), cell_idx_c(7) );
   for( cell_idx_t z = cells.zMin(); z <= cells.zMax(); ++z ) {
      for( cell_idx_t y = cells.yMin(); y <= cells.yMax(); ++y ) {
         for( cell_idx_t x = cells.xMin(); x <= cells.xMax(); ++x )
         {
            if( flagField_Ref(x,y,z) == 0 )
               flagField_Ref.addMask( x, y, z, domainFlag2 );
         }
      }
   }
   handling.fillWithDomain( domainFlag2, 0 );

   resetNearFlags( flagField_Ref, domain1, domain2, near, copy1, copy2, add );

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );

   // perform boundary handling

   boundarySweep( flagField_Ref, workField_Ref, factorField_Ref, near, copy1, copy2, add );
   handling();

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );
   WALBERLA_CHECK_EQUAL( workField_Ref, workField_BH );
   WALBERLA_CHECK_EQUAL( factorField_Ref, addBC.getFactorField() );

   WALBERLA_CHECK_EQUAL( copyBC.getBeforeCounter(), 3 );
   WALBERLA_CHECK_EQUAL( copyBC.getAfterCounter(), 3 );

   // forceBoundary - x,y,z

   cell = Cell(4,0,5);
   clear( cell.x(), cell.y(), cell.z(), flagField_Ref, copy1, copy2, unregisteredCells );
   flagField_Ref.addFlag( cell.x(), cell.y(), cell.z(), addFlag );
   factorField_Ref( cell.x(), cell.y(), cell.z() ) = uint_c(5);
   handling.forceBoundary( add, cell.x(), cell.y(), cell.z(), AddBoundary::Configuration( uint_c(5) ) );

   cell = Cell(8,1,4);
   clear( cell.x(), cell.y(), cell.z(), flagField_Ref, copy1, copy2, unregisteredCells );
   flagField_Ref.addFlag( cell.x(), cell.y(), cell.z(), copyFlag1 );
   registeredCells.push_back( cell.x(), cell.y(), cell.z() );
   handling.forceBoundary( copy1, cell.x(), cell.y(), cell.z() );

   cell = Cell(8,1,6);
   clear( cell.x(), cell.y(), cell.z(), flagField_Ref, copy1, copy2, unregisteredCells );
   flagField_Ref.addFlag( cell.x(), cell.y(), cell.z(), copyFlag2 );
   registeredCells.push_back( cell.x(), cell.y(), cell.z() );
   handling.forceBoundary( copyFlag2, cell.x(), cell.y(), cell.z() );

   // forceBoundary - interval

   cells = CellInterval( cell_idx_c(8), cell_idx_c(3), cell_idx_c(-1), cell_idx_c(8), cell_idx_c(3), cell_idx_c(8) );
   for( auto c = cells.begin(); c != cells.end(); ++c )
   {
      clear( c->x(), c->y(), c->z(), flagField_Ref, copy1, copy2, unregisteredCells );
      flagField_Ref.addFlag( c->x(), c->y(), c->z(), addFlag );
      factorField_Ref( c->x(), c->y(), c->z() ) = uint_c(5);
   }
   handling.forceBoundary( add, cells, AddBoundary::Configuration( uint_c(5) ) );

   cells = CellInterval( cell_idx_c(8), cell_idx_c(4), cell_idx_c(4), cell_idx_c(8), cell_idx_c(5), cell_idx_c(8) );
   for( auto c = cells.begin(); c != cells.end(); ++c )
   {
      clear( c->x(), c->y(), c->z(), flagField_Ref, copy1, copy2, unregisteredCells );
      flagField_Ref.addMask( c->x(), c->y(), c->z(), copyFlag1 );
      registeredCells.push_back( c->x(), c->y(), c->z() );
   }
   handling.forceBoundary( copyFlag1, cells );

   // forceBoundary - iterators

   cellVec.clear();
   for( uint_t i = 1; i != 4; ++i )
      cellVec.push_back( cell_idx_c(10), cell_idx_c(i), cell_idx_c(3) );
   for( auto c = cellVec.begin(); c != cellVec.end(); ++c )
   {
      clear( c->x(), c->y(), c->z(), flagField_Ref, copy1, copy2, unregisteredCells );
      flagField_Ref.addFlag( c->x(), c->y(), c->z(), addFlag );
      factorField_Ref( c->x(), c->y(), c->z() ) = uint_c(23);
   }
   handling.forceBoundary( add, cellVec.begin(), cellVec.end(), AddBoundary::Configuration( uint_c(23) ) );

   cellSet.clear();
   for( uint_t i = 1; i != 6; ++i )
      cellSet.insert( cell_idx_c(10), cell_idx_c(2), cell_idx_c(i) );
   for( auto c = cellSet.begin(); c != cellSet.end(); ++c )
   {
      clear( c->x(), c->y(), c->z(), flagField_Ref, copy1, copy2, unregisteredCells );
      flagField_Ref.addFlag( c->x(), c->y(), c->z(), copyFlag1 );
      registeredCells.push_back( c->x(), c->y(), c->z() );
   }
   handling.forceBoundary( copyFlag1, cellSet.begin(), cellSet.end() );

   resetNearFlags( flagField_Ref, domain1, domain2, near, copy1, copy2, add );

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );
   WALBERLA_CHECK( CellSet( registeredCells ) == CellSet( copyBC.getRegisteredCells() ) );
   WALBERLA_CHECK( CellSet( unregisteredCells ) == CellSet( copyBC.getUnregisteredCells() ) );
   WALBERLA_CHECK_EQUAL( factorField_Ref, addBC.getFactorField() );

   // perform boundary handling

   boundarySweep( flagField_Ref, workField_Ref, factorField_Ref, near, copy1, copy2, add );
   handling();

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );
   WALBERLA_CHECK_EQUAL( workField_Ref, workField_BH );
   WALBERLA_CHECK_EQUAL( factorField_Ref, addBC.getFactorField() );

   WALBERLA_CHECK_EQUAL( copyBC.getBeforeCounter(), 4 );
   WALBERLA_CHECK_EQUAL( copyBC.getAfterCounter(), 4 );

   // removeBoundary - x,y,z

   cell = Cell(0,0,0);
   remove( numeric_cast<flag_t>( copyFlag1 | copyFlag2 | addFlag ), cell.x(), cell.y(), cell.z(), flagField_Ref, copy1, copy2, add, unregisteredCells );
   handling.removeBoundary( cell.x(), cell.y(), cell.z() );

   cell = Cell(4,0,0);
   remove( numeric_cast<flag_t>( copyFlag1 | copyFlag2 | addFlag ), cell.x(), cell.y(), cell.z(), flagField_Ref, copy1, copy2, add, unregisteredCells );
   handling.removeBoundary( cell.x(), cell.y(), cell.z() );

   cell = Cell(4,0,3);
   remove( numeric_cast<flag_t>( copyFlag1 | copyFlag2 | addFlag ), cell.x(), cell.y(), cell.z(), flagField_Ref, copy1, copy2, add, unregisteredCells );
   handling.removeBoundary( cell.x(), cell.y(), cell.z() );

   cell = Cell(4,0,5);
   remove( copyFlag1, cell.x(), cell.y(), cell.z(), flagField_Ref, copy1, copy2, add, unregisteredCells );
   handling.removeBoundary( copy1, cell.x(), cell.y(), cell.z() );

   // removeBoundary - interval

   cells = CellInterval( cell_idx_c(4), cell_idx_c(-1), cell_idx_c(-1), cell_idx_c(4), cell_idx_c(2), cell_idx_c(5) );
   for( auto c = cells.begin(); c != cells.end(); ++c )
      remove( numeric_cast<flag_t>( copyFlag1 | copyFlag2 | addFlag ), c->x(), c->y(), c->z(), flagField_Ref, copy1, copy2, add, unregisteredCells );
   handling.removeBoundary( cells );

   cells = CellInterval( cell_idx_c(4), cell_idx_c(3), cell_idx_c(4), cell_idx_c(4), cell_idx_c(5), cell_idx_c(8) );
   for( auto c = cells.begin(); c != cells.end(); ++c )
      remove( numeric_cast<flag_t>( copyFlag1 | copyFlag2 ), c->x(), c->y(), c->z(), flagField_Ref, copy1, copy2, add, unregisteredCells );
   handling.removeBoundary( numeric_cast<flag_t>( copyFlag1 | copyFlag2 ), cells );

   // removeBoundary - iterators

   cellVec.clear();
   for( uint_t i = 1; i != 4; ++i )
      cellVec.push_back( cell_idx_c(6), cell_idx_c(2), cell_idx_c(i) );
   for( auto c = cellVec.begin(); c != cellVec.end(); ++c )
      remove( addFlag, c->x(), c->y(), c->z(), flagField_Ref, copy1, copy2, add, unregisteredCells );
   handling.removeBoundary( add, cellVec.begin(), cellVec.end() );

   cellSet.clear();
   for( uint_t i = 1; i != 4; ++i )
      cellSet.insert( cell_idx_c(6), cell_idx_c(3), cell_idx_c(i) );
   for( auto c = cellSet.begin(); c != cellSet.end(); ++c )
      remove( numeric_cast<flag_t>( copyFlag1 | copyFlag2 | addFlag ), c->x(), c->y(), c->z(), flagField_Ref, copy1, copy2, add, unregisteredCells );
   handling.removeBoundary( cellSet.begin(), cellSet.end() );

   resetNearFlags( flagField_Ref, domain1, domain2, near, copy1, copy2, add );

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );
   WALBERLA_CHECK( CellSet( registeredCells ) == CellSet( copyBC.getRegisteredCells() ) );
   WALBERLA_CHECK( CellSet( unregisteredCells ) == CellSet( copyBC.getUnregisteredCells() ) );
   WALBERLA_CHECK_EQUAL( factorField_Ref, addBC.getFactorField() );

   // perform boundary handling

   boundarySweep( flagField_Ref, workField_Ref, factorField_Ref, near, copy1, copy2, add );
   handling();

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );
   WALBERLA_CHECK_EQUAL( workField_Ref, workField_BH );
   WALBERLA_CHECK_EQUAL( factorField_Ref, addBC.getFactorField() );

   WALBERLA_CHECK_EQUAL( copyBC.getBeforeCounter(), 5 );
   WALBERLA_CHECK_EQUAL( copyBC.getAfterCounter(), 5 );

   // clear

   cell = Cell(8,1,4);
   clear( cell.x(), cell.y(), cell.z(), flagField_Ref, copy1, copy2, unregisteredCells );
   handling.clear( cell.x(), cell.y(), cell.z() );

   cell = Cell(8,1,5);
   clear( cell.x(), cell.y(), cell.z(), flagField_Ref, copy1, copy2, unregisteredCells );
   handling.clear( cell.x(), cell.y(), cell.z() );

   cell = Cell(8,1,6);
   clear( cell.x(), cell.y(), cell.z(), flagField_Ref, copy1, copy2, unregisteredCells );
   handling.clear( cell.x(), cell.y(), cell.z() );

   cells = CellInterval( cell_idx_c(8), cell_idx_c(3), cell_idx_c(4), cell_idx_c(8), cell_idx_c(3), cell_idx_c(8) );
   for( auto c = cells.begin(); c != cells.end(); ++c )
      clear( c->x(), c->y(), c->z(), flagField_Ref, copy1, copy2, unregisteredCells );
   handling.clear( cells );

   cells = CellInterval( cell_idx_c(8), cell_idx_c(5), cell_idx_c(4), cell_idx_c(8), cell_idx_c(5), cell_idx_c(8) );
   for( auto c = cells.begin(); c != cells.end(); ++c )
      clear( c->x(), c->y(), c->z(), flagField_Ref, copy1, copy2, unregisteredCells );
   handling.clear( cells );

   cellVec.clear();
   for( uint_t i = 0; i != 3; ++i )
      cellVec.push_back( cell_idx_c(10), cell_idx_c(i), cell_idx_c(3) );
   for( auto c = cellVec.begin(); c != cellVec.end(); ++c )
      clear( c->x(), c->y(), c->z(), flagField_Ref, copy1, copy2, unregisteredCells );
   handling.clear( cellVec.begin(), cellVec.end() );

   cellSet.clear();
   for( uint_t i = 0; i != 4; ++i )
      cellSet.insert( cell_idx_c(10), cell_idx_c(2), cell_idx_c(i) );
   for( auto c = cellSet.begin(); c != cellSet.end(); ++c )
      clear( c->x(), c->y(), c->z(), flagField_Ref, copy1, copy2, unregisteredCells );
   handling.clear( cellSet.begin(), cellSet.end() );

   resetNearFlags( flagField_Ref, domain1, domain2, near, copy1, copy2, add );

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );
   WALBERLA_CHECK( CellSet( registeredCells ) == CellSet( copyBC.getRegisteredCells() ) );
   WALBERLA_CHECK( CellSet( unregisteredCells ) == CellSet( copyBC.getUnregisteredCells() ) );
   WALBERLA_CHECK_EQUAL( factorField_Ref, addBC.getFactorField() );

   // perform boundary handling

   boundarySweep( flagField_Ref, workField_Ref, factorField_Ref, near, copy1, copy2, add );
   handling();

   WALBERLA_CHECK_EQUAL( flagField_Ref, flagField_BH );
   WALBERLA_CHECK_EQUAL( workField_Ref, workField_BH );
   WALBERLA_CHECK_EQUAL( factorField_Ref, addBC.getFactorField() );

   WALBERLA_CHECK_EQUAL( copyBC.getBeforeCounter(), 6 );
   WALBERLA_CHECK_EQUAL( copyBC.getAfterCounter(), 6 );

   // clear

   handling.clear();
   for( auto c = flagField_BH.begin(); c != flagField_BH.end(); ++c )
      WALBERLA_CHECK_EQUAL( *c, numeric_cast<flag_t>(0) );

   handling.clear( gl );
   for( auto c = flagField_BH.beginWithGhostLayer(); c != flagField_BH.end(); ++c )
      WALBERLA_CHECK_EQUAL( *c, numeric_cast<flag_t>(0) );

   WALBERLA_CHECK_EQUAL( copyBC.getRegisteredCells().size(), copyBC.getUnregisteredCells().size() );

   // INFO

   std::cout << "number of   registered boundary cells for CopyBoundary: " << copyBC.getRegisteredCells().size() << std::endl;
   std::cout << "number of unregistered boundary cells for CopyBoundary: " << copyBC.getUnregisteredCells().size() << std::endl;

   std::cout << handling << std::endl;

   return 0;
}



} // namespace walberla



int main(int argc, char **argv) {

   walberla::debug::enterTestMode();

   return walberla::main(argc,argv);
}
