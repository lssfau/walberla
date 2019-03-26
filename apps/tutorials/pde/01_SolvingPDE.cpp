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
//! \file 05_SolvingPDE.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/math/Constants.h"

#include "field/Field.h"
#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "stencil/D2Q5.h"

#include "gui/Gui.h"
#include "timeloop/SweepTimeloop.h"

#include <cmath>
#include <vector>


namespace walberla {

typedef GhostLayerField<real_t,1> ScalarField;
using Stencil_T = stencil::D2Q5;


// function to initialize the boundaries of the source and destination fields
// The values of the Dirichlet boundary conditions are stored in the ghost layers of the outer blocks.
void initBC( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & srcID, const BlockDataID & dstID )
{
   // iterate all blocks with an iterator 'block'
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      // The Dirichlet boundary condition is non-zero only at the top boundary.
      if( blocks->atDomainYMaxBorder( *block ) )
      {
         // get the field data out of the blocks
         auto src = block->getData< ScalarField >( srcID );
         auto dst = block->getData< ScalarField >( dstID );

         // obtain a CellInterval object that holds information about the number of cells in x,y,z direction of the field inlcuding ghost layers
         // Since src and dst have the same size, one object is enough.
         CellInterval xyz = src->xyzSizeWithGhostLayer();

         // To only loop over the top row of the cells, i.e., the ghost layer containing the domain boundary, restrict the range in y-direction.
         xyz.yMin() = xyz.yMax();

         // iterate all cells in the top boundary with the iterator 'cell'
         for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
         {
            // obtain the physical coordinate of the center of the current cell
            const Vector3< real_t > p = blocks->getBlockLocalCellCenter( *block, *cell );

            // set the values of the Dirichlet boundary condition given by the function u(x,1) = sin(2*M_PI*x)*sinh(2*M_PI) in the source and destination field
            src->get( *cell ) = std::sin( real_c(2) * math::M_PI * p[0] ) * std::sinh( real_c(2) * math::M_PI );
            dst->get( *cell ) = std::sin( real_c(2) * math::M_PI * p[0] ) * std::sinh( real_c(2) * math::M_PI );
         }
      }
   }
}



// function to initialize the field holding the right-hand side function f
void initRHS( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & rhsID )
{
   // iterate all blocks with an iterator 'block'
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      // get the field data out of the block
      auto f = block->getData< ScalarField >( rhsID );

      // obtain a CellInterval object that holds information about the number of cells in x,y,z direction of the field excluding ghost layers
      CellInterval xyz = f->xyzSize();

      // iterate all (inner) cells in the field
      for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
      {
         // obtain the physical coordinate of the center of the current cell
         const Vector3< real_t > p = blocks->getBlockLocalCellCenter( *block, *cell );

         // set the right-hand side, given by the function f(x,y) = 4*M_PI*PI*sin(2*M_PI*x)*sinh(2*M_PI*y)
         f->get( *cell ) = real_c(4) * math::M_PI * math::M_PI * std::sin( real_c(2) * math::M_PI * p[0] ) *
                           std::sinh( real_c(2) * math::M_PI * p[1] );
      }
   }
}



// class containing the Jacobi sweep
class JacobiSweep
{
public:

   // constructor
   JacobiSweep( const BlockDataID & srcID, const BlockDataID & dstID, const BlockDataID & rhsID, const real_t & dx, const real_t & dy )
      : srcID_( srcID ), dstID_( dstID ), rhsID_( rhsID ), dx_( dx ), dy_( dy ) {}

   // functor that carries out the Jacobi sweep
   void operator()(IBlock * const block);

private:

   const BlockDataID srcID_;
   const BlockDataID dstID_;
   const BlockDataID rhsID_;

   const real_t dx_;
   const real_t dy_;
};

void JacobiSweep::operator()( IBlock * const block )
{
   // get the source, destination, and rhs field data from the block
   auto src = block->getData< ScalarField >( srcID_ );
   auto dst = block->getData< ScalarField >( dstID_ );
   auto rhs = block->getData< ScalarField >( rhsID_ );

   // iterate all cells (excluding ghost layers) of the fields
   // This macro is given a field (here src) with the matching size.
   // Inside the macro, x, y and z can be used as indices to access the data.
   // Since the three fields have the same size, these indices are valid for all of them.
   WALBERLA_FOR_ALL_CELLS_XYZ(src,
      // carries out the sweep for the current cell (x,y,z) with the respective prefactors
      dst->get(x,y,z) =  rhs->get(x,y,z);
      dst->get(x,y,z) += ( real_c(1) / (dx_ * dx_) ) * src->get( x+1,  y , z );
      dst->get(x,y,z) += ( real_c(1) / (dx_ * dx_) ) * src->get( x-1,  y , z );
      dst->get(x,y,z) += ( real_c(1) / (dy_ * dy_) ) * src->get(  x , y+1, z );
      dst->get(x,y,z) += ( real_c(1) / (dy_ * dy_) ) * src->get(  x , y-1, z );
      dst->get(x,y,z) /= ( real_c(2) / (dx_ * dx_) + real_c(2)/(dy_ * dy_) + real_c(4) * math::M_PI * math::M_PI );
   )

   // swap source and destination fields
   src->swapDataPointers( dst );
}



// class containing the Jacobi sweep using the stencil concept
class JacobiSweepStencil
{
public:

   // constructor
   JacobiSweepStencil( const BlockDataID & srcID, const BlockDataID & dstID, const BlockDataID & rhsID, const std::vector< real_t > & weights )
      : srcID_( srcID ), dstID_( dstID ), rhsID_( rhsID ), weights_( weights )
   {
      // store the center weight in inverse form to avoid divisions in the sweep
      // weights_[ Stencil_T::idx[ stencil::C ] ] = real_c(1) / weights_[ Stencil_T::idx[ stencil::C ] ];
   }

   // functor that carries out the Jacobi sweep
   void operator()( IBlock * block );

private:

   const BlockDataID srcID_;
   const BlockDataID dstID_;
   const BlockDataID rhsID_;

   std::vector< real_t > weights_;
};

void JacobiSweepStencil::operator()( IBlock * block )
{
   // get the source, destination and rhs field data from the block
   auto src = block->getData< ScalarField >( srcID_ );
   auto dst = block->getData< ScalarField >( dstID_ );
   auto rhs = block->getData< ScalarField >( rhsID_ );

   // iterate all cells (excluding ghost layers) of the fields and carry out the Jacobi sweep
   WALBERLA_FOR_ALL_CELLS_XYZ(src,

      dst->get(x,y,z) =  rhs->get(x,y,z);

      // iterate the neighboring cells and multiply their value with the respective weight
      for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir )
         dst->get(x,y,z) -= weights_[ dir.toIdx() ] * src->getNeighbor(x,y,z,*dir);

      dst->get(x,y,z) /= weights_[ Stencil_T::idx[ stencil::C ] ];

      // When the center weight is stored in inverse form, a multiplication instead of a division can be applied here.
      // dst->get(x,y,z) *= weights_[ Stencil_T::idx[ stencil::C ] ];
   )

   // swap source and destination fields
   src->swapDataPointers( dst );
}



int main( int argc, char ** argv )
{
   walberla::Environment env( argc, argv );

   // number of blocks in x,y,z direction
   const uint_t xBlocks = uint_c(1);
   const uint_t yBlocks = uint_c(1);
   const uint_t zBlocks = uint_c(1);

   // number of cells per block in x,y,z direction
   // two dimensional layout, so only one cell in z direction
   const uint_t xCells = uint_c(50);
   const uint_t yCells = uint_c(25);
   const uint_t zCells = uint_c(1);

   // domain size is [0,2]x[0,1]
   const real_t xMin = real_c(0);
   const real_t xMax = real_c(2);
   const real_t yMin = real_c(0);
   const real_t yMax = real_c(1);

   // mesh spacings in x and y direction
   const real_t dx = (xMax - xMin)/real_c( xBlocks*xCells + uint_c(1) );
   const real_t dy = (yMax - yMin)/real_c( yBlocks*yCells + uint_c(1) );

   // create an axis-aligned bounding box to define domain
   // defines a rectangular domain by specifying (xMin, yMin, zMin, xMax, yMax, zMax)
   // care has to be taken regarding the cell-centered data layout in waLBerla
   // in z-direction only one cell layer with arbitrary thickness, here dx, is present
   auto aabb = math::AABB( xMin + real_c(0.5)*dx, yMin + real_c(0.5)*dy, real_c(0),
                           xMax - real_c(0.5)*dx, yMax - real_c(0.5)*dy, dx );

   // create blocks
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid (
             aabb,                               // axis-aligned bounding box
             xBlocks, yBlocks, zBlocks,          // number of blocks in x,y,z direction
             xCells,  yCells,  zCells,           // how many cells per block (x,y,z)
             false,                              // one block per process - "false" means all blocks to one process
             false, false, false );              // no periodicity


   // add fields with ghost layers to all blocks
   // source and destination fields for the unknowns u
   BlockDataID srcID = field::addToStorage< ScalarField >( blocks, "src", real_c(0), field::zyxf, uint_c(1));
   BlockDataID dstID = field::addToStorage< ScalarField >( blocks, "dst", real_c(0), field::zyxf, uint_c(1));
   // field to store the function f
   BlockDataID rhsID = field::addToStorage< ScalarField >( blocks, "rhs", real_c(0), field::zyxf, uint_c(1));

   // initialize the field
   initRHS( blocks, rhsID );

   // set the Dirichlet boundary conditions
   initBC( blocks, srcID, dstID );

   // create the communication scheme for the D2Q5 stencil
   blockforest::communication::UniformBufferedScheme< stencil::D2Q5 > myCommScheme( blocks );
   // add a PackInfo that packs/unpacks the source field
   // Since the RHS field is unchanged, there is no need to communicate it.
   myCommScheme.addPackInfo( make_shared< field::communication::PackInfo<ScalarField> >( srcID ) );

   // set up the timeloop object
   SweepTimeloop timeloop ( blocks, uint_c(1) );

   // register the communication function and the Jacobi sweep in the timeloop
   // either the hard coded version...
   /*
   timeloop.add() << BeforeFunction( myCommScheme, "Communication" )
                  << Sweep( JacobiSweep( srcID, dstID, rhsID, dx, dy ), "JacobiSweep" );
   */

   // ...or the variant using the stencil concept
   // set up the stencil weights
   std::vector< real_t > weights( Stencil_T::Size );
   weights[ Stencil_T::idx[ stencil::C ] ] = real_c(2) / ( dx * dx ) + real_c(2) / ( dy * dy ) + real_c(4) * math::M_PI * math::M_PI;
   weights[ Stencil_T::idx[ stencil::E ] ] = real_c(-1) / ( dx * dx );
   weights[ Stencil_T::idx[ stencil::W ] ] = real_c(-1) / ( dx * dx );
   weights[ Stencil_T::idx[ stencil::N ] ] = real_c(-1) / ( dy * dy );
   weights[ Stencil_T::idx[ stencil::S ] ] = real_c(-1) / ( dy * dy );

   timeloop.add() << BeforeFunction( myCommScheme, "Communication" )
                  << Sweep( JacobiSweepStencil( srcID, dstID, rhsID, weights ), "JacobiSweepStencil" );

   // start the GUI and run the simulation
   GUI gui ( timeloop, blocks, argc, argv );
   gui.run();

   return 0;
}
}

int main( int argc, char ** argv )
{
   walberla::main(argc, argv);
}
