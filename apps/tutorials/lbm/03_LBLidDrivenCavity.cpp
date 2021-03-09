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
//! \file LBLidDrivenCavity.cpp
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/BoundaryHandling.h"

#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"

#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/SimpleUBB.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/vtk/Density.h"
#include "lbm/vtk/Velocity.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/VTKOutput.h"

#include <cstdlib>


///////////
// USING //
///////////

namespace walberla {

//////////////
// TYPEDEFS //
//////////////

using LatticeModel_T = lbm::D3Q19<lbm::collision_model::SRT>;         // the LB lattice model - here: D3Q19, SRT, incompressible, no additional forces
using Stencil_T = LatticeModel_T::Stencil;              // just the D3Q19 stencil without LB specific information
using CommunicationStencil_T = LatticeModel_T::CommunicationStencil; // the stencil that is needed for the communication: This stencil determines which
                                                                         // neighbor blocks are involved during the communication.

using PdfField_T = lbm::PdfField<LatticeModel_T>; // type of the PDF field that stores the 19 distribution functions

using flag_t = walberla::uint8_t;      // each flag consists of an 8 bit value and therefore can distinguish between 8 different markers/flags
using FlagField_T = FlagField<flag_t>; // the flag field: used for marking cells as fluid or obstacle cells
                                          // (also used for distinguishing between different boundary conditions
                                          //  -> every boundary condition possesses its own, unique flag)

using NoSlip_T = lbm::NoSlip<LatticeModel_T, flag_t>; // no slip boundary condition
using UBB_T = lbm::SimpleUBB<LatticeModel_T, flag_t>;    // velocity bounce back boundary condition that internally works with one ...
                                                            // ... constant velocity that must be set during the setup phase

using BoundaryHandling_T = BoundaryHandling<FlagField_T, Stencil_T, NoSlip_T, UBB_T>; // the boundary handling, includes a collection of all boundary conditions

///////////
// FLAGS //
///////////

const FlagUID  Fluid_Flag( "fluid" );                // flag used for marking fluid cells as being fluid cells
const FlagUID    UBB_Flag( "velocity bounce back" ); // flag used for marking all cells at the top as being an obstacle (boundary condition = velocity bounce back)
const FlagUID NoSlip_Flag( "no slip" );              // flag used for marking all remaining domain borders as being an obstacle (boundary condition = no slip)



///////////////////////
// BOUNDARY HANDLING //
///////////////////////

//**********************************************************************************************************************
/*!
*   \brief Functor that is used to add a boundary handling object to each block
*
*   The member function "BoundaryHandling_T * operator()( IBlock* const block ) const" is called by the framework in
*   order to add a boundary handling object of type 'BoundaryHandling_T' to each block.
*/
//**********************************************************************************************************************

class MyBoundaryHandling
{
public:

   MyBoundaryHandling( const BlockDataID & flagFieldId, const BlockDataID & pdfFieldId, const real_t topVelocity ) :
      flagFieldId_( flagFieldId ), pdfFieldId_( pdfFieldId ), topVelocity_( topVelocity ) {}

   BoundaryHandling_T * operator()( IBlock* const block ) const;

private:

   const BlockDataID flagFieldId_;
   const BlockDataID  pdfFieldId_;

   const real_t topVelocity_;

}; // class MyBoundaryHandling

BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block ) const
{
   // Retrieve the pdf and flag field (they must be added to the block structure before the boundary handling is added!)

   FlagField_T * flagField = block->getData< FlagField_T >( flagFieldId_ );
   PdfField_T *   pdfField = block->getData< PdfField_T > (  pdfFieldId_ );

   const flag_t fluid = flagField->registerFlag( Fluid_Flag ); // register the fluid flag at the flag field

   // A new boundary handling instance that uses the just fetched flag field is created:
   // Additional, internal flags used by the boundary handling will be stored in this flag field.
   // The fluid flag is used in order to decide whether a certain cell is part of the domain (and therefore has to be
   // treated by the boundary handling if it is located next to a boundary). Also, the no slip and velocity bounce back
   // boundary conditions are initialized (Please note that the order in which these boundary conditions are initialized
   // must be identical to the order of the template arguments of 'BoundaryHandling_T'!).

   return new BoundaryHandling_T( "boundary handling", flagField, fluid,
                                  NoSlip_T( "no slip", NoSlip_Flag, pdfField ),
                                     UBB_T( "velocity bounce back", UBB_Flag, pdfField, topVelocity_, real_c(0), real_c(0) ) );
}



//**********************************************************************************************************************
/*!
*   \brief Initialization function used for setting flags in the flag field
*
*   Since the goal of this simulation is a lid driven cavity, all cells at the domain border need to be set to no slip,
*   except for the cells at the top: they have to be marked with the velocity bounce back flag.
*   Additionally, all the remaining cells have to be marked as fluid cells.
*/
//**********************************************************************************************************************

void setFlags( shared_ptr< StructuredBlockForest > & blocks, const BlockDataID & boundaryHandlingId )
{
   // First, let's retrieve a cell interval that spans the entire simulation domain. This cell interval stores
   // global cell coordinates: each cell anywhere in the simulation space possesses a unique global cell coordinate!
   // Please note that both cell interval boundaries (to lower as well as the upper boundary) are considered to be part
   // of the cell interval. This means the cell returned by calling min() as wall as the cell returned by calling max()
   // are part of the simulation domain.

   const CellInterval & domainBBInGlobalCellCoordinates = blocks->getDomainCellBB();

   // Iterate through all blocks that are allocated on this process and part of the block structure 'blocks'
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      // In order to mark cells as either domain or boundary, always do so by using member functions of the boundary handling.
      // Never ever set those flags directly in the flag field! Flags that are treated by the boundary handling must never be
      // set in the flag field manually.

      BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingId );

      // After retrieving the boundary handling of the block, let's convert the bounding box of the simulation domain from
      // global to block local cell coordinates.
      // When calling member functions of the boundary handling, always block local cell coordinates are expected!

      CellInterval domainBB( domainBBInGlobalCellCoordinates );
      blocks->transformGlobalToBlockLocalCellInterval( domainBB, *block );

      // Now let's create a cell interval that spans the western boundary of the simulation domain. In y- and z-direction, we
      // leave the initial cell interval unchanged. In x-direction, the newly created interval spans only one single cell,
      // therefore the x-min as wall as the x-max value are both set to the x-min value of the domain bounding box.
      // Note that not all cells of this interval have to be within the current block. Depending on the location of the block
      // in the block structure, maybe no cell of the interval 'west' intersects with the current block! This is fine, those
      // cells are just ignored - you do not have to take care about that.
      // Also note that we use the function "forceBoundary". This function ensures that after completion, all cells are marked
      // with the boundary that was specified via the corresponding flag (here: no slip). All boundary/domain flags that were
      // already set for those cells are safely overwritten by this call.

      CellInterval west( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMin(), domainBB.yMax(), domainBB.zMax() );
      boundaryHandling->forceBoundary( NoSlip_Flag, west );

      // no slip EAST
      CellInterval east( domainBB.xMax(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
      boundaryHandling->forceBoundary( NoSlip_Flag, east );

      // no slip SOUTH
      CellInterval south( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax() );
      boundaryHandling->forceBoundary( NoSlip_Flag, south );

      // no slip NORTH
      CellInterval north( domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
      boundaryHandling->forceBoundary( NoSlip_Flag, north );

      // no slip BOTTOM
      CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
      boundaryHandling->forceBoundary( NoSlip_Flag, bottom );

      // Finally, we mark the entire TOP plane with the velocity bounce back flag.
      // This finalizes our boundary setup.

      CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
      boundaryHandling->forceBoundary( UBB_Flag, top );

      // All the remaining cells need to be marked as being fluid. The 'fillWithDomain' call does just that.

      boundaryHandling->fillWithDomain( domainBB );
   }
}



/////////
// VTK //
/////////

//**********************************************************************************************************************
/*!
*   \brief Helper function for setting up a VTK output instance that generates output for the fluid field
*/
//**********************************************************************************************************************

shared_ptr< vtk::VTKOutput> createFluidFieldVTKWriter( shared_ptr< StructuredBlockForest > & blocks,
                                                       const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId )
{
   // First of all, we have to create a VTK output instance.
   // We configure the output to generate new files only every 20 time steps. We also configure the output to include the
   // ghost layer in the output.

   auto pdfFieldVTKWriter = vtk::createVTKOutput_BlockData( blocks, "fluid_field", uint_c(20), uint_c(1) );

   // Since we want to include the ghost layer in the output, we have to make sure that each cell in the ghost layer contains valid
   // information! For the LB method to work correctly, not all 19 (since we use a D3Q19 stencil) PDF values of a certain cell have to
   // be copied during the communication. For D3Q19, not even all cells need to be copied: the corner cells in the ghost layer are
   // never touched! However, for the visualization of the VTK files, all cells that have been written to file are considered.
   // As a consequence, before writing to file, we need to updated ALL values in ALL ghost layer cells. We can do so by configuring
   // a function that is called every time right before the output is triggered (-> 'addBeforeFunction').
   // For this purpose, we create a new communication scheme that considers all possible neighbors in 3D space (-> D3Q27 - the "center"
   // direction is always ignored). We serialize the data of our PDF field by using the generic 'FieldPackInfo' and by including all
   // 19 values of each cell in the serialization process [PdfField_T::F_SIZE is equal to 19 in our case].

   blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > pdfGhostLayerSync( blocks );
   pdfGhostLayerSync.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldId ) );
   pdfFieldVTKWriter->addBeforeFunction( pdfGhostLayerSync );

   // We cannot simple write fluid properties for each cell: some cells are boundary/obstacle cells and therefore should not be included
   // in the output. We can do so by adding a filter to our VTK output that tells the VTK output to only include (-> 'addCellInclusionFilter')
   // those cells into the output process which are marked by our fluid filter (-> 'fluidFilter'). The fluid filter is a 'FlagFieldCellFilter'
   // that is set up to filter only those cells marked with our fluid flag (-> 'fluidFlag').

   field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldId );
   fluidFilter.addFlag( Fluid_Flag );
   pdfFieldVTKWriter->addCellInclusionFilter( fluidFilter );

   // We are still left with telling the VTK output to actually output any data. We can do so by calling the 'addCellDataWriter' and
   // registering so-called call data writers. Here we use two writers already provided by the lbm module: one for writing the lattice
   // velocities and one for writing the lattice densities. You may notice the two "float" statements: Here we force the density and all
   // components of the velocity to be written as float (32 bit) values in order to save disk space/IO bandwidth, even if internally the
   // simulation uses double values (which is the default case!).

   auto velocityWriter = make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldId, "VelocityFromPDF" );
   auto  densityWriter = make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldId, "DensityFromPDF" );
   pdfFieldVTKWriter->addCellDataWriter( velocityWriter );
   pdfFieldVTKWriter->addCellDataWriter( densityWriter );

   return pdfFieldVTKWriter;
}



//////////
// MAIN //
//////////

int main( int argc, char ** argv )
{
   // We start our simulation by setting up the MPI environment

   mpi::Environment env( argc, argv );

   // This simulation is intended to work with one, four, or sixteen processes only.
   // If the number of currently active MPI processes does not match, we have to abort.

   int processes( MPIManager::instance()->numProcesses() );
   if( processes != 1 && processes != 4 && processes != 16 )
      WALBERLA_ABORT( "The number of processes must be equal to either 1, 4, or 16!" );

   // Our simulation domain spans 24x6x24 cells.
   // If we have only one process, we create just only blocks that spans the entire domain.
   // If we have four processes, we have 1x1x4 blocks with 24x6x6 cells on each block.
   // If we have sixteen processes, we have 4x1x4 blocks with 6x6x6 cells on each block.

   const uint_t xBlocks = ( processes == 16 ) ? uint_c(4) : uint_c(1);
   const uint_t yBlocks = uint_c(1);
   const uint_t zBlocks = ( processes ==  1 ) ? uint_c(1) : uint_c(4);

   // There are a number of different convenience functions that let us create a structured block storage, or more precisely: a structured
   // block forest (the structured block forest data structure is just one possible implementation of a structured block storage).
   // These setup convenience functions handle the domain partitioning, the parallel setup of the data structure and the initial load
   // balancing (= block process placement).
   // We simply choose the one function that best suits our needs: 'createUniformBlockGrid'. This function creates - as the name already
   // suggests - a uniform domain partitioning into equally sized blocks.

   auto blocks = blockforest::createUniformBlockGrid( xBlocks, yBlocks, zBlocks,  // number of blocks in x-, y-, and z-direction
                                                      uint_c(24) / xBlocks,       // number of cells per block (!) in x-direction
                                                      uint_c( 6),                 // number of cells per block (!) in y-direction
                                                      uint_c(24) / zBlocks,       // number of cells per block (!) in z-direction
                                                      real_c(1),                  // dx
                                                      true );                     // true: place exactly one block on each process

   // Next we have to specify all parameters that are required for our lattice model.
   // In our case, we have to provide only one parameter: 'omega' (used during the relaxation process in the lattice Boltzmann (LB) method).

   const real_t omega = real_c(1.4);
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::SRT( omega ) );

   // Now we can start adding data to our blocks.
   // First we add a field that can store the particle distribution functions for the LB method. For this purpose, there is a specialized
   // field class: the PdfField (derived from a ghost layer field and extended by some functionality that is especially useful for the LB method).
   // In order to add a PdfField we can use the 'addPdfFieldToStorage' function from the 'lbm' module. We just have to provide a name and an
   // instance of our lattice model. The PDF field will be default initialized to equilibrium (without any initial velocity and with an initial
   // lattice density equal to 1).

   BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( blocks, "pdf field", latticeModel );

   // After that, we add a flag field that is used for marking cells as either fluid or obstacle.
   // By using the 'addFlagFieldToStorage' convenience function of the 'field' module, the flag field is already default initialized to zero.

   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   // Next, we add a boundary handling instance to each block. Here we have to use the generic 'addBlockData' function of the block structure.
   // As "initialization function" we provide our 'MyBoundaryHandling' class which acts as a functor: in order to add data (in this case our
   // boundary handling) to each block, the function 'operator()( IBlock* const block )' of our 'MyBoundaryHandling' will be called. For each
   // block, this function creates a new instance of our boundary handling which then can be added to the block by 'addBlockData'.
   // For our boundary handling to work we have to provide a velocity for the lid at the top, the block data ID of the flag field, and the
   // block data ID of the PDF field (cf. constructor of 'MyBoundaryHandling').

   const real_t topVelocity = real_c(0.05);
   BlockDataID boundaryHandlingId = blocks->addBlockData< BoundaryHandling_T >( MyBoundaryHandling( flagFieldId, pdfFieldId, topVelocity ),
                                                                                "boundary handling" );

   // Once each block has its own boundary handling instance, we can call our 'setFlags' function which initializes the geometry of our simulation
   // by marking each cell as either fluid or boundary.

   setFlags( blocks, boundaryHandlingId );

   // Now that all data is added and initialized, we can create a time loop and configure it to run for 201 time steps.
   // [201 time steps since we configure our VTK output to write every 20 time steps. However, the first time step is always written by VTK. After
   // that, 19 are skipped and time step number 21 is again written to file. Then, 19 are skipped and in time step 41 another output is triggered.
   // And so on ... This means in time step 201, output number 11 is created - after that we finish the simulation.]

   const uint_t timeSteps = uint_c(201);
   SweepTimeloop timeloop( blocks->getBlockStorage(), timeSteps );

   // Since we are performing a parallel simulation (at least if we start the simulation with 4 or 16 processes), we need to set up
   // a communication scheme. The 'UniformBufferedScheme' together with the 'PdfFieldPackInfo' make sure that all necessary values/directions
   // of the outermost inner cells of the PDF field of one block are copied to the ghost layer of the neighboring block - thereby guaranteeing a
   // synchronization of all data relevant for the LB algorithm.

   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication( blocks );
   communication.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldId ) );

   // After setting up the communication, we schedule the communication to run right before the boundary handling is scheduled to run.

   timeloop.add()
      << BeforeFunction( communication, "communication" )
      << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingId ), "boundary handling" );

   // The LB algorithm is then scheduled to run right after the boundary handling.
   // For this simulation, we selected the cell-wise LB sweep which, because of providing flag field data, ignores cells not marked as fluid
   // cells and performs the LB stream and collide operation only on cells marked as being fluid.

   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, Fluid_Flag ) ), "LB stream & collide" );

   // Finally, we create a VTK output instance which is set up to generate VTK files that include information about the lattice velocity and
   // lattice density of our lid driven cavity simulation.
   // We schedule this VTK output to run at the end of each time step. Please note that the VTK output is indeed triggered at the end of each
   // time step. However, internally the VTK output is configured to actually generate new output only every 20th time step.

   auto pdfFieldVTKWriter = createFluidFieldVTKWriter( blocks, pdfFieldId, flagFieldId );

   timeloop.addFuncAfterTimeStep( vtk::writeFiles( pdfFieldVTKWriter ), "VTK (fluid field data)" );

   // This concludes the set up of our time loop.
   // Before we start our simulation, we generate two VTK outputs:
   // - The first VTK output visualizes the domain decomposition.
   // - The second VTK output shows the initial state of the flag field right before the simulation starts. However, since the geometry remains
   //   fixed during the entire simulation, the state of the flag field does not change and therefore generating only one VTK output of the flag
   //   field is sufficient.

   vtk::writeDomainDecomposition( blocks );

   field::createVTKOutput< FlagField_T >( flagFieldId, *blocks, "flag_field" )();

   // Finally, we can run our lid driven cavity simulation!

   timeloop.run();

   return EXIT_SUCCESS;
}
}

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
