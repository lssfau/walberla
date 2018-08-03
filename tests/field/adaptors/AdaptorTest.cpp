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
//! \file AdaptorTest.cpp
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "field/AddToStorage.h"
#include "field/adaptors/AdaptorCreators.h"
#include "field/adaptors/GhostLayerFieldAdaptor.h"
#include "field/vtk/VTKWriter.h"

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Vector3.h"

#include "gui/Gui.h"

#include "lbm/field/Adaptors.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/MacroscopicValueCalculation.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"

#include "timeloop/SweepTimeloop.h"


namespace walberla {

using LatticeModel_T = lbm::D3Q19<lbm::collision_model::SRT>;

using PdfField = lbm::PdfField<LatticeModel_T>;
typedef GhostLayerField<real_t,1 >          ScalarField;
typedef GhostLayerField<Vector3<real_t>,1 > VectorField;


template<typename Field_T>
class DoubledValueOfField
{
public:
   using basefield_t = Field_T;
   using basefield_iterator = typename Field_T::const_base_iterator;
   using value_type = typename Field_T::value_type;

   static_assert( basefield_t::F_SIZE >= 2, "Only for fields with F_SIZE >=2 " );

   static const uint_t F_SIZE = basefield_t::F_SIZE - 1 ;


   value_type operator() ( const basefield_t & baseField,
                           cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f = 0 ) const
   {
      return 2 * baseField( x,y,z,f );
   }

   value_type operator() ( const basefield_iterator & it ) const
   {
      return 2 * ( *it );
   }
};


void iteratorTest()
{
   const real_t MAGIC_SRC = 42.0;

   // field has two f's and two ghost layers
   GhostLayerField<real_t, 2> field ( 4, 4, 3, 2, MAGIC_SRC    );

   // adapter reduces field to have only one f and one ghost layer
   using Functor = DoubledValueOfField<decltype(field)>;
   field::GhostLayerFieldAdaptor<Functor, 1 > adaptor ( field );

   uint_t ctr = 0;
   for ( auto it = adaptor.beginGhostLayerOnly( stencil::T ); it != adaptor.end(); ++it ) {
      WALBERLA_CHECK_FLOAT_EQUAL( *it, 2 * MAGIC_SRC );
      ctr++;
   }
   WALBERLA_CHECK_EQUAL( ctr, 4*4 * 1 * 1 ); //xSize*ySize*nrGhostLayer*fSize  of adaptor


   ctr = 0;
   for ( auto it = adaptor.beginSliceBeforeGhostLayer( stencil::T, 2 ); it != adaptor.end(); ++it ) {
      WALBERLA_CHECK_FLOAT_EQUAL( *it, 2 * MAGIC_SRC );
      ctr++;
   }
   WALBERLA_CHECK_EQUAL( ctr, 4*4 * 2 * 1 ); //xSize*ySize*twoSlices*fSize  of adaptor

}


void simpleTest()
{
   typedef GhostLayerField<real_t, 3> MyField;

   MyField field( 2, 2, 1, 2, 42 );

   field::GhostLayerFieldAdaptor<DoubledValueOfField<MyField>, 1 > adaptor ( field );

   WALBERLA_CHECK_EQUAL( adaptor.nrOfGhostLayers(), 1 );
}


int main( int argc, char ** argv )
{
   walberla::Environment env( argc, argv );
   debug::enterTestMode();

   simpleTest();
   iteratorTest();


   using blockforest::createUniformBlockGrid;
   shared_ptr<StructuredBlockForest>
   blocks = createUniformBlockGrid( 1,1,1,
                                    4,4,4,
                                    1,
                                    true,               // one block per process
                                    false,false,false,  // periodicity
                                    false );            // do NOT keep global information


   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::SRT( real_t(1.4) ) );

   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage( blocks, "PdfField", latticeModel );

   BlockDataID densityID  = field::addFieldAdaptor<lbm::Adaptor<LatticeModel_T>::Density>       ( blocks, pdfFieldID, "DensityAdaptor" );
   BlockDataID velocityID = field::addFieldAdaptor<lbm::Adaptor<LatticeModel_T>::VelocityVector>( blocks, pdfFieldID, "VelocityAdaptor" );

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      PdfField * pdfField = blockIt->getData<PdfField> ( pdfFieldID );
      for( auto cellIt = pdfField->beginXYZ(); cellIt != pdfField->end(); ++cellIt )
      {
         Vector3<real_t> vel ( real_c(cellIt.x()), real_c(cellIt.y()), real_c(cellIt.z()) );
         lbm::setDensityAndVelocity( cellIt, latticeModel, vel, real_c( cellIt.x() ) );
      }
   }


   SweepTimeloop timeloop( blocks, 1 );

   auto writeDensity  = field::createVTKOutput<lbm::Adaptor<LatticeModel_T>::Density>         ( densityID,  *blocks, "density"  );
   auto writeVelocity = field::createVTKOutput<lbm::Adaptor<LatticeModel_T>::VelocityVector  >( velocityID, *blocks, "velocity" );
   writeDensity();
   writeVelocity();

   //GUI gui ( timeloop, blocks, argc, argv );
   //gui.run();


   return 0;
}
}

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
