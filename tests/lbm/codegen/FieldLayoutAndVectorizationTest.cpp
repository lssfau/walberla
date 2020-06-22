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
//! \file FieldLayoutAndVectorizationTest.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"
#include "core/logging/all.h"

#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"

#include "FieldLayoutAndVectorizationTest_FZYX_Vec_LatticeModel.h"
#include "FieldLayoutAndVectorizationTest_FZYX_NoVec_LatticeModel.h"
#include "FieldLayoutAndVectorizationTest_ZYXF_Vec_LatticeModel.h"
#include "FieldLayoutAndVectorizationTest_ZYXF_NoVec_LatticeModel.h"

namespace field_layout_vectorization_test {

using namespace walberla;

typedef walberla::uint8_t flag_t;
typedef FlagField<flag_t> FlagField_T;

template<typename LM1_T, typename LM2_T>
void checkEquivalence(const shared_ptr<StructuredBlockStorage> & blocks, BlockDataID pdfField1ID, BlockDataID pdfField2ID, std::string identifier)
{
   for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
      auto pdfField1 = blockIt->getData<lbm::PdfField<LM1_T> >(pdfField1ID);
      auto pdfField2 = blockIt->getData<lbm::PdfField<LM2_T> >(pdfField2ID);

      WALBERLA_FOR_ALL_CELLS_XYZ(pdfField1,
            for( uint_t f = 0; f < pdfField1->fSize(); ++f)
            {
               WALBERLA_CHECK_FLOAT_EQUAL(pdfField1->get(x,y,z,f), pdfField2->get(x,y,z,f),
                                          identifier << " - mismatch in " << x << " " << y << " " << z << " " << f);
            }
            )

   }

}

int main(int argc, char **argv) {

   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   walberla::Environment walberlaEnv(argc, argv);

   auto blocks = blockforest::createUniformBlockGrid( 1, 1, 1,
                                                      32, 32, 32, real_t(1),
                                                      0, false, false,
                                                      true, true, true, //periodicity
                                                      false );


   real_t omega = real_t(1.6);
   Vector3<real_t> initialVel(real_t(0.1), real_t(0), real_t(0));

   using LatticeModel1_T = lbm::FieldLayoutAndVectorizationTest_FZYX_Vec_LatticeModel;
   LatticeModel1_T latticeModel1(omega);
   BlockDataID pdfField1ID = lbm::addPdfFieldToStorage(blocks, "pdf field1", latticeModel1, initialVel, real_t(1), field::fzyx);
   for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
      auto sweep = LatticeModel1_T::Sweep(pdfField1ID);
      sweep(&(*blockIt));
   }

   using LatticeModel2_T = lbm::FieldLayoutAndVectorizationTest_FZYX_NoVec_LatticeModel;
   LatticeModel2_T latticeModel2(omega);
   BlockDataID pdfField2ID = lbm::addPdfFieldToStorage(blocks, "pdf field2", latticeModel2, initialVel, real_t(1), field::fzyx);
   for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
      auto sweep = LatticeModel2_T::Sweep(pdfField2ID);
      sweep(&(*blockIt));
   }

   checkEquivalence<LatticeModel1_T,LatticeModel2_T>(blocks, pdfField1ID, pdfField2ID, "1 vs 2");

   using LatticeModel3_T = lbm::FieldLayoutAndVectorizationTest_ZYXF_Vec_LatticeModel;
   LatticeModel3_T latticeModel3(omega);
   BlockDataID pdfField3ID = lbm::addPdfFieldToStorage(blocks, "pdf field3", latticeModel3, initialVel, real_t(1), field::zyxf);
   for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
      auto sweep = LatticeModel3_T::Sweep(pdfField3ID);
      sweep(&(*blockIt));
   }

   checkEquivalence<LatticeModel2_T,LatticeModel3_T>(blocks, pdfField2ID, pdfField3ID,  "2 vs 3");

   using LatticeModel4_T = lbm::FieldLayoutAndVectorizationTest_ZYXF_NoVec_LatticeModel;
   LatticeModel4_T latticeModel4(omega);
   BlockDataID pdfField4ID = lbm::addPdfFieldToStorage(blocks, "pdf field4", latticeModel4, initialVel, real_t(1), field::zyxf);
   for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt) {
      auto sweep = LatticeModel4_T::Sweep(pdfField4ID);
      sweep(&(*blockIt));
   }

   checkEquivalence<LatticeModel3_T,LatticeModel4_T>(blocks, pdfField3ID, pdfField4ID,  "3 vs 4");

   return EXIT_SUCCESS;
}

} //namespace field_layout_vectorization_test

int main( int argc, char **argv ){
   field_layout_vectorization_test::main(argc, argv);
}
