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
//! \file AccuracyEvaluationTest.cpp
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/debug/TestSubsystem.h"
#include "core/math/DistributedSample.h"
#include "core/math/Random.h"
#include "core/mpi/Environment.h"

#include "field/AccuracyEvaluation.h"
#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"

//#include <ctime>



namespace accuracy_evaluation_test {
   
using namespace walberla;

typedef field::GhostLayerField< real_t, 1 >           ScalarField_T;
typedef field::GhostLayerField< Vector3<real_t>, 1 >  VectorField_T;

const real_t          scalarValue( real_c(23) );
const Vector3<real_t> vectorValue( real_c(23), real_c(42), real_c(5) );


real_t scalarSolution( const Vector3<real_t> & )
{
   return scalarValue;
}

Vector3<real_t> vectorSolution( const Vector3<real_t> & )
{
   return vectorValue;
}


int main( int argc, char* argv[] )
{
   debug::enterTestMode();

   mpi::Environment mpiEnv( argc, argv );

   auto blocks = blockforest::createUniformBlockGrid( uint_t( 2), uint_t( 1), uint_t( 2), // blocks
                                                      uint_t(10), uint_t(10), uint_t(10), // cells
                                                      real_t(1), // dx
                                                      uint_t( 2), uint_t( 1), uint_t( 2) ); // number of processes

   //math::seedRandomGenerator( numeric_cast<std::mt19937::result_type>( std::time(0) ) );
   math::seedRandomGenerator( numeric_cast<std::mt19937::result_type>( MPIManager::instance()->rank() ) );

   auto sId = field::addToStorage< ScalarField_T >( blocks, "scalar field" );
   auto vId = field::addToStorage< VectorField_T >( blocks, "vector field" );

   math::DistributedSample sValues;
   math::DistributedSample s2Values;
   math::DistributedSample vValues;
   math::DistributedSample v2Values;

   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      ScalarField_T * sField = block->getData< ScalarField_T >( sId );
      VectorField_T * vField = block->getData< VectorField_T >( vId );

      for( auto cell = sField->begin(); cell != sField->end(); ++cell )
      {
         const real_t value = math::realRandom<real_t>();
         const real_t error = std::abs( value - scalarValue );
         sValues.insert( error );
         s2Values.insert( error * error );
         *cell = value;
      }

      for( auto cell = vField->begin(); cell != vField->end(); ++cell )
      {
         const Vector3<real_t> value( math::realRandom<real_t>(), math::realRandom<real_t>(), math::realRandom<real_t>() );
         const Vector3<real_t> diff = value - vectorValue;
         vValues.insert( diff.length() );
         v2Values.insert( diff.length() * diff.length() );
         *cell = value;
      }
   }

   sValues.mpiGatherRoot();
   s2Values.mpiGatherRoot();
   vValues.mpiGatherRoot();
   v2Values.mpiGatherRoot();

   auto sAccuracy = field::makeAccuracyEvaluation< ScalarField_T >( blocks, sId, &scalarSolution, uint_t(0), uint_t(1) );
   (*sAccuracy)();

   auto vAccuracy = field::makeAccuracyEvaluation< VectorField_T >( blocks, vId, &vectorSolution, uint_t(0), uint_t(1) );
   (*vAccuracy)();

   WALBERLA_ROOT_SECTION()
   {
      WALBERLA_CHECK_FLOAT_EQUAL( sAccuracy->L1(), sValues.sum() / real_c(4000) );
      WALBERLA_CHECK_FLOAT_EQUAL( sAccuracy->L2(), std::sqrt( s2Values.sum() / real_c(4000) ) );
      WALBERLA_CHECK_FLOAT_EQUAL( sAccuracy->Lmax(), sValues.max() );

      WALBERLA_CHECK_FLOAT_EQUAL( vAccuracy->L1(), vValues.sum() / real_c(4000) );
      WALBERLA_CHECK_FLOAT_EQUAL( vAccuracy->L2(), std::sqrt( v2Values.sum() / real_c(4000) ) );
      WALBERLA_CHECK_FLOAT_EQUAL( vAccuracy->Lmax(), vValues.max() );
   }

   return EXIT_SUCCESS;
}

}

int main( int argc, char* argv[] )
{
   return accuracy_evaluation_test::main( argc, argv );
}
