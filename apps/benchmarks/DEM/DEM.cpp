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
//! \file DEM.cpp
//! \brief demonstration of basic functionality of DEM
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "pe/basic.h"

#include <blockforest/Initialization.h>
#include <core/DataTypes.h>

#include <string>
#include <tuple>

namespace walberla {

namespace dem {
real_t calcCoefficientOfRestitution(const real_t k, const real_t gamma, const real_t meff)
{
   auto a = real_t(0.5) * gamma / meff;
   return std::exp(-a * math::pi / std::sqrt(k / meff - a*a));
}

real_t calcCollisionTime(const real_t k, const real_t gamma, const real_t meff)
{
   auto a = real_t(0.5) * gamma / meff;
   return math::pi / std::sqrt( k/meff - a*a);
}
}

int main( int argc, char** argv )
{
   using namespace walberla;
   using namespace walberla::pe;

   using BodyTuple = std::tuple<Sphere, Plane> ;

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   real_t dt         = real_c(0.0001);                                  //!< integration time
   real_t radius     = real_c(1);                                       //!< particle radius
   real_t density    = real_c(2707);                                    //!< particle density
   real_t m          = Sphere::calcMass( radius, density);              //!< particle mass
   std::string model = "kg";                                            //!< input model
   real_t k          = real_c(8.11e6);                                  //!< linear spring stiffness
   real_t gamma      = real_c(6.86e1);                                  //!< damper
   real_t e          = dem::calcCoefficientOfRestitution(k, gamma, m);  //!< coefficient of restitution
   real_t t          = dem::calcCollisionTime(k, gamma, m);             //!< collision time

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "-dt" )                      == 0 ) dt      = real_c( std::stod( argv[++i] ) );
      else if( std::strcmp( argv[i], "-radius" )             == 0 ) radius  = real_c( std::stod( argv[++i] ) );
      else if( std::strcmp( argv[i], "-density" )            == 0 ) density = real_c( std::stod( argv[++i] ) );
      else if( std::strcmp( argv[i], "-kg" )                 == 0 )
      {
         k     = real_c( std::stod( argv[++i] ) );
         gamma = real_c( std::stod( argv[++i] ) );
         model = "kg";
      }
      else if( std::strcmp( argv[i], "-et" )                 == 0 )
      {
         e = real_c( std::atof( argv[++i] ) );
         t = real_c( std::atof( argv[++i] ) );
         model = "et";
      }
      else WALBERLA_ABORT("Found invalid command line argument: \"" << argv[i] << "\" - aborting...");
   }

   m          = Sphere::calcMass( radius, density);
   if (model=="kg")
   {
      e          = dem::calcCoefficientOfRestitution(k, gamma, m);  //!< coefficient of restitution
      t          = dem::calcCollisionTime(k, gamma, m);             //!< collision time
   } else
      WALBERLA_ABORT("unsupported model");

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   // create blocks
   shared_ptr< StructuredBlockForest > forest = blockforest::createUniformBlockGrid(
                                                   math::AABB(-5,-5,0,5,5,10),
                                                   uint_c( 1), uint_c( 1), uint_c( 1), // number of blocks in x,y,z direction
                                                   uint_c( 1), uint_c( 1), uint_c( 1), // how many cells per block (x,y,z)
                                                   true,                               // max blocks per process
                                                   false, false, false,                // full periodicity
                                                   false);

   SetBodyTypeIDs<BodyTuple>::execute();

   auto storageID           = forest->addBlockData(createStorageDataHandling<BodyTuple>(), "Storage");
   auto hccdID              = forest->addBlockData(ccd::createHashGridsDataHandling( globalBodyStorage, storageID ), "HCCD");
   auto fcdID               = forest->addBlockData(fcd::createGenericFCDDataHandling<BodyTuple, fcd::AnalyticCollideFunctor>(), "FCD");
   cr::DEM cr(globalBodyStorage, forest->getBlockStoragePointer(), storageID, hccdID, fcdID);

   const real_t   static_cof  ( real_t(0.4) / real_t(2) );    // Coefficient of static friction. Roughly 0.85 with high variation depending on surface roughness for low stresses. Note: pe doubles the input coefficient of friction for material-material contacts.
   const real_t   dynamic_cof ( static_cof ); // Coefficient of dynamic friction. Similar to static friction for low speed friction.
   MaterialID material = createMaterial( "granular", density, e, static_cof, dynamic_cof, real_t( 0.5 ), 1, k, gamma, 0 );

   pe::createPlane( *globalBodyStorage, 0, Vec3(0, 0, 1), forest->getDomain().minCorner(), material );

   SphereID sp = pe::createSphere(
                    *globalBodyStorage,
                    forest->getBlockStorage(),
                    storageID,
                    999999999,
                    Vec3(0,0,1),
                    real_c(1.0),
                    material);
   WALBERLA_CHECK_NOT_NULLPTR(sp);
   sp->setLinearVel( Vec3(0,0,-1) );

   uint_t steps = 0;
   do
   {
      cr.timestep( dt );
      ++steps;
   } while (cr.getNumberOfContacts() != 0);

   WALBERLA_LOG_RESULT(std::setw(30) << "steps: "                      << steps );
   WALBERLA_LOG_RESULT(std::setw(30) << "final velocity: "             << sp->getLinearVel()[2]);
   WALBERLA_LOG_RESULT(std::setw(30) << "final position: "             << sp->getPosition()[2]);
   WALBERLA_LOG_RESULT(std::setw(30) << "integration time: "           << dt);
   WALBERLA_LOG_RESULT(std::setw(30) << "particle radius: "            << radius);
   WALBERLA_LOG_RESULT(std::setw(30) << "particle density: "           << density);
   WALBERLA_LOG_RESULT(std::setw(30) << "particle mass: "              << m);
   WALBERLA_LOG_RESULT(std::setw(30) << "linear spring stiffness: "    << k);
   WALBERLA_LOG_RESULT(std::setw(30) << "damper: "                     << gamma);
   WALBERLA_LOG_RESULT(std::setw(30) << "coefficient of restitution: " << e);
   WALBERLA_LOG_RESULT(std::setw(30) << "collision time: "             << t);

   return EXIT_SUCCESS;
}
}

int main( int argc, char** argv )
{
   return walberla::main(argc, argv);
}
