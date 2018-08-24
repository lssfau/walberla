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
//! \file PlainParMetisTest.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \brief Integration test for parmetis library
//
//======================================================================================================================

#include "core/mpi/MPIWrapper.h"
#include "core/load_balancing/ParMetisWrapper.h"
#include "core/debug/CheckFunctions.h"

#include <cstdlib>
#include <iostream>
#include <vector>

int main( int argc, char * argv[] )
{
   #ifdef WALBERLA_BUILD_WITH_MPI
   #ifdef WALBERLA_BUILD_WITH_PARMETIS

   using namespace walberla::core;

   MPI_Init(&argc, &argv);
   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

   const int64_t numVertices  = 25;

   std::vector< int64_t > vtxdist, xadj, adjncy;
   std::vector<double>    xyz;
   vtxdist.push_back( 0 );
   vtxdist.push_back( 5 );
   vtxdist.push_back( 10 );
   vtxdist.push_back( 15 );

   if (world_rank == 0)
   {
      xadj.push_back(  0 );
      xadj.push_back(  2 );
      xadj.push_back(  5 );
      xadj.push_back(  8 );
      xadj.push_back( 11 );
      xadj.push_back( 13 );

      xyz.push_back( 0 );
      xyz.push_back( 2 );
      xyz.push_back( 1 );
      xyz.push_back( 2 );
      xyz.push_back( 2 );
      xyz.push_back( 2 );
      xyz.push_back( 3 );
      xyz.push_back( 2 );
      xyz.push_back( 4 );
      xyz.push_back( 2 );

      adjncy.push_back(  1 );
      adjncy.push_back(  5 );
      adjncy.push_back(  0 );
      adjncy.push_back(  2 );
      adjncy.push_back(  6 );
      adjncy.push_back(  1 );
      adjncy.push_back(  3 );
      adjncy.push_back(  7 );
      adjncy.push_back(  2 );
      adjncy.push_back(  4 );
      adjncy.push_back(  8 );
      adjncy.push_back(  3 );
      adjncy.push_back(  9 );
   }

   if (world_rank == 1)
   {
      xadj.push_back(  0 );
      xadj.push_back(  3 );
      xadj.push_back(  7 );
      xadj.push_back( 11 );
      xadj.push_back( 15 );
      xadj.push_back( 18 );

      xyz.push_back( 0 );
      xyz.push_back( 1 );
      xyz.push_back( 1 );
      xyz.push_back( 1 );
      xyz.push_back( 2 );
      xyz.push_back( 1 );
      xyz.push_back( 3 );
      xyz.push_back( 1 );
      xyz.push_back( 4 );
      xyz.push_back( 1 );

      adjncy.push_back(  0 );
      adjncy.push_back(  6 );
      adjncy.push_back( 10 );
      adjncy.push_back(  1 );
      adjncy.push_back(  5 );
      adjncy.push_back(  7 );
      adjncy.push_back( 11 );
      adjncy.push_back(  2 );
      adjncy.push_back(  6 );
      adjncy.push_back(  8 );
      adjncy.push_back( 12 );
      adjncy.push_back(  3 );
      adjncy.push_back(  7 );
      adjncy.push_back(  9 );
      adjncy.push_back( 13 );
      adjncy.push_back(  4 );
      adjncy.push_back(  8 );
      adjncy.push_back( 14 );
   }

   if (world_rank == 2)
   {
      xadj.push_back(  0 );
      xadj.push_back(  2 );
      xadj.push_back(  5 );
      xadj.push_back(  8 );
      xadj.push_back( 11 );
      xadj.push_back( 13 );

      xyz.push_back( 0 );
      xyz.push_back( 0 );
      xyz.push_back( 1 );
      xyz.push_back( 0 );
      xyz.push_back( 2 );
      xyz.push_back( 0 );
      xyz.push_back( 3 );
      xyz.push_back( 0 );
      xyz.push_back( 4 );
      xyz.push_back( 0 );

      adjncy.push_back(  5 );
      adjncy.push_back( 11 );
      adjncy.push_back(  6 );
      adjncy.push_back( 10 );
      adjncy.push_back( 12 );
      adjncy.push_back(  7 );
      adjncy.push_back( 11 );
      adjncy.push_back( 13 );
      adjncy.push_back(  8 );
      adjncy.push_back( 12 );
      adjncy.push_back( 14 );
      adjncy.push_back(  9 );
      adjncy.push_back( 13 );
   }

   int64_t wgtflag = 0;
   int64_t numflag = 0;
   int64_t ndims   = 2;
   int64_t ncon    = int64_t( 1 );
   int64_t nparts  = int64_t( 5 );
   std::vector< double > tpwgts( static_cast<size_t>(nparts), 1.0 / static_cast<double>( nparts ) );
   std::vector< double > ubvec(  static_cast<size_t>(ncon), 1.05 );
   double  ipc2redist =  1.0;
   int64_t options[] = {0,0,0};
   int64_t edgecut;
   std::vector< int64_t > part( numVertices );
   MPI_Comm comm = MPI_COMM_WORLD;

   std::cout << "ParMETIS_V3_PartKway" << std::endl;
   WALBERLA_CHECK_EQUAL( ParMETIS_V3_PartKway( &(vtxdist.front()), &(xadj.front()), &(adjncy.front()), nullptr, nullptr, &wgtflag, &numflag, &ncon, &nparts,
                                &(tpwgts.front()), &(ubvec.front()), options, &edgecut, &(part.front()), &comm ),
                         METIS_OK );
   std::cout << "ParMETIS_V3_PartGeomKway" << std::endl;
   WALBERLA_CHECK_EQUAL( ParMETIS_V3_PartGeomKway( &(vtxdist.front()), &(xadj.front()), &(adjncy.front()), nullptr, nullptr, &wgtflag, &numflag, &ndims, &(xyz.front()), &ncon, &nparts,
                                &(tpwgts.front()), &(ubvec.front()), options, &edgecut, &(part.front()), &comm ),
                         METIS_OK );
   std::cout << "ParMETIS_V3_PartGeom" << std::endl;
   WALBERLA_CHECK_EQUAL( ParMETIS_V3_PartGeom( &(vtxdist.front()), &ndims, &(xyz.front()), &(part.front()), &comm ),
                         METIS_OK );
   std::cout << "ParMETIS_V3_AdaptiveRepart" << std::endl;
   WALBERLA_CHECK_EQUAL( ParMETIS_V3_AdaptiveRepart( &(vtxdist.front()), &(xadj.front()), &(adjncy.front()), nullptr, nullptr, nullptr, &wgtflag, &numflag, &ncon, &nparts,
                                &(tpwgts.front()), &(ubvec.front()), &ipc2redist, options, &edgecut, &(part.front()), &comm ),
                         METIS_OK );
   std::cout << "ParMETIS_V3_RefineKway" << std::endl;
   WALBERLA_CHECK_EQUAL( ParMETIS_V3_RefineKway( &(vtxdist.front()), &(xadj.front()), &(adjncy.front()), nullptr, nullptr, &wgtflag, &numflag, &ncon, &nparts,
                                &(tpwgts.front()), &(ubvec.front()), options, &edgecut, &(part.front()), &comm ),
                         METIS_OK );

   MPI_Finalize();

   #endif
   #endif

   return 0;
}
