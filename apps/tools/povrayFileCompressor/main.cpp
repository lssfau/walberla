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
//! \file main.cpp
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include <fstream>
#include <iomanip>
#include <iostream>


#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>


namespace fs = boost::filesystem;

bool verbose;
bool quiet;

#define PRINT_DEF(msg)  if(!quiet)  std::cout << msg;
#define PRINT_VER(msg)  if(verbose) std::cout << msg;
#define PRINT_ERR(msg) {            std::cerr << msg; return EXIT_FAILURE; }
#define PRINT_WAR(msg)              std::cerr << msg;

namespace walberla{
   namespace geometry{

      bool isGreaterMesh( const TriangleMesh& tm0, const TriangleMesh& tm1 ){
         return tm0.getNumTriangles() > tm1.getNumTriangles();
      }
   }
}

int main(int argc, char** argv)
{
   if( argc == 1 || boost::equals(argv[1], "-h" ) ){
      PRINT_DEF( "Usage: PovrayFileCompressor [OPTIONS] <inDir> <outDir>\n"
         << "compresses povray mesh2 files mode\n"
         << "OPTIONS (optional):\n -h help\n -q quiet: no output\n -v verbose: most output\n -n number of vertices per mesh\n -s number of unconnected meshes beginning with biggest mest\n -o maximum distance between two vertices to be merged" )
      return EXIT_SUCCESS;
   }

   if( argc < 3 || argv[argc-1][0] == '-' || argv[argc-2][0] == '-' )
      PRINT_ERR( "Usage: PovrayFileCompressor [-q|-v|-n|-s|-o ${count}] <inDir> <outDir>\n" )

   verbose = false;
   quiet   = false;
   size_t n = 9000u;
   size_t s =    0u;
   walberla::real_t o = walberla::real_t(0);

   if( argc > 3 ) {
      for( walberla::uint_t i = 1; i < walberla::uint_c(argc-2); ++i ) {
         if( boost::equals(argv[i], "-q" ) )
            quiet = true;
         else if( boost::equals(argv[i], "-v" ) )
            verbose = true;
         else if( boost::equals(argv[i], "-n" ) )
            n = walberla::numeric_cast< size_t >( atoi( argv[++i] ) );
         else if( boost::equals(argv[i], "-s" ) )
            s = walberla::numeric_cast< size_t >( atoi( argv[++i] ) );
         else if( boost::equals(argv[i], "-o" ) )
            o = walberla::real_c( atof( argv[++i] ) );
         else if( argv[i][0] != '-' )
            PRINT_ERR( "Usage: PovrayFileCompressor [-q|-v] <inDir> <outDir>\n" )
         else
            PRINT_WAR( "Ignore unknown option " << argv[i] << "\n" )
      }
   }

   if( quiet && verbose )
      PRINT_ERR( "PovrayFileCompressor can't be quiet (-q) and verbose (-v) at the same time\n" )

   fs::path inPath(argv[argc-2]);
   if( !fs::exists(inPath) || !fs::is_directory(inPath) )
      PRINT_ERR( "Path " << inPath << " does not exist or is not a directory!\n" );

   fs::path outPath(argv[argc-1]);
   if( !fs::exists(outPath) )
      fs::create_directories(outPath);

   PRINT_DEF( "Input Path: " << inPath << "\nOutput Path: " << outPath << "\n" )


   std::vector< std::vector<fs::path> > infiles;

   PRINT_DEF( "Collecting files to compress ..." )

   for(auto pit = fs::directory_iterator(inPath); pit != fs::directory_iterator(); ++pit)
   {
      if( !fs::is_directory( pit->status() ) )
         continue;

      PRINT_VER( "\nCollecting files to compress in: " << pit->path())
      std::vector<fs::path> pfiles;
      for(auto tit = fs::directory_iterator(pit->path()); tit != fs::directory_iterator(); ++tit)
      {
         static const boost::regex extensionExpression("\\.dat");
         if( boost::regex_match( tit->path().extension().string(), extensionExpression ) )
            pfiles.push_back(tit->path());
      }
      if( !pfiles.empty() )
         infiles.push_back(pfiles);
   }

   PRINT_VER( "\n" )
   if( infiles.empty() )
      PRINT_ERR( " found no files to compress\n" )
   else
      PRINT_DEF( " found " << infiles.size() * infiles[0].size() << " files.\n" )

   size_t processes = infiles.size();
   size_t timesteps = infiles[0].size();

   if( processes <= 1 )
      PRINT_ERR( "Files from only one process can't be compressed\n" )

   walberla::geometry::TriangleMesh mesh;

   for( size_t t = 0u; t < timesteps; ++t )
   {
      PRINT_VER( "Compress timestep " << t << "\n" )
      PRINT_VER( "Merge Splitted Meshes: \n" )
      PRINT_VER( "Merge Mesh: " << fs::absolute(infiles[0][t]).string().c_str() << "\n" )
      std::ifstream is( fs::absolute(infiles[0][t]).string().c_str() );
      walberla::geometry::readMeshPov(is, mesh, true);
      is.close();

      for( size_t p = 1u; p < processes; ++p ){
         PRINT_VER( "Merge Mesh: " << fs::absolute(infiles[p][t]).string().c_str() << "\n" )
         std::ifstream tis( fs::absolute(infiles[p][t]).string().c_str() );
         walberla::geometry::readMeshPov(tis, mesh, false);
         tis.close();
      }

      PRINT_VER( "Remove Duplicate Vertices ... \n" )
      size_t removed;
      if( o > walberla::real_t(0) )
         removed = mesh.removeDuplicateVertices( o );
      else
         removed = mesh.removeDuplicateVertices( );
      PRINT_VER( "Removed " << removed << " Duplicate Vertices \n" )

      if( s > 0u ){
         PRINT_VER( "Split meshes and use biggest " << s << " meshes \n" )
         std::vector< walberla::geometry::TriangleMesh > meshes;
         mesh.split( meshes );
         std::sort( meshes.begin(), meshes.end(), walberla::geometry::isGreaterMesh );
         mesh.clear();
         for( size_t index = size_t(0u); index < s; ++index ){
            PRINT_VER( "Merge Mesh: " << index << "\n" )
            mesh.merge( meshes[index] );
         }
      }

      fs::path oPath = fs::absolute(outPath) / infiles[0][t].filename();
      PRINT_VER( "Write New Mesh File: " << (oPath.string().c_str()) << "\n" )
      std::ofstream os( oPath.string().c_str() );
      walberla::geometry::writeMeshPov( os, mesh, 0u, n );
      os.close();

      float current = static_cast<float>(t+1) / static_cast<float>(timesteps) * 100.0f;
      if( verbose ) {
         PRINT_VER( " done. (" << std::fixed << std::setprecision(2) << current << "%)\n")
      } else {
         PRINT_DEF( "\r" << std::fixed << std::setprecision(2) << current << "%" << std::flush )
      }
   }

   return 0;
}
