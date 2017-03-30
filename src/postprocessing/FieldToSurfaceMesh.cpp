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
//! \file FieldToSurfaceMesh.cpp
//! \ingroup postprocessing
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "FieldToSurfaceMesh.h"

#include "core/mpi/MPIWrapper.h"
#include "core/mpi/Gatherv.h"


namespace walberla {
namespace postprocessing {


//*******************************************************************************************************************
/*! Gathers parts of a TriangleMesh on a single process and merges them into a single mesh
*
*   \param mesh           pointer to the mesh that is gathered. The coordinates of the mesh have to be already in
*                         global coordinates.
*   \param gatherNormals  By default only vertices and faces are gathered. If this parameter is true, the resulting
*                         mesh contains also normals
*   \param targetRank     rank of the process where the mesh is gathered
*   \param comm           MPI communicator used for the gather operation
*
*   \return null pointer on all processes with rank not equal targetRank. Pointer to gathered mesh on targetRank.
*
*/
//*******************************************************************************************************************
shared_ptr<TriangleMesh> gatherMesh( const shared_ptr<TriangleMesh> & mesh, bool gatherNormals,
                                     int targetRank, MPI_Comm comm )
{
   WALBERLA_NON_MPI_SECTION() {
      return mesh;
   }

   // Remove normal information, since if !calcNormals also no normal information from other processes
   // is received
   if ( ! gatherNormals ) {
      mesh->getVertexNormals().clear();
      mesh->getNormalIndices().clear();
   }


   mpi::SendBuffer sendBuffer;
   sendBuffer << mesh->getVertices();
   sendBuffer << mesh->getVertexIndices();
   if ( gatherNormals )
   {
      sendBuffer << mesh->getVertexNormals();
      sendBuffer << mesh->getNormalIndices();
   }

   mpi::RecvBuffer recvBuffer;
   mpi::gathervBuffer( sendBuffer, recvBuffer, targetRank, comm );

   if ( recvBuffer.size() > 0 )
   {
      while ( !recvBuffer.isEmpty() )
      {
         geometry::TriangleMesh receivedMesh;
         recvBuffer >> receivedMesh.getVertices();
         recvBuffer >> receivedMesh.getVertexIndices();
         if ( gatherNormals )
         {
            recvBuffer >> receivedMesh.getVertexNormals();
            recvBuffer >> receivedMesh.getNormalIndices();
         }

         mesh->merge( receivedMesh );
      }
      return mesh;
   }
   else
      return shared_ptr<geometry::TriangleMesh>();
}


//*******************************************************************************************************************
/*! Similar to gatherMesh but meshes is reduced into n pieces (instead of 1 piece )
*
*   Useful if meshes get too large to fit into a the memory of a single process.
*   \param pieceOut output parameter, different for each piece (can be used in output filename)
*   \return reduced mesh pieces, or null pointer on other processes
*
*/
//*******************************************************************************************************************
shared_ptr<TriangleMesh> gatherMeshIntoMultiplePieces( int & pieceOut, const shared_ptr<TriangleMesh> & mesh,
                                                       bool gatherNormals, uint_t pieces, MPI_Comm comm )
{
   WALBERLA_NON_MPI_SECTION() {
      return mesh;
   }

   int numProcesses = 0;
   int myRank = 0;
   MPI_Comm_size( comm, &numProcesses );
   MPI_Comm_rank( comm, &myRank );

   const int pieceSize =  numProcesses / int_c( pieces);
   const int myColor   =  myRank / pieceSize;

   pieceOut = myColor;
   MPI_Comm splittedComm;
   MPI_Comm_split( comm, myColor, myRank, &splittedComm );

   auto res = gatherMesh( mesh, gatherNormals, 0, splittedComm );

   MPI_Comm_free( &splittedComm );
   return res;
}



} // namespace postprocessing
} // namespace walberla




