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
//! \file FieldToSurfaceMesh.h
//! \ingroup postprocessing
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "MarchingCubes.h"
#include "geometry/mesh/TriangleMesh.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/FlagField.h"

#include <type_traits>

namespace walberla {
namespace postprocessing {



   shared_ptr<TriangleMesh> gatherMesh( const shared_ptr<TriangleMesh> & mesh, bool gatherNormals,
                                        int targetRank=0, MPI_Comm comm = MPI_COMM_WORLD  );


   shared_ptr<TriangleMesh> gatherMeshIntoMultiplePieces( int & pieceOut, const shared_ptr<TriangleMesh> & mesh, bool gatherNormals,
                                                          uint_t pieces=1, MPI_Comm comm = MPI_COMM_WORLD );





   //*******************************************************************************************************************
   /*! Generates isosurface mesh out of a flag field and gathers the mesh on a single process.
   *
   * All cells where on of the bits of the mask are set, are located inside the resulting mesh.
   *
   * Warning: Slow due to sub-optimal implementation - uses a temporary real field.
   */
   //*******************************************************************************************************************
   template<typename Field_T>
   shared_ptr<geometry::TriangleMesh> flagFieldToSurfaceMesh( const shared_ptr<StructuredBlockStorage> & bs,
                                                              ConstBlockDataID fieldID, typename Field_T::value_type mask,
                                                              bool calcNormals = false,
                                                              int targetRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
   {
      static_assert( std::is_unsigned<typename Field_T::value_type>::value, "Works only for FlagFields" );
      static_assert( Field_T::F_SIZE == 1, "Works only for FlagFields" );


      // Do not use ghost layer information in non-perdiodic directions
      CellInterval domainCellBB = bs->getDomainCellBB();
      for( size_t i=0; i<3; ++i )
         if ( ! bs->isPeriodic( i) )
         {
            domainCellBB.min()[i] ++;
            domainCellBB.max()[i] --;
         }


      auto mesh = make_shared<geometry::TriangleMesh> ();
      for( auto block = bs->begin(); block != bs->end(); ++block )
      {
         const Field_T * field = block->template getData<Field_T>( fieldID );

         CellInterval localDomainCellBB;
         bs->transformGlobalToBlockLocalCellInterval( localDomainCellBB, *block, domainCellBB );
         localDomainCellBB.intersect( field->xyzSize() );

         CellInterval blockCi = bs->getBlockCellBB( *block );
         Vector3<real_t> offset( real_c( blockCi.min()[0]) * bs->dx(),
                                 real_c( blockCi.min()[1]) * bs->dy(),
                                 real_c( blockCi.min()[2]) * bs->dz() );

         GhostLayerField<real_t,1> tmpField ( field->xSize(), field->ySize(), field->zSize(), field->nrOfGhostLayers() );
         auto realFieldIt = tmpField.beginWithGhostLayerXYZ();
         auto flagFieldIt = field->beginWithGhostLayerXYZ();
         while ( realFieldIt != tmpField.end() )
         {
            if( isPartOfMaskSet( flagFieldIt, mask) )
               *realFieldIt = real_t( 1.0 );
            else
               *realFieldIt = real_t( 0.0 );

            ++realFieldIt;
            ++flagFieldIt;
         }

         generateIsoSurface( tmpField, 0.5, *mesh, Vector3<real_t>( bs->dx(), bs->dy(), bs->dz() ),
                             0, offset, localDomainCellBB, calcNormals );
      }

      return gatherMesh( mesh, calcNormals, targetRank, comm );
   }



   //*******************************************************************************************************************
   /*! Generates isosurface mesh out of a real valued field and gathers the mesh on a single process.
   *
   */
   //*******************************************************************************************************************
   template<typename Field_T>
   shared_ptr<geometry::TriangleMesh> realFieldToSurfaceMesh( const shared_ptr<StructuredBlockStorage> & bs,
                                                              ConstBlockDataID fieldID, real_t threshold,
                                                              uint_t fCoord = 0, bool calcNormals = false,
                                                              int targetRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
   {
      static_assert( (std::is_same< typename Field_T::value_type, real_t >::value), "Function works only for fields of real" );

      auto mesh = make_shared<geometry::TriangleMesh> ();


      // Do not use ghost layer information in non-perdiodic directions
      CellInterval domainCellBB = bs->getDomainCellBB();
      for( size_t i=0; i<3; ++i )
         if ( ! bs->isPeriodic( i) )
         {
            domainCellBB.min()[i] ++;
            domainCellBB.max()[i] --;
         }


      for( auto block = bs->begin(); block != bs->end(); ++block )
      {
         const Field_T * field = block->template getData<Field_T>( fieldID );

         CellInterval localDomainCellBB;
         bs->transformGlobalToBlockLocalCellInterval( localDomainCellBB, *block, domainCellBB );
         localDomainCellBB.intersect( field->xyzSize() );

         CellInterval blockCi = bs->getBlockCellBB( *block );
         Vector3<real_t> offset( real_c( blockCi.min()[0]) * bs->dx(),
                                 real_c( blockCi.min()[1]) * bs->dy(),
                                 real_c( blockCi.min()[2]) * bs->dz() );

         generateIsoSurface( *field, threshold, *mesh, Vector3<real_t>( bs->dx(), bs->dy(), bs->dz() ),
                             fCoord, offset, localDomainCellBB, calcNormals );
      }

      return gatherMesh( mesh, calcNormals, targetRank, comm );
   }




} // namespace postprocessing
} // namespace walberla


