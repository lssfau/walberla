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
//! \file DistributedVTKMeshWriter.h
//! \ingroup mesh
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "VTKMeshWriter.h"

#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/MPITextFile.h"

#include <sstream>
#include <string>

namespace walberla {
namespace mesh {

template< typename MeshType >
class DistributedVTKMeshWriter : public VTKMeshWriter<MeshType> {
public:
  
   DistributedVTKMeshWriter( const shared_ptr<const MeshType> & mesh, const std::string & identifier,
                             const uint_t writeFrequency, const std::string & baseFolder = "vtk_out" );

   void operator()();

};


template< typename MeshType >
DistributedVTKMeshWriter<MeshType>::DistributedVTKMeshWriter( const shared_ptr<const MeshType> & mesh, const std::string & identifier, 
                                                              const uint_t writeFrequency, const std::string & baseFolder )
   : VTKMeshWriter<MeshType>( mesh, identifier, writeFrequency, baseFolder )
{

}


template< typename MeshType >
void DistributedVTKMeshWriter<MeshType>::operator()()
{
   if( this->isWriteScheduled() )
   {
      const int rank         = MPIManager::instance()->rank();
      const int numProcesses = MPIManager::instance()->numProcesses();

      std::ostringstream oss;
      if( rank == 0 )
      {
         this->writePrefix(oss);
      }

      this->writePiece(oss);

      if( rank == numProcesses - 1 )
      {
         this->writePostfix(oss);
      }

      std::ostringstream filePathVtp;
      filePathVtp << this->baseFolder_ << '/' << this->identifier_ << '/' << this->identifier_ << '_' << this->timestep_ << ".vtp";

      mpi::writeMPITextFile( filePathVtp.str(), oss.str() );

      if( rank == 0 )
      {
         std::ostringstream filePathPvd;
         filePathPvd << this->baseFolder_ << '/' << this->identifier_ << ".pvd";
         std::ofstream ofsPvd( filePathPvd.str().c_str() );
         this->writePVD( ofsPvd );
      }
   }

   this->incrementTimeStep();
}


} // namespace mesh
} // namespace walberla
