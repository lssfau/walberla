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
//! \file ErrorChecking.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Abort.h"

#include <sstream>

#include "gpu/DeviceWrapper.h"

namespace walberla {
namespace gpu {


#define WALBERLA_GPU_CHECK(ans) { ::walberla::gpu::checkForError((ans), __FILE__, __LINE__); }
#define WALBERLA_GPU_CHECK_LAST_ERROR() { ::walberla::gpu::checkForLastError(__FILE__, __LINE__); }



inline void checkForError( gpuError_t code, const std::string & callerPath, const int line )
{
   // Oftentimes CUDA functions return an error code (if error has occurred) This function converts the error string in human-readable output.
   // For general error checking use checkForLastError
  if(code != gpuSuccess)
  {
    std::stringstream ss;
    ss << "GPU Error: " << code << " " << gpuGetErrorName(code) << ": " << gpuGetErrorString( code );
    Abort::instance()->abort( ss.str(), callerPath, line );
  }
}

#ifndef NDEBUG
inline void checkForLastError( const std::string & callerPath, const int line )
{
   // Forces immediate checking with a synchronizing. This breaks asynchrony/concurrency structure. Thus, only in debug mode executed.
   gpuError_t code = gpuGetLastError();
   if(code != gpuSuccess)
   {
      std::stringstream ss;
      ss << "CUDA Error: " << code << " " << gpuGetErrorName(code) << ": " << gpuGetErrorString( code );
      Abort::instance()->abort( ss.str(), callerPath, line );
   }
}
#else
inline void checkForLastError( const std::string & /*callerPath*/, const int /*line*/ ){}
#endif


} // namespace gpu
} // namespace walberla


