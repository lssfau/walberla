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
//! \file NVTX.h
//! \ingroup cuda
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"

#include <string>

#include <nvToolsExt.h>
#include <nvToolsExtCuda.h>
#include <nvToolsExtCudaRt.h>

namespace walberla{
namespace cuda {

inline void nvtxMarker(const std::string& name, const uint32_t color=0xaaaaaa)
{
    nvtxEventAttributes_t eventAttrib;
    memset(&eventAttrib, 0, NVTX_EVENT_ATTRIB_STRUCT_SIZE);
    eventAttrib.version = NVTX_VERSION;
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
    eventAttrib.colorType = NVTX_COLOR_ARGB;
    eventAttrib.color = 0xFF000000 | color;
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
    eventAttrib.message.ascii = name.c_str();
    nvtxMarkEx(&eventAttrib);
}

inline void nameStream(const cudaStream_t & stream, const std::string & name)
{
    nvtxNameCudaStreamA(stream, name.c_str());
}

class NvtxRange
{
public:
    NvtxRange(const std::string & name, const uint32_t color=0xaaaaaa)
    {
        memset(&eventAttrib, 0, NVTX_EVENT_ATTRIB_STRUCT_SIZE);
        eventAttrib.version = NVTX_VERSION;
        eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
        eventAttrib.colorType = NVTX_COLOR_ARGB;
        eventAttrib.color = 0xFF000000 | color;
        eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
        eventAttrib.message.ascii = name.c_str();
        nvtxRangePushEx(&eventAttrib);
    }
    ~NvtxRange()
    {
        nvtxRangePop();
    }
private:
    nvtxEventAttributes_t eventAttrib;
};


} // namespace cuda
} // namespace walberla