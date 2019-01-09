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
//! \\file UniformGridGPU_LbKernel.h
//! \\author pystencils
//======================================================================================================================

#include "core/DataTypes.h"

#include "cuda/GPUField.h"
#include "cuda/ParallelStreams.h"
#include "field/SwapableCompare.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"

#include <set>

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

namespace walberla {
namespace pystencils {


class UniformGridGPU_LbKernel
{
public:
    UniformGridGPU_LbKernel( BlockDataID pdfsID_, double omega_)
        : pdfsID(pdfsID_), omega(omega_)
    {};

    
    ~UniformGridGPU_LbKernel() {  
        for(auto p: cache_pdfs_) {
            delete p;
        }
     }



    void operator() ( IBlock * block , cudaStream_t stream = 0 );

    void inner( IBlock * block , cudaStream_t stream = 0 );
    void outer( IBlock * block , cudaStream_t stream = 0 );

    void setOuterPriority(int priority ) {
        
        parallelStreams_.setStreamPriority(priority);
        
    }
private:
    BlockDataID pdfsID;
    double omega;

    std::set< cuda::GPUField<double> *, field::SwapableCompare< cuda::GPUField<double> * > > cache_pdfs_;

    
    cuda::ParallelStreams parallelStreams_;
    
};


} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif