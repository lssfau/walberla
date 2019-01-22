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
//! \\file UniformGridGPU_UBB.h
//! \\author pystencils
//======================================================================================================================


#include "core/DataTypes.h"

#include "cuda/GPUField.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "blockforest/StructuredBlockForest.h"
#include "field/FlagField.h"

#include <set>
#include <vector>

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

namespace walberla {
namespace lbm {


class UniformGridGPU_UBB
{
public:
    struct IndexInfo { 
        int32_t x;
        int32_t y;
        int32_t z;
        int32_t dir;
        IndexInfo(int32_t x_, int32_t y_, int32_t z_, int32_t dir_) : x(x_), y(y_), z(z_), dir(dir_) {}
        bool operator==(const IndexInfo & o) const {
            return x == o.x && y == o.y && z == o.z && dir == o.dir;
        }
    };



    class IndexVectors
    {
    public:
        using CpuIndexVector = std::vector<IndexInfo>;

        enum Type {
            ALL = 0,
            INNER = 1,
            OUTER = 2,
            NUM_TYPES = 3
        };

        IndexVectors() : cpuVectors_(NUM_TYPES)  {}
        bool operator==(IndexVectors & other) { return other.cpuVectors_ == cpuVectors_; }

        ~IndexVectors() {
            for( auto & gpuVec: gpuVectors_)
                cudaFree( gpuVec );
        }
        

        CpuIndexVector & indexVector(Type t) { return cpuVectors_[t]; }
        IndexInfo * pointerCpu(Type t)  { return &(cpuVectors_[t][0]); }

        IndexInfo * pointerGpu(Type t)  { return gpuVectors_[t]; }
        

        void syncGPU()
        {
            gpuVectors_.resize( cpuVectors_.size() );
            for(size_t i=0; i < size_t(NUM_TYPES); ++i )
            {
                auto & gpuVec = gpuVectors_[i];
                auto & cpuVec = cpuVectors_[i];
                cudaFree( gpuVec );
                cudaMalloc( &gpuVec, sizeof(IndexInfo) * cpuVec.size() );
                cudaMemcpy( gpuVec, &cpuVec[0], sizeof(IndexInfo) * cpuVec.size(), cudaMemcpyHostToDevice );
            }
        }

    private:
        std::vector<CpuIndexVector> cpuVectors_;

        using GpuIndexVector = IndexInfo *;
        std::vector<GpuIndexVector> gpuVectors_;
        
    };


    UniformGridGPU_UBB( const shared_ptr<StructuredBlockForest> & blocks,
                   BlockDataID pdfsID_ )
        : pdfsID(pdfsID_)
    {
        auto createIdxVector = []( IBlock * const , StructuredBlockStorage * const ) { return new IndexVectors(); };
        indexVectorID = blocks->addStructuredBlockData< IndexVectors >( createIdxVector, "IndexField_UniformGridGPU_UBB");
    };

    void operator() ( IBlock * block , cudaStream_t stream = 0 );
    void inner( IBlock * block , cudaStream_t stream = 0 );
    void outer( IBlock * block , cudaStream_t stream = 0 );


    template<typename FlagField_T>
    void fillFromFlagField( const shared_ptr<StructuredBlockForest> & blocks, ConstBlockDataID flagFieldID,
                            FlagUID boundaryFlagUID, FlagUID domainFlagUID)
    {
        for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
            fillFromFlagField<FlagField_T>( &*blockIt, flagFieldID, boundaryFlagUID, domainFlagUID );
    }


    template<typename FlagField_T>
    void fillFromFlagField( IBlock * block, ConstBlockDataID flagFieldID,
                            FlagUID boundaryFlagUID, FlagUID domainFlagUID )
    {
        auto * indexVectors = block->getData< IndexVectors > ( indexVectorID );
        auto & indexVectorAll = indexVectors->indexVector(IndexVectors::ALL);
        auto & indexVectorInner = indexVectors->indexVector(IndexVectors::INNER);
        auto & indexVectorOuter = indexVectors->indexVector(IndexVectors::OUTER);


        auto * flagField = block->getData< FlagField_T > ( flagFieldID );

        auto boundaryFlag = flagField->getFlag(boundaryFlagUID);
        auto domainFlag = flagField->getFlag(domainFlagUID);

        auto inner = flagField->xyzSize();
        inner.expand( cell_idx_t(-1) );


        indexVectorAll.clear();
        indexVectorInner.clear();
        indexVectorOuter.clear();

        for( auto it = flagField->begin(); it != flagField->end(); ++it )
        {
            if( ! isFlagSet(it, domainFlag) )
                continue;
            if ( isFlagSet( it.neighbor(0, 0, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  0 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, 1, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  1 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, -1, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  2 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(-1, 0, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  3 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(1, 0, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  4 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, 0, 1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  5 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, 0, -1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  6 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(-1, 1, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  7 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(1, 1, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  8 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(-1, -1, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  9 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(1, -1, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  10 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, 1, 1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  11 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, -1, 1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  12 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(-1, 0, 1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  13 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(1, 0, 1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  14 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, 1, -1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  15 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, -1, -1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  16 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(-1, 0, -1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  17 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(1, 0, -1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  18 );
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
        }

        indexVectors->syncGPU();
    }

private:
    void run( IBlock * block, IndexVectors::Type type, cudaStream_t stream = 0 );

    BlockDataID indexVectorID;

    BlockDataID pdfsID;
};



} // namespace lbm
} // namespace walberla