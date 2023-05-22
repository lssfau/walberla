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
//! \\file FreeSlipD3Q19.h
//! \\author pystencils
//======================================================================================================================

#pragma once
#include "core/DataTypes.h"

#include "field/GhostLayerField.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "blockforest/StructuredBlockForest.h"
#include "field/FlagField.h"
#include "core/debug/Debug.h"

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


class FreeSlipD3Q19
{
public:
    struct IndexInfo { 
        int32_t x;
        int32_t y;
        int32_t z;
        int32_t dir;
        int32_t wnx;
        int32_t wny;
        int32_t wnz;
        int32_t ref_dir;
        IndexInfo(int32_t x_, int32_t y_, int32_t z_, int32_t dir_) : x(x_), y(y_), z(z_), dir(dir_), wnx(), wny(), wnz(), ref_dir() {}
        bool operator==(const IndexInfo & o) const {
            return x == o.x && y == o.y && z == o.z && dir == o.dir && wnx == o.wnx && wny == o.wny && wnz == o.wnz && ref_dir == o.ref_dir;
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

        IndexVectors() = default;
        bool operator==(IndexVectors const &other) const { return other.cpuVectors_ == cpuVectors_; }

        CpuIndexVector & indexVector(Type t) { return cpuVectors_[t]; }
        IndexInfo * pointerCpu(Type t)  { return cpuVectors_[t].data(); }

        void syncGPU()
        {
            
        }

    private:
        std::vector<CpuIndexVector> cpuVectors_{NUM_TYPES};

        
    };

    FreeSlipD3Q19( const shared_ptr<StructuredBlockForest> & blocks,
                   BlockDataID pdfsID_)
        : pdfsID(pdfsID_)
    {
        auto createIdxVector = []( IBlock * const , StructuredBlockStorage * const ) { return new IndexVectors(); };
        indexVectorID = blocks->addStructuredBlockData< IndexVectors >( createIdxVector, "IndexField_FreeSlipD3Q19");
    };

    void run (IBlock * block);

    void operator() (IBlock * block)
    {
        run(block);
    }

    void inner (IBlock * block);

    void outer (IBlock * block);

    std::function<void (IBlock *)> getSweep()
    {
        return [this]
               (IBlock * b)
               { this->run(b); };
    }

    std::function<void (IBlock *)> getInnerSweep()
    {
        return [this]
               (IBlock * b)
               { this->inner(b); };
    }

    std::function<void (IBlock *)> getOuterSweep()
    {
        return [this]
               (IBlock * b)
               { this->outer(b); };
    }

    template<typename FlagField_T>
    void fillFromFlagField( const shared_ptr<StructuredBlockForest> & blocks, ConstBlockDataID flagFieldID,
                            FlagUID boundaryFlagUID, FlagUID domainFlagUID)
    {
        for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
            fillFromFlagField<FlagField_T>(&*blockIt, flagFieldID, boundaryFlagUID, domainFlagUID );
    }


    template<typename FlagField_T>
    void fillFromFlagField(IBlock * block, ConstBlockDataID flagFieldID,
                            FlagUID boundaryFlagUID, FlagUID domainFlagUID )
    {
        auto * indexVectors = block->getData< IndexVectors > ( indexVectorID );
        auto & indexVectorAll = indexVectors->indexVector(IndexVectors::ALL);
        auto & indexVectorInner = indexVectors->indexVector(IndexVectors::INNER);
        auto & indexVectorOuter = indexVectors->indexVector(IndexVectors::OUTER);

        auto * flagField = block->getData< FlagField_T > ( flagFieldID );
        

        if( !(flagField->flagExists(boundaryFlagUID) && flagField->flagExists(domainFlagUID) ))
            return;

        auto boundaryFlag = flagField->getFlag(boundaryFlagUID);
        auto domainFlag = flagField->getFlag(domainFlagUID);

        auto inner = flagField->xyzSize();
        inner.expand( cell_idx_t(-1) );

        indexVectorAll.clear();
        indexVectorInner.clear();
        indexVectorOuter.clear();

        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(0, 0, 0 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  0 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(0, 0, 0);
                int32_t ref_dir = 0; // dir: 0
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + 0, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = 0;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + 0, n.z(), domainFlag ) )
                {
                   element.wny = 0;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + 0, domainFlag ) )
                {
                   element.wnz = 0;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = 0;
                   element.wny = 0;
                   element.wnz = 0;
                   ref_dir = 0;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(0, 1, 0 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  1 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(0, 1, 0);
                int32_t ref_dir = 2; // dir: 1
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + 0, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = 0;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + -1, n.z(), domainFlag ) )
                {
                   element.wny = -1;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + 0, domainFlag ) )
                {
                   element.wnz = 0;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = 0;
                   element.wny = -1;
                   element.wnz = 0;
                   ref_dir = 1;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(0, -1, 0 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  2 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(0, -1, 0);
                int32_t ref_dir = 1; // dir: 2
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + 0, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = 0;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + 1, n.z(), domainFlag ) )
                {
                   element.wny = 1;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + 0, domainFlag ) )
                {
                   element.wnz = 0;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = 0;
                   element.wny = 1;
                   element.wnz = 0;
                   ref_dir = 2;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(-1, 0, 0 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  3 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(-1, 0, 0);
                int32_t ref_dir = 4; // dir: 3
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + 1, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = 1;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + 0, n.z(), domainFlag ) )
                {
                   element.wny = 0;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + 0, domainFlag ) )
                {
                   element.wnz = 0;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = 1;
                   element.wny = 0;
                   element.wnz = 0;
                   ref_dir = 3;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(1, 0, 0 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  4 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(1, 0, 0);
                int32_t ref_dir = 3; // dir: 4
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + -1, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = -1;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + 0, n.z(), domainFlag ) )
                {
                   element.wny = 0;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + 0, domainFlag ) )
                {
                   element.wnz = 0;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = -1;
                   element.wny = 0;
                   element.wnz = 0;
                   ref_dir = 4;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(0, 0, 1 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  5 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(0, 0, 1);
                int32_t ref_dir = 6; // dir: 5
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + 0, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = 0;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + 0, n.z(), domainFlag ) )
                {
                   element.wny = 0;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + -1, domainFlag ) )
                {
                   element.wnz = -1;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = 0;
                   element.wny = 0;
                   element.wnz = -1;
                   ref_dir = 5;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(0, 0, -1 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  6 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(0, 0, -1);
                int32_t ref_dir = 5; // dir: 6
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + 0, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = 0;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + 0, n.z(), domainFlag ) )
                {
                   element.wny = 0;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + 1, domainFlag ) )
                {
                   element.wnz = 1;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = 0;
                   element.wny = 0;
                   element.wnz = 1;
                   ref_dir = 6;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(-1, 1, 0 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  7 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(-1, 1, 0);
                int32_t ref_dir = 10; // dir: 7
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + 1, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = 1;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + -1, n.z(), domainFlag ) )
                {
                   element.wny = -1;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + 0, domainFlag ) )
                {
                   element.wnz = 0;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = 1;
                   element.wny = -1;
                   element.wnz = 0;
                   ref_dir = 7;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(1, 1, 0 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  8 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(1, 1, 0);
                int32_t ref_dir = 9; // dir: 8
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + -1, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = -1;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + -1, n.z(), domainFlag ) )
                {
                   element.wny = -1;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + 0, domainFlag ) )
                {
                   element.wnz = 0;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = -1;
                   element.wny = -1;
                   element.wnz = 0;
                   ref_dir = 8;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(-1, -1, 0 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  9 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(-1, -1, 0);
                int32_t ref_dir = 8; // dir: 9
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + 1, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = 1;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + 1, n.z(), domainFlag ) )
                {
                   element.wny = 1;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + 0, domainFlag ) )
                {
                   element.wnz = 0;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = 1;
                   element.wny = 1;
                   element.wnz = 0;
                   ref_dir = 9;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(1, -1, 0 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  10 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(1, -1, 0);
                int32_t ref_dir = 7; // dir: 10
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + -1, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = -1;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + 1, n.z(), domainFlag ) )
                {
                   element.wny = 1;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + 0, domainFlag ) )
                {
                   element.wnz = 0;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = -1;
                   element.wny = 1;
                   element.wnz = 0;
                   ref_dir = 10;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(0, 1, 1 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  11 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(0, 1, 1);
                int32_t ref_dir = 16; // dir: 11
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + 0, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = 0;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + -1, n.z(), domainFlag ) )
                {
                   element.wny = -1;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + -1, domainFlag ) )
                {
                   element.wnz = -1;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = 0;
                   element.wny = -1;
                   element.wnz = -1;
                   ref_dir = 11;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(0, -1, 1 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  12 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(0, -1, 1);
                int32_t ref_dir = 15; // dir: 12
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + 0, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = 0;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + 1, n.z(), domainFlag ) )
                {
                   element.wny = 1;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + -1, domainFlag ) )
                {
                   element.wnz = -1;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = 0;
                   element.wny = 1;
                   element.wnz = -1;
                   ref_dir = 12;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(-1, 0, 1 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  13 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(-1, 0, 1);
                int32_t ref_dir = 18; // dir: 13
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + 1, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = 1;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + 0, n.z(), domainFlag ) )
                {
                   element.wny = 0;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + -1, domainFlag ) )
                {
                   element.wnz = -1;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = 1;
                   element.wny = 0;
                   element.wnz = -1;
                   ref_dir = 13;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(1, 0, 1 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  14 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(1, 0, 1);
                int32_t ref_dir = 17; // dir: 14
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + -1, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = -1;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + 0, n.z(), domainFlag ) )
                {
                   element.wny = 0;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + -1, domainFlag ) )
                {
                   element.wnz = -1;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = -1;
                   element.wny = 0;
                   element.wnz = -1;
                   ref_dir = 14;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(0, 1, -1 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  15 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(0, 1, -1);
                int32_t ref_dir = 12; // dir: 15
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + 0, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = 0;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + -1, n.z(), domainFlag ) )
                {
                   element.wny = -1;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + 1, domainFlag ) )
                {
                   element.wnz = 1;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = 0;
                   element.wny = -1;
                   element.wnz = 1;
                   ref_dir = 15;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(0, -1, -1 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  16 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(0, -1, -1);
                int32_t ref_dir = 11; // dir: 16
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + 0, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = 0;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + 1, n.z(), domainFlag ) )
                {
                   element.wny = 1;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + 1, domainFlag ) )
                {
                   element.wnz = 1;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = 0;
                   element.wny = 1;
                   element.wnz = 1;
                   ref_dir = 16;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(-1, 0, -1 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  17 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(-1, 0, -1);
                int32_t ref_dir = 14; // dir: 17
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + 1, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = 1;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + 0, n.z(), domainFlag ) )
                {
                   element.wny = 0;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + 1, domainFlag ) )
                {
                   element.wnz = 1;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = 1;
                   element.wny = 0;
                   element.wnz = 1;
                   ref_dir = 17;
                }
                element.ref_dir = ref_dir;
              indexVectorAll.push_back( element );
              if( inner.contains( it.x(), it.y(), it.z() ) )
                 indexVectorInner.push_back( element );
              else
                 indexVectorOuter.push_back( element );
           }
        }
        
        for( auto it = flagField->beginWithGhostLayerXYZ( cell_idx_c( flagField->nrOfGhostLayers() - 1 ) ); it != flagField->end(); ++it )
        {
           if( ! isFlagSet(it, domainFlag) )
              continue;

           if ( isFlagSet( it.neighbor(1, 0, -1 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  18 );
              const int32_t x_axis_mirrored_stencil_dir [] = { 0,1,2,4,3,5,6,8,7,10,9,11,12,14,13,15,16,18,17 };
                const int32_t y_axis_mirrored_stencil_dir [] = { 0,2,1,3,4,5,6,9,10,7,8,12,11,13,14,16,15,17,18 };
                const int32_t z_axis_mirrored_stencil_dir [] = { 0,1,2,3,4,6,5,7,8,9,10,15,16,17,18,11,12,13,14 };
                const Cell n = it.cell() + Cell(1, 0, -1);
                int32_t ref_dir = 13; // dir: 18
                element.wnx = 0; // compute discrete normal vector of free slip wall
                element.wny = 0;
                if( flagField->isPartOfMaskSet( n.x() + -1, n.y(), n.z(), domainFlag ) )
                {
                   element.wnx = -1;
                   ref_dir = x_axis_mirrored_stencil_dir[ ref_dir ];
                }
                if( flagField->isPartOfMaskSet( n.x(), n.y() + 0, n.z(), domainFlag ) )
                {
                   element.wny = 0;
                   ref_dir = y_axis_mirrored_stencil_dir[ ref_dir ];
                }
                element.wnz = 0;
                if( flagField->isPartOfMaskSet( n.x(), n.y(), n.z() + 1, domainFlag ) )
                {
                   element.wnz = 1;
                   ref_dir = z_axis_mirrored_stencil_dir[ ref_dir ];
                }
                // concave corner (neighbors are non-fluid)
                if( element.wnx == 0 && element.wny == 0 && element.wnz == 0 )
                {
                   element.wnx = -1;
                   element.wny = 0;
                   element.wnz = 1;
                   ref_dir = 18;
                }
                element.ref_dir = ref_dir;
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
    void run_impl(IBlock * block, IndexVectors::Type type);

    BlockDataID indexVectorID;
    
public:
    BlockDataID pdfsID;
};



} // namespace lbm
} // namespace walberla