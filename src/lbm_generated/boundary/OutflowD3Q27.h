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
//! \\file OutflowD3Q27.h
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


class OutflowD3Q27
{
public:
    struct IndexInfo { 
        int32_t x;
        int32_t y;
        int32_t z;
        int32_t dir;
        double pdf;
        double pdf_nd;
        IndexInfo(int32_t x_, int32_t y_, int32_t z_, int32_t dir_) : x(x_), y(y_), z(z_), dir(dir_), pdf(), pdf_nd() {}
        bool operator==(const IndexInfo & o) const {
            return x == o.x && y == o.y && z == o.z && dir == o.dir && floatIsEqual(pdf, o.pdf) && floatIsEqual(pdf_nd, o.pdf_nd);
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

    OutflowD3Q27( const shared_ptr<StructuredBlockForest> & blocks,
                   BlockDataID pdfsID_)
        : pdfsID(pdfsID_)
    {
        auto createIdxVector = []( IBlock * const , StructuredBlockStorage * const ) { return new IndexVectors(); };
        indexVectorID = blocks->addStructuredBlockData< IndexVectors >( createIdxVector, "IndexField_OutflowD3Q27");
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
        auto pdfs = block->getData< field::GhostLayerField<real_t, 27> >(pdfsID); 

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

           if ( isFlagSet( it.neighbor(1, 0, 0 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  4 );
              element.pdf = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(0), it.z() + cell_idx_c(0), 3);
                element.pdf_nd = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(0), it.z() + cell_idx_c(0), 3);
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
              element.pdf = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(1), it.z() + cell_idx_c(0), 9);
                element.pdf_nd = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(1), it.z() + cell_idx_c(0), 9);
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
              element.pdf = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(-1), it.z() + cell_idx_c(0), 7);
                element.pdf_nd = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(-1), it.z() + cell_idx_c(0), 7);
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
              element.pdf = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(0), it.z() + cell_idx_c(1), 17);
                element.pdf_nd = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(0), it.z() + cell_idx_c(1), 17);
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
              element.pdf = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(0), it.z() + cell_idx_c(-1), 13);
                element.pdf_nd = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(0), it.z() + cell_idx_c(-1), 13);
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

           if ( isFlagSet( it.neighbor(1, 1, 1 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  19 );
              element.pdf = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(1), it.z() + cell_idx_c(1), 26);
                element.pdf_nd = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(1), it.z() + cell_idx_c(1), 26);
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

           if ( isFlagSet( it.neighbor(1, -1, 1 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  21 );
              element.pdf = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(-1), it.z() + cell_idx_c(1), 24);
                element.pdf_nd = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(-1), it.z() + cell_idx_c(1), 24);
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

           if ( isFlagSet( it.neighbor(1, 1, -1 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  23 );
              element.pdf = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(1), it.z() + cell_idx_c(-1), 22);
                element.pdf_nd = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(1), it.z() + cell_idx_c(-1), 22);
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

           if ( isFlagSet( it.neighbor(1, -1, -1 , 0 ), boundaryFlag ) )
           {
              auto element = IndexInfo(it.x(), it.y(),  it.z(),  25 );
              element.pdf = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(-1), it.z() + cell_idx_c(-1), 20);
                element.pdf_nd = pdfs->get(it.x() + cell_idx_c(0), it.y() + cell_idx_c(-1), it.z() + cell_idx_c(-1), 20);
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
    BlockDataID pdfsCPUID;
public:
    BlockDataID pdfsID;
};



} // namespace lbm
} // namespace walberla