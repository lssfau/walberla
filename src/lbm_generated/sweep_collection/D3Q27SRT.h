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
//! \\file D3Q27SRT.h
//! \\author pystencils
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/logging/Logging.h"
#include "core/Macros.h"



#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/SwapableCompare.h"
#include "field/GhostLayerField.h"

#include <set>
#include <cmath>



using namespace std::placeholders;

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
#   pragma GCC diagnostic ignored "-Wreorder"
#endif

namespace walberla {
namespace lbm {


class D3Q27SRT
{
public:
  enum Type { ALL = 0, INNER = 1, OUTER = 2 };

   D3Q27SRT(const shared_ptr< StructuredBlockStorage > & blocks, BlockDataID pdfsID_, BlockDataID densityID_, BlockDataID velocityID_, double omega, const Cell & outerWidth=Cell(1, 1, 1))
     : blocks_(blocks), pdfsID(pdfsID_), densityID(densityID_), velocityID(velocityID_), omega_(omega), outerWidth_(outerWidth)
   {
      

      for (auto& iBlock : *blocks)
      {
         if (int_c(blocks->getNumberOfXCells(iBlock)) <= outerWidth_[0] * 2 ||
             int_c(blocks->getNumberOfYCells(iBlock)) <= outerWidth_[1] * 2 ||
             int_c(blocks->getNumberOfZCells(iBlock)) <= outerWidth_[2] * 2)
          WALBERLA_ABORT_NO_DEBUG_INFO("innerOuterSplit too large - make it smaller or increase cellsPerBlock")
      }
   };

   
    ~D3Q27SRT() {  
        for(auto p: cache_pdfs_) {
            delete p;
        }
     }


   /*************************************************************************************
   *                Internal Function Definitions with raw Pointer
   *************************************************************************************/
   static void streamCollide (field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 27> * pdfs_tmp, double omega, const cell_idx_t ghost_layers = 0);
   static void streamCollideCellInterval (field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 27> * pdfs_tmp, double omega, const CellInterval & ci);
   
   static void collide (field::GhostLayerField<double, 27> * pdfs, double omega, const cell_idx_t ghost_layers = 0);
   static void collideCellInterval (field::GhostLayerField<double, 27> * pdfs, double omega, const CellInterval & ci);
   
   static void stream (field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 27> * pdfs_tmp, const cell_idx_t ghost_layers = 0);
   static void streamCellInterval (field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 27> * pdfs_tmp, const CellInterval & ci);
   
   static void streamOnlyNoAdvancement (field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 27> * pdfs_tmp, const cell_idx_t ghost_layers = 0);
   static void streamOnlyNoAdvancementCellInterval (field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 27> * pdfs_tmp, const CellInterval & ci);
   
   static void initialise (field::GhostLayerField<double, 1> * density, field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 3> * velocity, const cell_idx_t ghost_layers = 0);
   static void initialiseCellInterval (field::GhostLayerField<double, 1> * density, field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 3> * velocity, const CellInterval & ci);
   
   static void calculateMacroscopicParameters (field::GhostLayerField<double, 1> * density, field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 3> * velocity, const cell_idx_t ghost_layers = 0);
   static void calculateMacroscopicParametersCellInterval (field::GhostLayerField<double, 1> * density, field::GhostLayerField<double, 27> * pdfs, field::GhostLayerField<double, 3> * velocity, const CellInterval & ci);
   

   /*************************************************************************************
   *                Function Definitions for external Usage
   *************************************************************************************/

   std::function<void (IBlock *)> streamCollide()
   {
      return [this](IBlock* block) { streamCollide(block); };
   }

   std::function<void (IBlock *)> streamCollide(Type type)
   {
      switch (type)
      {
         case Type::INNER:
            return [this](IBlock* block) { streamCollideInner(block); };
         case Type::OUTER:
            return [this](IBlock* block) { streamCollideOuter(block); };
         default:
            return [this](IBlock* block) { streamCollide(block); };
      }
   }

   std::function<void (IBlock *)> streamCollide(Type type, const cell_idx_t ghost_layers)
   {
      switch (type)
      {
         case Type::INNER:
            return [this](IBlock* block) { streamCollideInner(block); };
         case Type::OUTER:
            return [this](IBlock* block) { streamCollideOuter(block); };
         default:
            return [this, ghost_layers](IBlock* block) { streamCollide(block, ghost_layers); };
      }
   }

   

   void streamCollide(IBlock * block)
   {
      const cell_idx_t ghost_layers = 0;
      

      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      auto & omega = this->omega_;
      
      streamCollide(pdfs, pdfs_tmp, omega, ghost_layers);
      pdfs->swapDataPointers(pdfs_tmp);

   }

   void streamCollide(IBlock * block, const cell_idx_t ghost_layers)
   {
      

      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      auto & omega = this->omega_;
      
      streamCollide(pdfs, pdfs_tmp, omega, ghost_layers);
      pdfs->swapDataPointers(pdfs_tmp);

   }

   

   void streamCollideCellInterval(IBlock * block, const CellInterval & ci)
   {
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      auto & omega = this->omega_;
      
      streamCollideCellInterval(pdfs, pdfs_tmp, omega, ci);
      pdfs->swapDataPointers(pdfs_tmp);

   }

   void streamCollideInner(IBlock * block)
   {
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      auto & omega = this->omega_;
      

      CellInterval inner = pdfs->xyzSize();
      inner.expand(Cell(-outerWidth_[0], -outerWidth_[1], -outerWidth_[2]));

      streamCollideCellInterval(pdfs, pdfs_tmp, omega, inner);
   }

   void streamCollideOuter(IBlock * block)
   {

      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      auto & omega = this->omega_;
      

      if( layers_.empty() )
      {
         CellInterval ci;

         pdfs->getSliceBeforeGhostLayer(stencil::T, ci, outerWidth_[2], false);
         layers_.push_back(ci);
         pdfs->getSliceBeforeGhostLayer(stencil::B, ci, outerWidth_[2], false);
         layers_.push_back(ci);

         pdfs->getSliceBeforeGhostLayer(stencil::N, ci, outerWidth_[1], false);
         ci.expand(Cell(0, 0, -outerWidth_[2]));
         layers_.push_back(ci);
         pdfs->getSliceBeforeGhostLayer(stencil::S, ci, outerWidth_[1], false);
         ci.expand(Cell(0, 0, -outerWidth_[2]));
         layers_.push_back(ci);

         pdfs->getSliceBeforeGhostLayer(stencil::E, ci, outerWidth_[0], false);
         ci.expand(Cell(0, -outerWidth_[1], -outerWidth_[2]));
         layers_.push_back(ci);
         pdfs->getSliceBeforeGhostLayer(stencil::W, ci, outerWidth_[0], false);
         ci.expand(Cell(0, -outerWidth_[1], -outerWidth_[2]));
         layers_.push_back(ci);
      }

    
      for( auto & ci: layers_ )
      {
         streamCollideCellInterval(pdfs, pdfs_tmp, omega, ci);
      }
    

    pdfs->swapDataPointers(pdfs_tmp);

   }
   

   std::function<void (IBlock *)> collide()
   {
      return [this](IBlock* block) { collide(block); };
   }

   std::function<void (IBlock *)> collide(Type type)
   {
      switch (type)
      {
         case Type::INNER:
            return [this](IBlock* block) { collideInner(block); };
         case Type::OUTER:
            return [this](IBlock* block) { collideOuter(block); };
         default:
            return [this](IBlock* block) { collide(block); };
      }
   }

   std::function<void (IBlock *)> collide(Type type, const cell_idx_t ghost_layers)
   {
      switch (type)
      {
         case Type::INNER:
            return [this](IBlock* block) { collideInner(block); };
         case Type::OUTER:
            return [this](IBlock* block) { collideOuter(block); };
         default:
            return [this, ghost_layers](IBlock* block) { collide(block, ghost_layers); };
      }
   }

   

   void collide(IBlock * block)
   {
      const cell_idx_t ghost_layers = 0;
      

      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);

      auto & omega = this->omega_;
      
      collide(pdfs, omega, ghost_layers);
      
   }

   void collide(IBlock * block, const cell_idx_t ghost_layers)
   {
      

      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);

      auto & omega = this->omega_;
      
      collide(pdfs, omega, ghost_layers);
      
   }

   

   void collideCellInterval(IBlock * block, const CellInterval & ci)
   {
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);

      auto & omega = this->omega_;
      
      collideCellInterval(pdfs, omega, ci);
      
   }

   void collideInner(IBlock * block)
   {
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);

      auto & omega = this->omega_;
      

      CellInterval inner = pdfs->xyzSize();
      inner.expand(Cell(-outerWidth_[0], -outerWidth_[1], -outerWidth_[2]));

      collideCellInterval(pdfs, omega, inner);
   }

   void collideOuter(IBlock * block)
   {

      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);

      auto & omega = this->omega_;
      

      if( layers_.empty() )
      {
         CellInterval ci;

         pdfs->getSliceBeforeGhostLayer(stencil::T, ci, outerWidth_[2], false);
         layers_.push_back(ci);
         pdfs->getSliceBeforeGhostLayer(stencil::B, ci, outerWidth_[2], false);
         layers_.push_back(ci);

         pdfs->getSliceBeforeGhostLayer(stencil::N, ci, outerWidth_[1], false);
         ci.expand(Cell(0, 0, -outerWidth_[2]));
         layers_.push_back(ci);
         pdfs->getSliceBeforeGhostLayer(stencil::S, ci, outerWidth_[1], false);
         ci.expand(Cell(0, 0, -outerWidth_[2]));
         layers_.push_back(ci);

         pdfs->getSliceBeforeGhostLayer(stencil::E, ci, outerWidth_[0], false);
         ci.expand(Cell(0, -outerWidth_[1], -outerWidth_[2]));
         layers_.push_back(ci);
         pdfs->getSliceBeforeGhostLayer(stencil::W, ci, outerWidth_[0], false);
         ci.expand(Cell(0, -outerWidth_[1], -outerWidth_[2]));
         layers_.push_back(ci);
      }

    
      for( auto & ci: layers_ )
      {
         collideCellInterval(pdfs, omega, ci);
      }
    

    
   }
   

   std::function<void (IBlock *)> stream()
   {
      return [this](IBlock* block) { stream(block); };
   }

   std::function<void (IBlock *)> stream(Type type)
   {
      switch (type)
      {
         case Type::INNER:
            return [this](IBlock* block) { streamInner(block); };
         case Type::OUTER:
            return [this](IBlock* block) { streamOuter(block); };
         default:
            return [this](IBlock* block) { stream(block); };
      }
   }

   std::function<void (IBlock *)> stream(Type type, const cell_idx_t ghost_layers)
   {
      switch (type)
      {
         case Type::INNER:
            return [this](IBlock* block) { streamInner(block); };
         case Type::OUTER:
            return [this](IBlock* block) { streamOuter(block); };
         default:
            return [this, ghost_layers](IBlock* block) { stream(block, ghost_layers); };
      }
   }

   

   void stream(IBlock * block)
   {
      const cell_idx_t ghost_layers = 0;
      

      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      
      
      stream(pdfs, pdfs_tmp, ghost_layers);
      pdfs->swapDataPointers(pdfs_tmp);

   }

   void stream(IBlock * block, const cell_idx_t ghost_layers)
   {
      

      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      
      
      stream(pdfs, pdfs_tmp, ghost_layers);
      pdfs->swapDataPointers(pdfs_tmp);

   }

   

   void streamCellInterval(IBlock * block, const CellInterval & ci)
   {
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      
      
      streamCellInterval(pdfs, pdfs_tmp, ci);
      pdfs->swapDataPointers(pdfs_tmp);

   }

   void streamInner(IBlock * block)
   {
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      
      

      CellInterval inner = pdfs->xyzSize();
      inner.expand(Cell(-outerWidth_[0], -outerWidth_[1], -outerWidth_[2]));

      streamCellInterval(pdfs, pdfs_tmp, inner);
   }

   void streamOuter(IBlock * block)
   {

      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      
      

      if( layers_.empty() )
      {
         CellInterval ci;

         pdfs->getSliceBeforeGhostLayer(stencil::T, ci, outerWidth_[2], false);
         layers_.push_back(ci);
         pdfs->getSliceBeforeGhostLayer(stencil::B, ci, outerWidth_[2], false);
         layers_.push_back(ci);

         pdfs->getSliceBeforeGhostLayer(stencil::N, ci, outerWidth_[1], false);
         ci.expand(Cell(0, 0, -outerWidth_[2]));
         layers_.push_back(ci);
         pdfs->getSliceBeforeGhostLayer(stencil::S, ci, outerWidth_[1], false);
         ci.expand(Cell(0, 0, -outerWidth_[2]));
         layers_.push_back(ci);

         pdfs->getSliceBeforeGhostLayer(stencil::E, ci, outerWidth_[0], false);
         ci.expand(Cell(0, -outerWidth_[1], -outerWidth_[2]));
         layers_.push_back(ci);
         pdfs->getSliceBeforeGhostLayer(stencil::W, ci, outerWidth_[0], false);
         ci.expand(Cell(0, -outerWidth_[1], -outerWidth_[2]));
         layers_.push_back(ci);
      }

    
      for( auto & ci: layers_ )
      {
         streamCellInterval(pdfs, pdfs_tmp, ci);
      }
    

    pdfs->swapDataPointers(pdfs_tmp);

   }
   

   std::function<void (IBlock *)> streamOnlyNoAdvancement()
   {
      return [this](IBlock* block) { streamOnlyNoAdvancement(block); };
   }

   std::function<void (IBlock *)> streamOnlyNoAdvancement(Type type)
   {
      switch (type)
      {
         case Type::INNER:
            return [this](IBlock* block) { streamOnlyNoAdvancementInner(block); };
         case Type::OUTER:
            return [this](IBlock* block) { streamOnlyNoAdvancementOuter(block); };
         default:
            return [this](IBlock* block) { streamOnlyNoAdvancement(block); };
      }
   }

   std::function<void (IBlock *)> streamOnlyNoAdvancement(Type type, const cell_idx_t ghost_layers)
   {
      switch (type)
      {
         case Type::INNER:
            return [this](IBlock* block) { streamOnlyNoAdvancementInner(block); };
         case Type::OUTER:
            return [this](IBlock* block) { streamOnlyNoAdvancementOuter(block); };
         default:
            return [this, ghost_layers](IBlock* block) { streamOnlyNoAdvancement(block, ghost_layers); };
      }
   }

   

   void streamOnlyNoAdvancement(IBlock * block)
   {
      const cell_idx_t ghost_layers = 0;
      

      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      
      
      streamOnlyNoAdvancement(pdfs, pdfs_tmp, ghost_layers);
      
   }

   void streamOnlyNoAdvancement(IBlock * block, const cell_idx_t ghost_layers)
   {
      

      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      
      
      streamOnlyNoAdvancement(pdfs, pdfs_tmp, ghost_layers);
      
   }

   

   void streamOnlyNoAdvancementCellInterval(IBlock * block, const CellInterval & ci)
   {
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      
      
      streamOnlyNoAdvancementCellInterval(pdfs, pdfs_tmp, ci);
      
   }

   void streamOnlyNoAdvancementInner(IBlock * block)
   {
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      
      

      CellInterval inner = pdfs->xyzSize();
      inner.expand(Cell(-outerWidth_[0], -outerWidth_[1], -outerWidth_[2]));

      streamOnlyNoAdvancementCellInterval(pdfs, pdfs_tmp, inner);
   }

   void streamOnlyNoAdvancementOuter(IBlock * block)
   {

      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      field::GhostLayerField<double, 27> * pdfs_tmp;
      {
          // Getting temporary field pdfs_tmp
          auto it = cache_pdfs_.find( pdfs );
          if( it != cache_pdfs_.end() )
          {
              pdfs_tmp = *it;
          }
          else
          {
              pdfs_tmp = pdfs->cloneUninitialized();
              cache_pdfs_.insert(pdfs_tmp);
          }
      }

      
      

      if( layers_.empty() )
      {
         CellInterval ci;

         pdfs->getSliceBeforeGhostLayer(stencil::T, ci, outerWidth_[2], false);
         layers_.push_back(ci);
         pdfs->getSliceBeforeGhostLayer(stencil::B, ci, outerWidth_[2], false);
         layers_.push_back(ci);

         pdfs->getSliceBeforeGhostLayer(stencil::N, ci, outerWidth_[1], false);
         ci.expand(Cell(0, 0, -outerWidth_[2]));
         layers_.push_back(ci);
         pdfs->getSliceBeforeGhostLayer(stencil::S, ci, outerWidth_[1], false);
         ci.expand(Cell(0, 0, -outerWidth_[2]));
         layers_.push_back(ci);

         pdfs->getSliceBeforeGhostLayer(stencil::E, ci, outerWidth_[0], false);
         ci.expand(Cell(0, -outerWidth_[1], -outerWidth_[2]));
         layers_.push_back(ci);
         pdfs->getSliceBeforeGhostLayer(stencil::W, ci, outerWidth_[0], false);
         ci.expand(Cell(0, -outerWidth_[1], -outerWidth_[2]));
         layers_.push_back(ci);
      }

    
      for( auto & ci: layers_ )
      {
         streamOnlyNoAdvancementCellInterval(pdfs, pdfs_tmp, ci);
      }
    

    
   }
   

   std::function<void (IBlock *)> initialise()
   {
      return [this](IBlock* block) { initialise(block); };
   }

   std::function<void (IBlock *)> initialise(Type type)
   {
      switch (type)
      {
         case Type::INNER:
            return [this](IBlock* block) { initialiseInner(block); };
         case Type::OUTER:
            return [this](IBlock* block) { initialiseOuter(block); };
         default:
            return [this](IBlock* block) { initialise(block); };
      }
   }

   std::function<void (IBlock *)> initialise(Type type, const cell_idx_t ghost_layers)
   {
      switch (type)
      {
         case Type::INNER:
            return [this](IBlock* block) { initialiseInner(block); };
         case Type::OUTER:
            return [this](IBlock* block) { initialiseOuter(block); };
         default:
            return [this, ghost_layers](IBlock* block) { initialise(block, ghost_layers); };
      }
   }

   

   void initialise(IBlock * block)
   {
      const cell_idx_t ghost_layers = 0;
      

      auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      auto density = block->getData< field::GhostLayerField<double, 1> >(densityID);

      
      
      initialise(density, pdfs, velocity, ghost_layers);
      
   }

   void initialise(IBlock * block, const cell_idx_t ghost_layers)
   {
      

      auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      auto density = block->getData< field::GhostLayerField<double, 1> >(densityID);

      
      
      initialise(density, pdfs, velocity, ghost_layers);
      
   }

   

   void initialiseCellInterval(IBlock * block, const CellInterval & ci)
   {
      auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      auto density = block->getData< field::GhostLayerField<double, 1> >(densityID);

      
      
      initialiseCellInterval(density, pdfs, velocity, ci);
      
   }

   void initialiseInner(IBlock * block)
   {
      auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      auto density = block->getData< field::GhostLayerField<double, 1> >(densityID);

      
      

      CellInterval inner = density->xyzSize();
      inner.expand(Cell(-outerWidth_[0], -outerWidth_[1], -outerWidth_[2]));

      initialiseCellInterval(density, pdfs, velocity, inner);
   }

   void initialiseOuter(IBlock * block)
   {

      auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      auto density = block->getData< field::GhostLayerField<double, 1> >(densityID);

      
      

      if( layers_.empty() )
      {
         CellInterval ci;

         density->getSliceBeforeGhostLayer(stencil::T, ci, outerWidth_[2], false);
         layers_.push_back(ci);
         density->getSliceBeforeGhostLayer(stencil::B, ci, outerWidth_[2], false);
         layers_.push_back(ci);

         density->getSliceBeforeGhostLayer(stencil::N, ci, outerWidth_[1], false);
         ci.expand(Cell(0, 0, -outerWidth_[2]));
         layers_.push_back(ci);
         density->getSliceBeforeGhostLayer(stencil::S, ci, outerWidth_[1], false);
         ci.expand(Cell(0, 0, -outerWidth_[2]));
         layers_.push_back(ci);

         density->getSliceBeforeGhostLayer(stencil::E, ci, outerWidth_[0], false);
         ci.expand(Cell(0, -outerWidth_[1], -outerWidth_[2]));
         layers_.push_back(ci);
         density->getSliceBeforeGhostLayer(stencil::W, ci, outerWidth_[0], false);
         ci.expand(Cell(0, -outerWidth_[1], -outerWidth_[2]));
         layers_.push_back(ci);
      }

    
      for( auto & ci: layers_ )
      {
         initialiseCellInterval(density, pdfs, velocity, ci);
      }
    

    
   }
   

   std::function<void (IBlock *)> calculateMacroscopicParameters()
   {
      return [this](IBlock* block) { calculateMacroscopicParameters(block); };
   }

   std::function<void (IBlock *)> calculateMacroscopicParameters(Type type)
   {
      switch (type)
      {
         case Type::INNER:
            return [this](IBlock* block) { calculateMacroscopicParametersInner(block); };
         case Type::OUTER:
            return [this](IBlock* block) { calculateMacroscopicParametersOuter(block); };
         default:
            return [this](IBlock* block) { calculateMacroscopicParameters(block); };
      }
   }

   std::function<void (IBlock *)> calculateMacroscopicParameters(Type type, const cell_idx_t ghost_layers)
   {
      switch (type)
      {
         case Type::INNER:
            return [this](IBlock* block) { calculateMacroscopicParametersInner(block); };
         case Type::OUTER:
            return [this](IBlock* block) { calculateMacroscopicParametersOuter(block); };
         default:
            return [this, ghost_layers](IBlock* block) { calculateMacroscopicParameters(block, ghost_layers); };
      }
   }

   

   void calculateMacroscopicParameters(IBlock * block)
   {
      const cell_idx_t ghost_layers = 0;
      

      auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      auto density = block->getData< field::GhostLayerField<double, 1> >(densityID);

      
      
      calculateMacroscopicParameters(density, pdfs, velocity, ghost_layers);
      
   }

   void calculateMacroscopicParameters(IBlock * block, const cell_idx_t ghost_layers)
   {
      

      auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      auto density = block->getData< field::GhostLayerField<double, 1> >(densityID);

      
      
      calculateMacroscopicParameters(density, pdfs, velocity, ghost_layers);
      
   }

   

   void calculateMacroscopicParametersCellInterval(IBlock * block, const CellInterval & ci)
   {
      auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      auto density = block->getData< field::GhostLayerField<double, 1> >(densityID);

      
      
      calculateMacroscopicParametersCellInterval(density, pdfs, velocity, ci);
      
   }

   void calculateMacroscopicParametersInner(IBlock * block)
   {
      auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      auto density = block->getData< field::GhostLayerField<double, 1> >(densityID);

      
      

      CellInterval inner = density->xyzSize();
      inner.expand(Cell(-outerWidth_[0], -outerWidth_[1], -outerWidth_[2]));

      calculateMacroscopicParametersCellInterval(density, pdfs, velocity, inner);
   }

   void calculateMacroscopicParametersOuter(IBlock * block)
   {

      auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
      auto pdfs = block->getData< field::GhostLayerField<double, 27> >(pdfsID);
      auto density = block->getData< field::GhostLayerField<double, 1> >(densityID);

      
      

      if( layers_.empty() )
      {
         CellInterval ci;

         density->getSliceBeforeGhostLayer(stencil::T, ci, outerWidth_[2], false);
         layers_.push_back(ci);
         density->getSliceBeforeGhostLayer(stencil::B, ci, outerWidth_[2], false);
         layers_.push_back(ci);

         density->getSliceBeforeGhostLayer(stencil::N, ci, outerWidth_[1], false);
         ci.expand(Cell(0, 0, -outerWidth_[2]));
         layers_.push_back(ci);
         density->getSliceBeforeGhostLayer(stencil::S, ci, outerWidth_[1], false);
         ci.expand(Cell(0, 0, -outerWidth_[2]));
         layers_.push_back(ci);

         density->getSliceBeforeGhostLayer(stencil::E, ci, outerWidth_[0], false);
         ci.expand(Cell(0, -outerWidth_[1], -outerWidth_[2]));
         layers_.push_back(ci);
         density->getSliceBeforeGhostLayer(stencil::W, ci, outerWidth_[0], false);
         ci.expand(Cell(0, -outerWidth_[1], -outerWidth_[2]));
         layers_.push_back(ci);
      }

    
      for( auto & ci: layers_ )
      {
         calculateMacroscopicParametersCellInterval(density, pdfs, velocity, ci);
      }
    

    
   }
   

   

   private:
      shared_ptr< StructuredBlockStorage > blocks_;
      BlockDataID pdfsID;
    BlockDataID densityID;
    BlockDataID velocityID;
    double omega_;

    private: std::set< field::GhostLayerField<double, 27> *, field::SwapableCompare< field::GhostLayerField<double, 27> * > > cache_pdfs_;

      Cell outerWidth_;
      std::vector<CellInterval> layers_;

      
};


} // namespace lbm
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif