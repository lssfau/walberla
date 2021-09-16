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
//! \file GhostLayerField.impl.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Definitions of Field members
//
//======================================================================================================================

#include "field/iterators/IteratorMacros.h"
#include "field/GhostRegions.h"

#include "core/debug/Debug.h"


namespace walberla {
namespace field {


   //===================================================================================================================
   //
   //  CONSTRUCTION
   //
   //===================================================================================================================


   //*******************************************************************************************************************
   /*!\brief Creates a field of zero size
    *
    * To use this field call the init function
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   GhostLayerField<T,fSize_>::GhostLayerField( )
      : gl_(0)
   {
   }



   //*******************************************************************************************************************
   /*!\brief Creates an uninitialized field of given size
    *
    * \param _xSize size of x dimension without ghost layers
    * \param _ySize size of y dimension without ghost layers
    * \param _zSize size of z dimension without ghost layers
    * \param gl     number of ghost layers
    * \param l      memory layout of the field (see Layout)
    * \param alloc  class that describes how to allocate memory for the field, see FieldAllocator
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   GhostLayerField<T,fSize_>::GhostLayerField( uint_t _xSize, uint_t _ySize, uint_t _zSize, uint_t gl,
                                               const Layout & l, const shared_ptr<FieldAllocator<T> > &alloc)
   {
      init( _xSize, _ySize, _zSize, gl, l, alloc );
   }


   //*******************************************************************************************************************
   /*!\brief Creates a field and initializes it with constant value
    *
    * \param _xSize  size of x dimension without ghost layers
    * \param _ySize  size of y dimension without ghost layers
    * \param _zSize  size of z dimension without ghost layers
    * \param gl      number of ghost layers
    * \param initVal every element of the field is set to initVal (also ghost layers)
    * \param l       memory layout of the field (see Layout)
    * \param alloc  class that describes how to allocate memory for the field, see FieldAllocator
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   GhostLayerField<T,fSize_>::GhostLayerField( uint_t _xSize, uint_t _ySize, uint_t _zSize, uint_t gl,
                                               const T & initVal, const Layout & l,
                                               const shared_ptr<FieldAllocator<T> > &alloc    )
   {
      init( _xSize, _ySize, _zSize, gl, l, alloc );
      setWithGhostLayer( initVal );
   }


   //*******************************************************************************************************************
   /*!\brief Creates a field and initializes f coordinate with vector values
    *
    * \param _xSize  size of x dimension
    * \param _ySize  size of y dimension
    * \param _zSize  size of z dimension
    * \param gl      number of ghost layers
    * \param fValues initializes f coordinate with values from vector (see set(std::vector&) ) also ghost layers
    * \param l       memory layout of the field (see Layout)
    * \param alloc   class that describes how to allocate memory for the field, see FieldAllocator
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   GhostLayerField<T,fSize_>::GhostLayerField( uint_t _xSize, uint_t _ySize, uint_t _zSize, uint_t gl,
                                               const std::vector<T> & fValues, const Layout & l,
                                               const shared_ptr<FieldAllocator<T> > &alloc)
   {
      init( _xSize, _ySize, _zSize, gl, l, alloc );
      setWithGhostLayer( fValues );
   }



   //*******************************************************************************************************************
   /*!\brief Initializes a field, must be called exactly once.
    *
    * Is automatically called by constructors that take at least one argument
    *
    * \param _xSize size of x dimension without ghost layers
    * \param _ySize size of y dimension without ghost layers
    * \param _zSize size of z dimension without ghost layers
    * \param gl     number of ghost layers
    * \param l      memory layout of the field (see Layout)
    * \param alloc  class that describes how to allocate memory for the field, see FieldAllocator
    *******************************************************************************************************************/
    template<typename T, uint_t fSize_>
    void GhostLayerField<T,fSize_>::init( uint_t _xSize, uint_t _ySize, uint_t _zSize, uint_t gl,
                                          const Layout & l, const shared_ptr<FieldAllocator<T> > &alloc)
    {
       gl_ = gl;
       uint_t innerGhostLayerSize = ( l == fzyx ) ? gl : uint_t(0);
       Field<T,fSize_>::init( _xSize + 2*gl ,
                              _ySize + 2*gl,
                              _zSize + 2*gl, l, alloc,
                              innerGhostLayerSize );

       Field<T,fSize_>::setOffsets( gl, _xSize,
                                    gl, _ySize,
                                    gl, _zSize );
    }



   //*******************************************************************************************************************
   /*!\brief Deletes all stored data, and resizes the field
    *
    *  The resized field is uninitialized.
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   void GhostLayerField<T,fSize_>::resize( uint_t _xSize, uint_t _ySize, uint_t _zSize )
   {
      if ( _xSize == this->xSize() && _ySize == this->ySize() && _zSize == this->zSize()  )
         return;

      Field<T,fSize_>::resize( _xSize+2*gl_, _ySize+2*gl_, _zSize+2*gl_);
      Field<T,fSize_>::setOffsets( gl_, _xSize, gl_, _ySize, gl_, _zSize );
   }



   //*******************************************************************************************************************
   /*!\brief Deletes all stored data, and resizes the field
    *
    *  The resized field is uninitialized.
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   void GhostLayerField<T,fSize_>::resize( uint_t _xSize, uint_t _ySize, uint_t _zSize, uint_t _gl )
   {
      if ( _xSize == this->xSize() && _ySize == this->ySize() && _zSize == this->zSize() && _gl == gl_ )
         return;

      gl_ = _gl;
      Field<T,fSize_>::resize( _xSize+2*gl_, _ySize+2*gl_, _zSize+2*gl_);
      Field<T,fSize_>::setOffsets( gl_, _xSize, gl_, _ySize, gl_, _zSize );
   }






   //===================================================================================================================
   //
   //  ELEMENT ACCESS
   //
   //===================================================================================================================


   //*******************************************************************************************************************
   /*!\brief Sets all entries (including the ghost layer) of the field to given value
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   void GhostLayerField<T,fSize_>::setWithGhostLayer (const T & value)
   {
#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wundefined-bool-conversion"
#endif
      // take care of proper thread<->memory assignment (first-touch allocation policy !)
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( this,

         for( uint_t f = uint_t(0); f < fSize_; ++f )
            this->get(x,y,z,f) = value;
      
      ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif
   }

   //*******************************************************************************************************************
   /*!\brief Initializes the f coordinate to values from vector, in all cells including the ghost layers
    * Sets the entry (x,y,z,f) to fValues[f]
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   void GhostLayerField<T,fSize_>::setWithGhostLayer (const std::vector<T> & fValues)
   {
      WALBERLA_ASSERT(fValues.size() == fSize_);
          
#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wundefined-bool-conversion"
#endif
      // take care of proper thread<->memory assignment (first-touch allocation policy !)
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( this,

         for( uint_t f = uint_t(0); f < fSize_; ++f )
            this->get(x,y,z,f) = fValues[f];
      
      ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif
   }






   //===================================================================================================================
   //
   //  ITERATORS
   //
   //===================================================================================================================



   //*******************************************************************************************************************
   /*!\brief Iterator over all cells, including the ghost layers
    *
    * same as begin() , but with ghost layer
    *
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::iterator
   GhostLayerField<T,fSize_>::beginWithGhostLayer( )
   {
      return beginWithGhostLayer( cell_idx_c( gl_ ) );
   }

   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::iterator
   GhostLayerField<T,fSize_>::beginWithGhostLayer( cell_idx_t numGhostLayers )
   {
      WALBERLA_ASSERT_LESS_EQUAL( numGhostLayers, cell_idx_c( gl_ )  );
      const uint_t xs = Field<T,fSize_>::xSize() + 2 * uint_c( numGhostLayers );
      const uint_t ys = Field<T,fSize_>::ySize() + 2 * uint_c( numGhostLayers );
      const uint_t zs = Field<T,fSize_>::zSize() + 2 * uint_c( numGhostLayers );

      return iterator( this,
                      -numGhostLayers,-numGhostLayers,-numGhostLayers,0,
                       xs, ys, zs, fSize_ );
   }



   //*******************************************************************************************************************
   /*!\brief Returns const_iterator, see beginWithGhostLayer()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::const_iterator
   GhostLayerField<T,fSize_>::beginWithGhostLayer( ) const
   {
      return beginWithGhostLayer( cell_idx_c( gl_ ) );
   }

   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::const_iterator
   GhostLayerField<T,fSize_>::beginWithGhostLayer( cell_idx_t numGhostLayers ) const
   {
      const uint_t xs = Field<T,fSize_>::xSize() + 2 * uint_c( numGhostLayers );
      const uint_t ys = Field<T,fSize_>::ySize() + 2 * uint_c( numGhostLayers );
      const uint_t zs = Field<T,fSize_>::zSize() + 2 * uint_c( numGhostLayers );

      return const_iterator(  this,
                              -numGhostLayers,-numGhostLayers,-numGhostLayers,0,
                              xs, ys, zs, fSize_ );

   }



   //*******************************************************************************************************************
   /*!\brief Iterates only over all cells including ghost layers of XYZ coordinate, f is always 0
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::iterator
   GhostLayerField<T,fSize_>::beginWithGhostLayerXYZ( )
   {
      return beginWithGhostLayerXYZ( cell_idx_c( gl_ ) );
   }

   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::iterator
   GhostLayerField<T,fSize_>::beginWithGhostLayerXYZ( cell_idx_t numGhostLayers )
   {
      const uint_t xs = Field<T,fSize_>::xSize() + 2 * uint_c( numGhostLayers );
      const uint_t ys = Field<T,fSize_>::ySize() + 2 * uint_c( numGhostLayers );
      const uint_t zs = Field<T,fSize_>::zSize() + 2 * uint_c( numGhostLayers );

      return iterator( this,
                        -numGhostLayers,-numGhostLayers,-numGhostLayers,0,
                        xs, ys, zs, 1 );
   }

   //*******************************************************************************************************************
   /*!\brief Const version of beginWithGhostLayerXYZ()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::const_iterator
   GhostLayerField<T,fSize_>::beginWithGhostLayerXYZ( ) const
   {
      return beginWithGhostLayerXYZ( cell_idx_c( gl_ ) );
   }

   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::const_iterator
   GhostLayerField<T,fSize_>::beginWithGhostLayerXYZ( cell_idx_t numGhostLayers ) const
   {
      const uint_t xs = Field<T,fSize_>::xSize() + 2 * uint_c( numGhostLayers );
      const uint_t ys = Field<T,fSize_>::ySize() + 2 * uint_c( numGhostLayers );
      const uint_t zs = Field<T,fSize_>::zSize() + 2 * uint_c( numGhostLayers );

      return const_iterator( this,
                             -numGhostLayers,-numGhostLayers,-numGhostLayers,0,
                             xs, ys, zs, 1 );
   }

   template<typename T, uint_t fSize_>
   void GhostLayerField<T,fSize_>::getGhostRegion(stencil::Direction d, CellInterval & ci,
                                                  cell_idx_t thickness, bool fullSlice ) const
   {
      ci = field::getGhostRegion( *this, d, thickness, fullSlice );
   }

   template<typename T, uint_t fSize_>
   void GhostLayerField<T,fSize_>::getSliceBeforeGhostLayer(stencil::Direction d, CellInterval & ci,
                                                            cell_idx_t thickness, bool fullSlice ) const
   {
      ci = field::getSliceBeforeGhostLayer( *this, d, thickness, fullSlice );
   }

   //*******************************************************************************************************************
   /*!\brief Checks if a given cell is in the inner part of the field ( not in ghost region or outside )
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   bool GhostLayerField<T,fSize_>::isInInnerPart( const Cell & cell ) const
   {
      return !(cell[0] < 0 ||
           cell[1] < 0 ||
           cell[2] < 0 ||
           cell[0] >= cell_idx_c( this->xSize() ) ||
           cell[1] >= cell_idx_c( this->ySize() ) ||
           cell[2] >= cell_idx_c( this->zSize() ));
   }

   //*******************************************************************************************************************
   /*!\brief Iterates only over ghost layers of a given direction
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::iterator
   GhostLayerField<T,fSize_>::beginGhostLayerOnly( stencil::Direction dir, bool fullSlice )
   {
      CellInterval ci;
      getGhostRegion(dir,ci, cell_idx_c(gl_), fullSlice );

      return ForwardFieldIterator<T,fSize_>( this,
                                      ci.xMin(),ci.yMin(), ci.zMin(), 0,
                                      ci.xSize(), ci.ySize(), ci.zSize(), fSize_ );
   }


   //*******************************************************************************************************************
   /*!\brief Const version of beginGhostLayersOnly()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::const_iterator
   GhostLayerField<T,fSize_>::beginGhostLayerOnly( stencil::Direction dir, bool fullSlice ) const
   {
      CellInterval ci;
      getGhostRegion(dir,ci, cell_idx_c(gl_), fullSlice );

      return ForwardFieldIterator<const T,fSize_>( this,
                                            ci.xMin(),ci.yMin(), ci.zMin(), 0,
                                            ci.xSize(), ci.ySize(), ci.zSize(), fSize_ );
   }


   //*******************************************************************************************************************
   /*!\brief Iterates only over specified number of ghost layers of a given direction
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::iterator
   GhostLayerField<T,fSize_>::beginGhostLayerOnly( uint_t thickness, stencil::Direction dir, bool fullSlice )
   {
      CellInterval ci;
      getGhostRegion( dir, ci, cell_idx_c(thickness), fullSlice );

      return ForwardFieldIterator<T,fSize_>( this,
                                      ci.xMin(),ci.yMin(), ci.zMin(), 0,
                                      ci.xSize(), ci.ySize(), ci.zSize(), fSize_ );
   }


   //*******************************************************************************************************************
   /*!\brief Const version of beginGhostLayersOnly(uint_t thickness, stencil::Direction)
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::const_iterator
   GhostLayerField<T,fSize_>::beginGhostLayerOnly(  uint_t thickness, stencil::Direction dir, bool fullSlice ) const
   {
      CellInterval ci;
      getGhostRegion( dir, ci, cell_idx_c(thickness), fullSlice );

      return ForwardFieldIterator<const T,fSize_>( this,
                                            ci.xMin(),ci.yMin(), ci.zMin(), 0,
                                            ci.xSize(), ci.ySize(), ci.zSize(), fSize_ );
   }



   //*******************************************************************************************************************
   /*!\brief Iterates only over ghost layers of a given direction, only over xyz coordinates, f is fixed
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::iterator
   GhostLayerField<T,fSize_>::beginGhostLayerOnlyXYZ( stencil::Direction dir, cell_idx_t f, bool fullSlice )
   {
      CellInterval ci;
      getGhostRegion( dir, ci, cell_idx_c(gl_) , fullSlice );

      return ForwardFieldIterator<T,fSize_>( this,
                                      ci.xMin(),ci.yMin(), ci.zMin(), f,
                                      ci.xSize(), ci.ySize(), ci.zSize(), uint_c(f+1) );
   }



   //*******************************************************************************************************************
   /*!\brief Const version of beginGhostLayersOnlyXYZ()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::const_iterator
   GhostLayerField<T,fSize_>::beginGhostLayerOnlyXYZ( stencil::Direction dir, cell_idx_t f, bool fullSlice ) const
   {
      CellInterval ci;
      getGhostRegion(dir, ci, cell_idx_c(gl_), fullSlice );

      return ForwardFieldIterator<const T,fSize_>( this,
                                                   ci.xMin(),ci.yMin(), ci.zMin(), f,
                                                   ci.xSize(), ci.ySize(), ci.zSize(), uint_c( f + cell_idx_t(1) ) );
   }

   //*******************************************************************************************************************
   /*!\brief Iterates only over ghost layers of a given direction, only over xyz coordinates, f is fixed
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::iterator
   GhostLayerField<T,fSize_>::beginGhostLayerOnlyXYZ( uint_t thickness, stencil::Direction dir, cell_idx_t f, bool fullSlice )
   {
      CellInterval ci;
      getGhostRegion( dir, ci, cell_idx_c(thickness) , fullSlice );

      return ForwardFieldIterator<T,fSize_>( this,
                                      ci.xMin(),ci.yMin(), ci.zMin(), f,
                                      ci.xSize(), ci.ySize(), ci.zSize(), uint_c( f + cell_idx_t(1) ) );
   }



   //*******************************************************************************************************************
   /*!\brief Const version of beginGhostLayersOnlyXYZ()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::const_iterator
   GhostLayerField<T,fSize_>::beginGhostLayerOnlyXYZ( uint_t thickness, stencil::Direction dir, cell_idx_t f, bool fullSlice ) const
   {
      CellInterval ci;
      getGhostRegion(dir, ci, cell_idx_c(thickness), fullSlice );

      return ForwardFieldIterator<const T,fSize_>( this,
                                                   ci.xMin(),ci.yMin(), ci.zMin(), f,
                                                   ci.xSize(), ci.ySize(), ci.zSize(), uint_c( f + cell_idx_t(1) ) );
   }


   //*******************************************************************************************************************
   /*!\brief Iterates only over the last slice before ghost layer
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::iterator
   GhostLayerField<T,fSize_>::beginSliceBeforeGhostLayer( stencil::Direction dir, cell_idx_t thickness, bool fullSlice )
   {
      CellInterval ci;
      getSliceBeforeGhostLayer(dir, ci, thickness, fullSlice );

      return ForwardFieldIterator<T,fSize_>( this,
                                      ci.xMin(),ci.yMin(), ci.zMin(), 0,
                                      ci.xSize(), ci.ySize(), ci.zSize(), fSize_ );
   }


   //*******************************************************************************************************************
   /*!\brief Const version of beginSliceBeforeGhostLayer()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::const_iterator
   GhostLayerField<T,fSize_>::beginSliceBeforeGhostLayer( stencil::Direction dir, cell_idx_t thickness, bool fullSlice ) const
   {
      CellInterval ci;
      getSliceBeforeGhostLayer( dir, ci, thickness, fullSlice );

      return ForwardFieldIterator<const T,fSize_>( this,
                                            ci.xMin(),ci.yMin(), ci.zMin(), 0,
                                            ci.xSize(), ci.ySize(), ci.zSize(), fSize_ );
   }


   //*******************************************************************************************************************
   /*!\brief Iterates only over the last slice before ghost layer, only in XYZ direction, f is fixed
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::iterator
   GhostLayerField<T,fSize_>::beginSliceBeforeGhostLayerXYZ( stencil::Direction dir, cell_idx_t thickness,
                                                             cell_idx_t f, bool fullSlice  )
   {
      CellInterval ci;
      getSliceBeforeGhostLayer(dir, ci, thickness, fullSlice );

      return ForwardFieldIterator<T,fSize_>( this,
                                      ci.xMin(),ci.yMin(), ci.zMin(), f,
                                      ci.xSize(), ci.ySize(), ci.zSize(), uint_c( f + cell_idx_t(1) ) );
   }


   //*******************************************************************************************************************
   /*!\brief Const version of beginSliceBeforeGhostLayerXYZ()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::const_iterator
   GhostLayerField<T,fSize_>::beginSliceBeforeGhostLayerXYZ( stencil::Direction dir,  cell_idx_t thickness,
                                                             cell_idx_t f, bool fullSlice ) const
   {
      CellInterval ci;
      getSliceBeforeGhostLayer( dir, ci, thickness, fullSlice );

      return ForwardFieldIterator<const T,fSize_>( this,
                                            ci.xMin(),ci.yMin(), ci.zMin(), f,
                                            ci.xSize(), ci.ySize(), ci.zSize(), uint_c( f + cell_idx_t(1) ) );
   }

   //*******************************************************************************************************************
   /*!\brief Returns the x/y/z Size of the field with ghost layers
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline CellInterval GhostLayerField<T,fSize_>::xyzSizeWithGhostLayer() const
   {
      CellInterval ci = Field<T,fSize_>::xyzSize();
      for( uint_t i=0; i<3; ++i ) {
         ci.min()[i] -= cell_idx_c( gl_ );
         ci.max()[i] += cell_idx_c( gl_ );
      }
      return ci;
   }








   //===================================================================================================================
   //
   //  REVERSE ITERATORS
   //
   //===================================================================================================================



   //*******************************************************************************************************************
   /*!\brief Reverse Iterator over all cells, including the ghost layers
    *
    * same as rbegin() , but with ghost layer
    *
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::reverse_iterator
   GhostLayerField<T,fSize_>::rbeginWithGhostLayer()
   {
      const uint_t xs = Field<T,fSize_>::xSize() + 2*gl_;
      const uint_t ys = Field<T,fSize_>::ySize() + 2*gl_;
      const uint_t zs = Field<T,fSize_>::zSize() + 2*gl_;

      return reverse_iterator( this,
                              -cell_idx_c(gl_),-cell_idx_c(gl_),-cell_idx_c(gl_),0,
                               xs, ys, zs, fSize_ );
   }


   //*******************************************************************************************************************
   /*!\brief Returns const_iterator, see beginWithGhostLayer()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::const_reverse_iterator
   GhostLayerField<T,fSize_>::rbeginWithGhostLayer() const
   {
      const uint_t xs = Field<T,fSize_>::xSize() + 2*gl_;
      const uint_t ys = Field<T,fSize_>::ySize() + 2*gl_;
      const uint_t zs = Field<T,fSize_>::zSize() + 2*gl_;

      return const_reverse_iterator(  this,
                                      -cell_idx_c(gl_),-cell_idx_c(gl_),-cell_idx_c(gl_),0,
                                      xs, ys, zs, fSize_ );

   }


   //*******************************************************************************************************************
   /*!\brief Iterates only over all cells including ghost layers of XYZ coordinate, f is always 0
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::reverse_iterator
   GhostLayerField<T,fSize_>::rbeginWithGhostLayerXYZ()
   {
      const uint_t xs = Field<T,fSize_>::xSize() + 2*gl_;
      const uint_t ys = Field<T,fSize_>::ySize() + 2*gl_;
      const uint_t zs = Field<T,fSize_>::zSize() + 2*gl_;

      return reverse_iterator( this,
                              -cell_idx_c(gl_),-cell_idx_c(gl_),-cell_idx_c(gl_),0,
                               xs, ys, zs, 1 );
   }

   //*******************************************************************************************************************
   /*!\brief Const version of beginWithGhostLayerXYZ()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename GhostLayerField<T,fSize_>::const_reverse_iterator
   GhostLayerField<T,fSize_>::rbeginWithGhostLayerXYZ() const
   {
      const uint_t xs = Field<T,fSize_>::xSize() + 2*gl_;
      const uint_t ys = Field<T,fSize_>::ySize() + 2*gl_;
      const uint_t zs = Field<T,fSize_>::zSize() + 2*gl_;

      return const_reverse_iterator( this,
                                     -cell_idx_c(gl_),-cell_idx_c(gl_),-cell_idx_c(gl_),0,
                                     xs, ys, zs, 1 );
   }







   //===================================================================================================================
   //
   //  SLICING AND CLONING
   //
   //===================================================================================================================

   //*******************************************************************************************************************
   /*!\brief Private copy constructor, which does a shallow copy
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   GhostLayerField<T,fSize_>::GhostLayerField(const GhostLayerField<T,fSize_> & other)
      : Field<T,fSize_>::Field(other),
        gl_( other.gl_ )
   {
   }

   //*******************************************************************************************************************
   /*!\brief Private copy constructor, which does a flattened shallow copy
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   template <typename T2, uint_t fSize2>
   GhostLayerField<T,fSize_>::GhostLayerField(const GhostLayerField<T2,fSize2> & other)
      : Field<T,fSize_>::Field(other),
        gl_( other.gl_ )
   {
   }

   //*******************************************************************************************************************
   /*!\brief See Field::cloneShallowCopyInternal()
    *  Has to be re-implemented because a new GhostLayerField is created
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   Field<T,fSize_> * GhostLayerField<T,fSize_>::cloneShallowCopyInternal() const
   {
      return new GhostLayerField<T,fSize_>(*this);
   }

   //*******************************************************************************************************************
   /*!\brief See Field::flattenedShallowCopyInternal()
    *  Has to be re-implemented because a new GhostLayerField is created
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   typename Field<T,fSize_>::FlattenedField * GhostLayerField<T,fSize_>::flattenedShallowCopyInternal() const
   {
      return new GhostLayerField<T,fSize_>::FlattenedField(*this);
   }


   template<typename T, uint_t fSize_>
   GhostLayerField<T,fSize_> * GhostLayerField<T,fSize_>::clone() const
   {
      return dynamic_cast<GhostLayerField<T,fSize_>* > (Field<T,fSize_>::clone() );
   }

   template<typename T, uint_t fSize_>
   GhostLayerField<T,fSize_> * GhostLayerField<T,fSize_>::cloneUninitialized() const
   {
      return dynamic_cast<GhostLayerField<T,fSize_>* > (Field<T,fSize_>::cloneUninitialized() );
   }

   template<typename T, uint_t fSize_>
   GhostLayerField<T,fSize_> * GhostLayerField<T,fSize_>::cloneShallowCopy() const
   {
      return dynamic_cast<GhostLayerField<T,fSize_>* > (Field<T,fSize_>::cloneShallowCopy() );
   }

   template<typename T, uint_t fSize_>
   typename GhostLayerField<T,fSize_>::FlattenedField * GhostLayerField<T,fSize_>::flattenedShallowCopy() const
   {
      return dynamic_cast<GhostLayerField<T,fSize_>::FlattenedField* > (Field<T,fSize_>::flattenedShallowCopy() );
   }

   template<typename T, uint_t fSize_>
   void GhostLayerField<T,fSize_>::slice( const CellInterval & interval )
   {
      Field<T,fSize_>::slice( interval );

      // Assert that there is still space for ghost-layers after slicing

      WALBERLA_ASSERT_GREATER_EQUAL( uint_c(this->xOff()), gl_ );
      WALBERLA_ASSERT_GREATER_EQUAL( this->xAllocSize() - uint_c(this->xOff()) - this->xSize(), gl_ );
      WALBERLA_ASSERT_GREATER_EQUAL( uint_c(this->yOff()), gl_ );
      WALBERLA_ASSERT_GREATER_EQUAL( this->yAllocSize() - uint_c(this->yOff()) - this->ySize(), gl_ );
      WALBERLA_ASSERT_GREATER_EQUAL( uint_c(this->zOff()), gl_ );
      WALBERLA_ASSERT_GREATER_EQUAL( this->zAllocSize() - uint_c(this->zOff()) - this->zSize(), gl_ );

   }

   template<typename T, uint_t fSize_>
   void GhostLayerField<T,fSize_>::shiftCoordinates( cell_idx_t cx, cell_idx_t cy, cell_idx_t cz )
   {
      Field<T,fSize_>::shiftCoordinates( cx, cy, cz );

      // Assert that there is still space for ghost-layers after slicing
      WALBERLA_ASSERT_GREATER_EQUAL( uint_c(this->xOff()), gl_ );
      WALBERLA_ASSERT_GREATER_EQUAL( this->xAllocSize() - uint_c(this->xOff()) - this->xSize(), gl_ );
      WALBERLA_ASSERT_GREATER_EQUAL( uint_c(this->yOff()), gl_ );
      WALBERLA_ASSERT_GREATER_EQUAL( this->yAllocSize() - uint_c(this->yOff()) - this->ySize(), gl_ );
      WALBERLA_ASSERT_GREATER_EQUAL( uint_c(this->zOff()), gl_ );
      WALBERLA_ASSERT_GREATER_EQUAL( this->zAllocSize() - uint_c(this->zOff()) - this->zSize(), gl_ );
   }



   template<typename T, uint_t fSize_>
   GhostLayerField<T,fSize_> * GhostLayerField<T,fSize_>::getSlicedField( const CellInterval & ci ) const
   {
      return dynamic_cast<GhostLayerField<T,fSize_> *>( Field<T,fSize_>::getSlicedField(ci) );
   }

}
}
