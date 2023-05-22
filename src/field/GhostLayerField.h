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
//! \file GhostLayerField.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief GhostLayerField class extends Field with ghost layer information
//
//======================================================================================================================

#pragma once

#include "Field.h"
#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "stencil/Directions.h"


namespace walberla {
namespace field {


   //*******************************************************************************************************************
   /*! Extends the Field with ghost-layer information
   *
   * All variants of the begin() function exist also in the "WithGhostLayer" variant
   * which iterate over the field including the ghost layers. There are also iterators that
   * go only over the ghost layer or, only over the last inner slice, which is useful when doing
   * ghost layer based communication.
   *
   * \ingroup field
   *
   * See also \ref fieldPage
   *
   */
   //*******************************************************************************************************************
   template<typename T, uint_t fSize_>
   class GhostLayerField : public Field<T,fSize_>
   {
   public:
      //** Type Definitions  *******************************************************************************************
      /*! \name Type Definitions */
      //@{
      using value_type = typename Field<T, fSize_>::value_type;

      using iterator = typename Field<T, fSize_>::iterator;
      using const_iterator = typename Field<T, fSize_>::const_iterator;

      using reverse_iterator = typename Field<T, fSize_>::reverse_iterator;
      using const_reverse_iterator = typename Field<T, fSize_>::const_reverse_iterator;

      using base_iterator = typename Field<T, fSize_>::base_iterator;
      using const_base_iterator = typename Field<T, fSize_>::const_base_iterator;

      using Ptr = typename Field<T, fSize_>::Ptr;
      using ConstPtr = typename Field<T, fSize_>::ConstPtr;

      using FlattenedField = typename std::conditional<VectorTrait<T>::F_SIZE != 0, GhostLayerField<typename VectorTrait<T>::OutputType, VectorTrait<T>::F_SIZE * fSize_>, GhostLayerField<T, fSize_>>::type;
      //@}
      //****************************************************************************************************************


      //**Construction & Destruction************************************************************************************
      /*! \name Construction & Destruction */
      //@{


      GhostLayerField( uint_t xSize, uint_t ySize, uint_t zSize, uint_t gl,
                      const Layout & layout = fzyx,
                      const shared_ptr<FieldAllocator<T> > &alloc = shared_ptr<FieldAllocator<T> >() );
      GhostLayerField( uint_t xSize, uint_t ySize, uint_t zSize, uint_t gl,
                       const T & initValue, const Layout & layout = fzyx,
                       const shared_ptr<FieldAllocator<T> > &alloc = shared_ptr<FieldAllocator<T> >() );
      GhostLayerField( uint_t xSize, uint_t ySize, uint_t zSize, uint_t gl,
                       const std::vector<T> & fValues, const Layout & layout = fzyx,
                       const shared_ptr<FieldAllocator<T> > &alloc = shared_ptr<FieldAllocator<T> >() );

      ~GhostLayerField() override = default;



      void init( uint_t xSizeWithoutGhostLayer,
                 uint_t ySizeWithoutGhostLayer,
                 uint_t zSizeWithoutGhostLayer,
                 uint_t nrGhostLayers,
                 const Layout & layout = fzyx,
                 const shared_ptr<FieldAllocator<T> > &alloc = shared_ptr<FieldAllocator<T> >() );


      void resize( uint_t xSize, uint_t ySize, uint_t zSize ) override;
              void resize( uint_t xSize, uint_t ySize, uint_t zSize, uint_t gl );

      using Field<T,fSize_>::resize;

      inline GhostLayerField<T,fSize_> * clone()              const;
      inline GhostLayerField<T,fSize_> * cloneUninitialized() const;
      inline GhostLayerField<T,fSize_> * cloneShallowCopy()   const;
      inline FlattenedField * flattenedShallowCopy() const;
      //@}
      //****************************************************************************************************************


      //** Size Information ********************************************************************************************
      /*! \name Size Information */
      //@{
      inline uint_t       xSizeWithGhostLayer()        const  { return Field<T,fSize_>::xSize() + uint_c(2)*gl_; }
      inline uint_t       ySizeWithGhostLayer()        const  { return Field<T,fSize_>::ySize() + uint_c(2)*gl_; }
      inline uint_t       zSizeWithGhostLayer()        const  { return Field<T,fSize_>::zSize() + uint_c(2)*gl_; }
      inline uint_t       sizeWithGhostLayer(uint_t i) const  { return i==3 ? fSize_ :
                                                                              Field<T,fSize_>::size(i) + uint_c(2)*gl_; }
      inline uint_t       nrOfGhostLayers()            const  { return gl_; }
      inline CellInterval xyzSizeWithGhostLayer()      const;
      //@}
      //****************************************************************************************************************


      //** Element Access **********************************************************************************************
      /*! \name Element Access */
      //@{
      void setWithGhostLayer (const T & value);
      void setWithGhostLayer (const std::vector<T> & fValues);
      //@}
      //****************************************************************************************************************


      //** Iterators  **************************************************************************************************
      /*! \name Iterators */
      //@{
            iterator beginWithGhostLayer();
      const_iterator beginWithGhostLayer() const;

            iterator beginWithGhostLayer( cell_idx_t numGhostLayers );
      const_iterator beginWithGhostLayer( cell_idx_t numGhostLayers ) const;

            iterator beginWithGhostLayerXYZ();
      const_iterator beginWithGhostLayerXYZ() const;

            iterator beginWithGhostLayerXYZ( cell_idx_t numGhostLayers );
      const_iterator beginWithGhostLayerXYZ( cell_idx_t numGhostLayers ) const;

            iterator beginGhostLayerOnly( stencil::Direction dir, bool fullSlice = false );
      const_iterator beginGhostLayerOnly( stencil::Direction dir, bool fullSlice = false ) const;

            iterator beginGhostLayerOnly( uint_t thickness, stencil::Direction dir, bool fullSlice = false );
      const_iterator beginGhostLayerOnly( uint_t thickness, stencil::Direction dir, bool fullSlice = false ) const;

            iterator beginGhostLayerOnlyXYZ( stencil::Direction dir, cell_idx_t f = 0, bool fullSlice = false );
      const_iterator beginGhostLayerOnlyXYZ( stencil::Direction dir, cell_idx_t f = 0, bool fullSlice = false ) const;

            iterator beginGhostLayerOnlyXYZ( uint_t thickness, stencil::Direction dir, cell_idx_t f = 0, bool fullSlice = false );
      const_iterator beginGhostLayerOnlyXYZ( uint_t thickness, stencil::Direction dir, cell_idx_t f = 0, bool fullSlice = false ) const;

            iterator beginSliceBeforeGhostLayer( stencil::Direction dir, cell_idx_t thickness = 1,
                                                 bool fullSlice = false );
      const_iterator beginSliceBeforeGhostLayer( stencil::Direction dir, cell_idx_t thickness = 1,
                                                 bool fullSlice = false ) const;

            iterator beginSliceBeforeGhostLayerXYZ( stencil::Direction dir, cell_idx_t thickness = 1,
                                                    cell_idx_t f = 0, bool fullSlice = false );
      const_iterator beginSliceBeforeGhostLayerXYZ( stencil::Direction dir, cell_idx_t thickness = 1,
                                                    cell_idx_t f = 0, bool fullSlice = false ) const;


      void getGhostRegion( stencil::Direction dir, CellInterval & ghostAreaOut, cell_idx_t thickness,
                           bool fullSlice = false ) const;
      void getSliceBeforeGhostLayer( stencil::Direction d, CellInterval & ci, cell_idx_t thickness=1,
                                     bool fullSlice = false ) const;
      bool isInInnerPart( const Cell & cell ) const;
      //@}
      //****************************************************************************************************************

      //** Reverse Iterators *******************************************************************************************
      /*! \name Reverse Iterators */
      //@{
            reverse_iterator rbeginWithGhostLayer();
      const_reverse_iterator rbeginWithGhostLayer() const;

            reverse_iterator rbeginWithGhostLayerXYZ();
      const_reverse_iterator rbeginWithGhostLayerXYZ() const;
      //@}
      //****************************************************************************************************************


      //** Slicing  ****************************************************************************************************
      /*! \name Slicing */
      //@{
      GhostLayerField<T,fSize_> * getSlicedField( const CellInterval & interval ) const;
      void slice           ( const CellInterval & interval ) override;
      void shiftCoordinates( cell_idx_t cx, cell_idx_t cy, cell_idx_t cz ) override;
      //@}
      //****************************************************************************************************************

      //** TimestepInformation *****************************************************************************************
      /*! \name TimestepCounter */
      //@{
      inline uint8_t advanceTimestep()
      {
         timestepCounter_ = (timestepCounter_ + 1) & 1;
         return timestepCounter_;
      }
      inline uint8_t getTimestep() const { return timestepCounter_; }
      inline uint8_t getTimestepPlusOne() const { return (timestepCounter_ + 1) & 1; }
      inline bool isEvenTimeStep() const {return (((timestepCounter_) &1) ^ 1); }
      //@}
      //****************************************************************************************************************

   protected:
      GhostLayerField( );


      uint_t gl_; ///< Number of ghost layers

      //** Shallow Copy ************************************************************************************************
      /*! \name Shallow Copy */
      //@{
      Field<T,fSize_> * cloneShallowCopyInternal()   const override;
      typename Field<T,fSize_>::FlattenedField * flattenedShallowCopyInternal() const override;
      GhostLayerField(const GhostLayerField<T,fSize_> & other);
      template <typename T2, uint_t fSize2>
      GhostLayerField(const GhostLayerField<T2, fSize2> & other);
      //@}
      //****************************************************************************************************************

      template <typename T2, uint_t fSize2>
      friend class GhostLayerField;

      uint8_t timestepCounter_;
   };

} // namespace field
} // namespace walberla

#include "GhostLayerField.impl.h"


//======================================================================================================================
//
//  EXPORTS
//
//======================================================================================================================

namespace walberla {
   // Export ghost layer field class to walberla namespace
   using field::GhostLayerField;
}
