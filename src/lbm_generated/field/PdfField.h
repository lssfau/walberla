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
//! \file PdfField.h
//! \ingroup lbm_generated
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "field/GhostLayerField.h"
#include "field/SwapableCompare.h"


namespace walberla::lbm_generated {

template< typename LatticeStorageSpecification_T >
class PdfField : public GhostLayerField< typename LatticeStorageSpecification_T::value_type, LatticeStorageSpecification_T::Stencil::Size >
{
public:

   //** Type Definitions  **********************************************************************************************
   /*! \name Type Definitions */
   //@{
   using LatticeStorageSpecification = LatticeStorageSpecification_T;
   using Stencil = typename LatticeStorageSpecification_T::Stencil;

   using value_type = typename LatticeStorageSpecification_T::value_type;

   using Ptr = typename GhostLayerField<value_type, Stencil::Size>::Ptr;
   using ConstPtr = typename GhostLayerField<value_type, Stencil::Size>::ConstPtr;
   //@}
   //*******************************************************************************************************************

   PdfField( const uint_t _xSize, const uint_t _ySize, const uint_t _zSize,
            const LatticeStorageSpecification_T & storageSpecification,
             const uint_t ghostLayers = uint_t(1), const field::Layout & _layout = field::zyxf,
             const shared_ptr< field::FieldAllocator<value_type> > & alloc = shared_ptr< field::FieldAllocator<value_type> >() );

   ~PdfField() override = default;

   inline PdfField * clone()              const;
   inline PdfField * cloneUninitialized() const;
   inline PdfField * cloneShallowCopy()   const;


   /////////////////////////////////////////////////
   // Access functions (with stencil::Direction!) //
   /////////////////////////////////////////////////

   using GhostLayerField< value_type, Stencil::Size >::get;

         value_type & get( cell_idx_t x, cell_idx_t y, cell_idx_t z, stencil::Direction d )       { return get( x, y, z, Stencil::idx[d] ); }
   const value_type & get( cell_idx_t x, cell_idx_t y, cell_idx_t z, stencil::Direction d ) const { return get( x, y, z, Stencil::idx[d] ); }
         value_type & get( const Cell & c, stencil::Direction d )       { return get( c.x(), c.y(), c.z(), Stencil::idx[d] ); }
   const value_type & get( const Cell & c, stencil::Direction d ) const { return get( c.x(), c.y(), c.z(), Stencil::idx[d] ); }

   using GhostLayerField< value_type, Stencil::Size >::operator();

         value_type & operator()( cell_idx_t x, cell_idx_t y, cell_idx_t z, stencil::Direction d )       { return get( x, y, z, Stencil::idx[d] ); }
   const value_type & operator()( cell_idx_t x, cell_idx_t y, cell_idx_t z, stencil::Direction d ) const { return get( x, y, z, Stencil::idx[d] ); }
         value_type & operator()( const Cell & c, stencil::Direction d )       { return get( c.x(), c.y(), c.z(), Stencil::idx[d] ); }
   const value_type & operator()( const Cell & c, stencil::Direction d ) const { return get( c.x(), c.y(), c.z(), Stencil::idx[d] ); }


protected:
   //** Shallow Copy ***************************************************************************************************
   /*! \name Shallow Copy */
   //@{
   inline PdfField( const PdfField< LatticeStorageSpecification_T > & other );
   Field< value_type, Stencil::Size > * cloneShallowCopyInternal() const override { return new PdfField< LatticeStorageSpecification_T >( *this ); }
   //@}
   //*******************************************************************************************************************

   LatticeStorageSpecification_T storageSpecification_;
};



template< typename LatticeStorageSpecification_T >
PdfField< LatticeStorageSpecification_T >::PdfField( const uint_t _xSize, const uint_t _ySize, const uint_t _zSize,
                                                    const LatticeStorageSpecification_T & storageSpecification,
                                      const uint_t ghostLayers, const field::Layout & _layout,
                                      const shared_ptr< field::FieldAllocator<value_type> > & alloc ) :

   GhostLayerField< value_type, Stencil::Size >( _xSize, _ySize, _zSize, ghostLayers, _layout, alloc ),
      storageSpecification_( storageSpecification )

{
#ifdef _OPENMP
   // take care of proper thread<->memory assignment (first-touch allocation policy !)
   this->setWithGhostLayer( value_type(0) );
#endif
   this->setWithGhostLayer( value_type(0) );
}



template< typename LatticeStorageSpecification_T >
inline PdfField< LatticeStorageSpecification_T > * PdfField< LatticeStorageSpecification_T >::clone() const
{
   return dynamic_cast< PdfField * >( GhostLayerField< value_type, Stencil::Size >::clone() );
}

template< typename LatticeStorageSpecification_T >
inline PdfField< LatticeStorageSpecification_T > * PdfField< LatticeStorageSpecification_T >::cloneUninitialized() const
{
   return dynamic_cast< PdfField * >( GhostLayerField< value_type, Stencil::Size >::cloneUninitialized() );
}

template< typename LatticeStorageSpecification_T >
inline PdfField< LatticeStorageSpecification_T > * PdfField< LatticeStorageSpecification_T >::cloneShallowCopy() const
{
   return dynamic_cast< PdfField * >( GhostLayerField< value_type, Stencil::Size >::cloneShallowCopy() );
}

template< typename LatticeStorageSpecification_T >
inline PdfField< LatticeStorageSpecification_T >::PdfField( const PdfField< LatticeStorageSpecification_T > & other )
   : GhostLayerField< value_type, Stencil::Size >::GhostLayerField( other )
{
}

} // namespace lbm
