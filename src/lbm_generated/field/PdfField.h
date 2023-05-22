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
class PdfField : public GhostLayerField< real_t, LatticeStorageSpecification_T::Stencil::Size >
{
public:

   //** Type Definitions  **********************************************************************************************
   /*! \name Type Definitions */
   //@{
   using LatticeStorageSpecification = LatticeStorageSpecification_T;
   using Stencil = typename LatticeStorageSpecification_T::Stencil;

   using value_type = typename GhostLayerField<real_t, Stencil::Size>::value_type;

   using Ptr = typename GhostLayerField<real_t, Stencil::Size>::Ptr;
   using ConstPtr = typename GhostLayerField<real_t, Stencil::Size>::ConstPtr;
   //@}
   //*******************************************************************************************************************

   PdfField( const uint_t _xSize, const uint_t _ySize, const uint_t _zSize,
            const LatticeStorageSpecification_T & storageSpecification,
             const uint_t ghostLayers = uint_t(1), const field::Layout & _layout = field::zyxf,
             const shared_ptr< field::FieldAllocator<real_t> > & alloc = shared_ptr< field::FieldAllocator<real_t> >() );

   ~PdfField() override = default;

   inline PdfField * clone()              const;
   inline PdfField * cloneUninitialized() const;
   inline PdfField * cloneShallowCopy()   const;


   /////////////////////////////////////////////////
   // Access functions (with stencil::Direction!) //
   /////////////////////////////////////////////////

   using GhostLayerField< real_t, Stencil::Size >::get;

         real_t & get( cell_idx_t x, cell_idx_t y, cell_idx_t z, stencil::Direction d )       { return get( x, y, z, Stencil::idx[d] ); }
   const real_t & get( cell_idx_t x, cell_idx_t y, cell_idx_t z, stencil::Direction d ) const { return get( x, y, z, Stencil::idx[d] ); }
         real_t & get( const Cell & c, stencil::Direction d )       { return get( c.x(), c.y(), c.z(), Stencil::idx[d] ); }
   const real_t & get( const Cell & c, stencil::Direction d ) const { return get( c.x(), c.y(), c.z(), Stencil::idx[d] ); }

   using GhostLayerField< real_t, Stencil::Size >::operator();

         real_t & operator()( cell_idx_t x, cell_idx_t y, cell_idx_t z, stencil::Direction d )       { return get( x, y, z, Stencil::idx[d] ); }
   const real_t & operator()( cell_idx_t x, cell_idx_t y, cell_idx_t z, stencil::Direction d ) const { return get( x, y, z, Stencil::idx[d] ); }
         real_t & operator()( const Cell & c, stencil::Direction d )       { return get( c.x(), c.y(), c.z(), Stencil::idx[d] ); }
   const real_t & operator()( const Cell & c, stencil::Direction d ) const { return get( c.x(), c.y(), c.z(), Stencil::idx[d] ); }


protected:
   //** Shallow Copy ***************************************************************************************************
   /*! \name Shallow Copy */
   //@{
   inline PdfField( const PdfField< LatticeStorageSpecification_T > & other );
   Field< real_t, Stencil::Size > * cloneShallowCopyInternal() const override { return new PdfField< LatticeStorageSpecification_T >( *this ); }
   //@}
   //*******************************************************************************************************************

   LatticeStorageSpecification_T storageSpecification_;
};



template< typename LatticeStorageSpecification_T >
PdfField< LatticeStorageSpecification_T >::PdfField( const uint_t _xSize, const uint_t _ySize, const uint_t _zSize,
                                                    const LatticeStorageSpecification_T & storageSpecification,
                                      const uint_t ghostLayers, const field::Layout & _layout,
                                      const shared_ptr< field::FieldAllocator<real_t> > & alloc ) :

   GhostLayerField< real_t, Stencil::Size >( _xSize, _ySize, _zSize, ghostLayers, _layout, alloc ),
      storageSpecification_( storageSpecification )

{
#ifdef _OPENMP
   // take care of proper thread<->memory assignment (first-touch allocation policy !)
   this->setWithGhostLayer( real_t(0) );
#endif
   this->setWithGhostLayer( real_t(0) );
}



template< typename LatticeStorageSpecification_T >
inline PdfField< LatticeStorageSpecification_T > * PdfField< LatticeStorageSpecification_T >::clone() const
{
   return dynamic_cast< PdfField * >( GhostLayerField< real_t, Stencil::Size >::clone() );
}

template< typename LatticeStorageSpecification_T >
inline PdfField< LatticeStorageSpecification_T > * PdfField< LatticeStorageSpecification_T >::cloneUninitialized() const
{
   return dynamic_cast< PdfField * >( GhostLayerField< real_t, Stencil::Size >::cloneUninitialized() );
}

template< typename LatticeStorageSpecification_T >
inline PdfField< LatticeStorageSpecification_T > * PdfField< LatticeStorageSpecification_T >::cloneShallowCopy() const
{
   return dynamic_cast< PdfField * >( GhostLayerField< real_t, Stencil::Size >::cloneShallowCopy() );
}

template< typename LatticeStorageSpecification_T >
inline PdfField< LatticeStorageSpecification_T >::PdfField( const PdfField< LatticeStorageSpecification_T > & other )
   : GhostLayerField< real_t, Stencil::Size >::GhostLayerField( other )
{
}

} // namespace lbm
