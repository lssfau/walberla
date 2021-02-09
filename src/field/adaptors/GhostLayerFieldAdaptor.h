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
//! \file GhostLayerFieldAdaptor.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "AdaptorIterator.h"
#include "field/GhostLayerField.h"


namespace walberla {
namespace field {


//**********************************************************************************************************************
/*! Adaptor for walberla::field::GhostLayerField
*
* \ingroup field
*
* Motivation:
*  - assume you have a LBM PDF field and you want to print / view / process the velocities which can
*    be (easily) computed from the values of the PDF field
*  - one could create a second field, filled with the velocities, that is updated every time the pdfs changed
*  - the storage saving solution would be to create an adaptor that behaves just as a field, but holds no data.
*    Instead the data is computed when it is requested.
*
* Features:
*  - The GhostLayerFieldAdaptor implements (almost) the same members as a GhostLayerField
*     -> it can be passed to any templated function that expects a GhostLayerField
*  - The GhostLayerFieldAdaptor wraps one base-field, which in the example above would have been the pdf field
*  - the adaptor has automatically the same x,y,z size as the base-field.
*  - the adaptor can have a different fSize and a smaller or equal number of ghost layers
*
* Usage
*  - The customization of the adaptor (i.e. what it should compute) is done using a AdaptorFunctor that
*    has to implement the following concept:
*    \code
         class AdaptorFunctor
         {
            // required: type of the basis-field
            typedef Field_T  basefield_t;

            // required: f-size of the adaptor, can be different than basefield_t::F_SIZE
            static const uint_t F_SIZE = 1 ;

            // required: function that takes an iterator of the base_field and returns a transformed value
            typedef typename Field_T::value_type value_type;
            value_type operator() ( const basefield_iterator & it ) const {
               // for example: return lbm::calcRho( it );
            }

            // required: same as above, but taking coordinates instead of iterator
            value_type operator() ( const basefield_t & baseField,
                                    cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f = 0 ) const
            {
               // can be implemented like this: (but more efficient version may be possible)
               return (*this) ( basefield_t::Ptr( baseField, x,y,z,f ) );
            }
         };
*    \endcode
*
*/
//**********************************************************************************************************************
template< typename Functor, uint_t glDecrease >
class GhostLayerFieldAdaptor
{
public:
   //** Type Definitions  **********************************************************************************************
   /*! \name Type Definitions */
   //@{
   typedef Functor                       functor_t;
   typedef typename Functor::basefield_t basefield_t;
   typedef typename Functor::value_type  value_type;
   typedef typename Functor::value_type  T;

   static const uint_t F_SIZE = Functor::F_SIZE;

   typedef typename basefield_t::const_base_iterator     adapted_base_iterator;
   typedef typename basefield_t::const_iterator          adapted_iterator;
   typedef typename basefield_t::const_reverse_iterator  adapted_reverse_iterator;

   typedef GhostLayerFieldAdaptor<Functor,glDecrease> OwnType;
   typedef AdaptorIterator<adapted_base_iterator,   OwnType > const_base_iterator;
   typedef AdaptorIterator<adapted_iterator,        OwnType > const_iterator;
   typedef AdaptorIterator<adapted_reverse_iterator,OwnType > const_reverse_iterator;

   typedef FieldPointer<OwnType, const OwnType, const T > ConstPtr;
   //@}
   //*******************************************************************************************************************


   GhostLayerFieldAdaptor( const basefield_t & field, const Functor & func = Functor() )
      : functor_(func), glField_ ( field )
   {}

   // since Field/GhostLayerField is polymorphic, this class also has to be polymorphic
   // example: dynamic_cast< field_t > would fail if field_t = GhostLayerFieldAdaptor
   virtual ~GhostLayerFieldAdaptor() = default;


   const functor_t   & getFunctor()     const { return functor_; }
   const basefield_t & getAdaptedField()const { return glField_; }



   //** Equality Checks ************************************************************************************************
   /*! \name Equality Checks */
   //@{
   inline bool operator==      ( const OwnType & other ) const { return glField_ == other.glField_ ;                 }
   inline bool hasSameAllocSize( const OwnType & other ) const { return glField_.hasSameAllocSize( other.glField_ ); }
   inline bool hasSameSize     ( const OwnType & other) const  { return glField_.hasSameSize( other.glField_ );      }
   //@}
   //*******************************************************************************************************************



   //** Iterators  *****************************************************************************************************
   /*! \name Iterators */
   //@{
   //inline const_iterator begin() const;
   inline const_iterator beginXYZ() const;
   inline const_iterator beginSlice( cell_idx_t xBeg, cell_idx_t yBeg, cell_idx_t zBeg, cell_idx_t fBeg,
                                     cell_idx_t xEnd, cell_idx_t yEnd, cell_idx_t zEnd, cell_idx_t fEnd ) const;
   inline const_iterator beginSliceXYZ ( const CellInterval & interval, cell_idx_t f = 0 ) const;

   inline const const_iterator  end() const;

   inline const_iterator beginWithGhostLayer() const;
   inline const_iterator beginWithGhostLayerXYZ() const;

   inline const_iterator beginGhostLayerOnly( stencil::Direction dir ) const;
   inline const_iterator beginGhostLayerOnlyXYZ( stencil::Direction dir, cell_idx_t f = 0 ) const;
   inline const_iterator beginSliceBeforeGhostLayer( stencil::Direction dir, cell_idx_t width = 1 ) const;
   inline const_iterator beginSliceBeforeGhostLayerXYZ( stencil::Direction dir, cell_idx_t width = 1, cell_idx_t f = 0 ) const;

   inline void getGhostRegion( stencil::Direction dir, CellInterval & ci ) const;
   inline void getSliceBeforeGhostLayer( stencil::Direction d, CellInterval & ci, cell_idx_t width=1 ) const;

   //@}
   //*******************************************************************************************************************


   //** Adaptors for Field        **************************************************************************************
   /*! \name Adaptors for Field  */
   //@{
   inline uint_t  xSize() const  { return glField_.xSize(); }
   inline uint_t  ySize() const  { return glField_.ySize(); }
   inline uint_t  zSize() const  { return glField_.zSize(); }
   inline uint_t  fSize() const  { return Functor::F_SIZE;  }

   inline uint_t  xAllocSize() const  { return glField_.xAllocSize(); }
   inline uint_t  yAllocSize() const  { return glField_.yAllocSize(); }
   inline uint_t  zAllocSize() const  { return glField_.zAllocSize(); }
   inline uint_t  fAllocSize() const  { return Functor::F_SIZE;       }
   inline uint_t  allocSize()  const  { return glField_.allocSize();  }

   inline CellInterval xyzSize()      const { return glField_.xyzSize();       }
   inline CellInterval xyzAllocSize() const { return glField_.xyzAllocSize();  }

   bool coordinatesValid( cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f ) const;

   inline Layout layout() const { return glField_.layout(); }

   inline T operator()( cell_idx_t x, cell_idx_t y, cell_idx_t z) const               { return functor_( glField_, x,y,z);   }
   inline T operator()( cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f) const { return functor_( glField_, x,y,z,f); }
   inline T operator()( cell_idx_t x, cell_idx_t y, cell_idx_t z, uint_t f) const     { return functor_( glField_, x,y,z,cell_idx_c(f)); }
   inline T operator()( const Cell & c ) const                                        { return functor_( glField_, c[0], c[1], c[2] ); }
   inline T operator()( const Cell & c, cell_idx_t f ) const                          { return functor_( glField_, c[0], c[1], c[2], f ); }
   inline T operator()( const Cell & c, uint_t f ) const                              { return functor_( glField_, c[0], c[1], c[2], cell_idx_c(f) ); }


   inline T get( cell_idx_t x, cell_idx_t y, cell_idx_t z) const                      { return functor_( glField_, x, y, z    ); }
   inline T get( cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f ) const       { return functor_( glField_, x, y, z, f ); }
   inline T get( cell_idx_t x, cell_idx_t y, cell_idx_t z, uint_t f ) const           { return functor_( glField_, x, y, z, cell_idx_c(f) ); }
   inline T get( const Cell & c ) const                                               { return functor_( glField_, c[0], c[1], c[2] ); }
   inline T get( const Cell & c, cell_idx_t f ) const                                 { return functor_( glField_, c[0], c[1], c[2], f ); }
   inline T get( const Cell & c, uint_t f ) const                                     { return functor_( glField_, c[0], c[1], c[2], cell_idx_c(f) ); }
   //@}
   //*******************************************************************************************************************


   //** Adaptors for GhostLayerField        ************************************************************************
   /*! \name Adaptors for GhostLayerField  */
   //@{
   inline uint_t       xSizeWithGhostLayer()   const  { return glField_.xSizeWithGhostLayer() - 2 * glDecrease;   }
   inline uint_t       ySizeWithGhostLayer()   const  { return glField_.ySizeWithGhostLayer() - 2 * glDecrease;   }
   inline uint_t       zSizeWithGhostLayer()   const  { return glField_.zSizeWithGhostLayer() - 2 * glDecrease;   }
   inline uint_t       nrOfGhostLayers()       const  { return glField_.nrOfGhostLayers() - glDecrease;           }

   inline CellInterval xyzSizeWithGhostLayer() const;

   //@}
   //*******************************************************************************************************************

protected:
   const Functor functor_;
   const basefield_t & glField_;
};










//======================================================================================================================
//
//  Implementation
//
//======================================================================================================================

/*template< typename Functor, uint_t glD >
typename GhostLayerFieldAdaptor<Functor,glD>::const_iterator GhostLayerFieldAdaptor<Functor,glD>::begin() const
{
   return const_iterator( glField_.beginSlice(0,0,0,0,xSize(),ySize(),zSize(),F_SIZE), this );
}*/

template< typename Functor, uint_t glD >
typename GhostLayerFieldAdaptor<Functor,glD>::const_iterator GhostLayerFieldAdaptor<Functor,glD>::beginXYZ() const
{
   return const_iterator( glField_.beginXYZ(), this );
}

template< typename Functor, uint_t glD >
typename GhostLayerFieldAdaptor<Functor,glD>::const_iterator GhostLayerFieldAdaptor<Functor,glD>::beginSlice(
                           cell_idx_t xBeg, cell_idx_t yBeg, cell_idx_t zBeg, cell_idx_t fBeg,
                           cell_idx_t xEnd, cell_idx_t yEnd, cell_idx_t zEnd, cell_idx_t fEnd ) const
{
   WALBERLA_ASSERT_LESS_EQUAL( fEnd, cell_idx_c(F_SIZE) );
   return const_iterator( glField_.beginSlice(xBeg,yBeg,zBeg,fBeg, xEnd, yEnd, zEnd, fEnd ), this );
}

template< typename Functor, uint_t glD >
typename GhostLayerFieldAdaptor<Functor,glD>::const_iterator GhostLayerFieldAdaptor<Functor,glD>::beginSliceXYZ (
         const CellInterval & interval, cell_idx_t f ) const
{
   WALBERLA_ASSERT_LESS_EQUAL( f, cell_idx_c(F_SIZE) );
   return const_iterator( glField_.beginSliceXYZ(interval,f), this );
}

template< typename Functor, uint_t glD >
const typename  GhostLayerFieldAdaptor<Functor,glD>::const_iterator  GhostLayerFieldAdaptor<Functor,glD>::end() const {
   static const const_iterator cached_end ( glField_.end(), this );
   return cached_end;
}

template< typename Functor, uint_t glD >
typename GhostLayerFieldAdaptor<Functor,glD>::const_iterator GhostLayerFieldAdaptor<Functor,glD>::beginWithGhostLayer() const
{
   const uint_t gl = cell_idx_c( glField_.nrOfGhostLayers() - glD );
   return const_iterator( glField_.beginSlice( -gl, -gl,-gl,0,
                                               cell_idx_c( xSize() ) + gl,
                                               cell_idx_c( ySize() ) + gl,
                                               cell_idx_c( zSize() ) + gl,
                                               cell_idx_c( F_SIZE ) ),
                          this );
}

template< typename Functor, uint_t glD >
typename GhostLayerFieldAdaptor<Functor,glD>::const_iterator GhostLayerFieldAdaptor<Functor,glD>::beginWithGhostLayerXYZ() const
{
   const cell_idx_t gl = cell_idx_c( glField_.nrOfGhostLayers() - glD);
   return const_iterator( glField_.beginSlice( -gl, -gl,-gl,0,
                                               cell_idx_c( xSize() ) + gl,
                                               cell_idx_c( ySize() ) + gl,
                                               cell_idx_c( zSize() ) + gl,
                                               1),
                          this );
}

template< typename Functor, uint_t glD >
typename GhostLayerFieldAdaptor<Functor,glD>::const_iterator GhostLayerFieldAdaptor<Functor,glD>::beginGhostLayerOnly(
         stencil::Direction dir ) const
{
   CellInterval ci;
   getGhostRegion(dir,ci);

   return const_iterator( glField_.beginSlice( ci.xMin()  ,  ci.yMin(),   ci.zMin()  , 0,
                                               ci.xMax()+1,  ci.yMax()+1, ci.zMax()+1, F_SIZE ),
                          this);
}

template< typename Functor, uint_t glD >
typename GhostLayerFieldAdaptor<Functor,glD>::const_iterator GhostLayerFieldAdaptor<Functor,glD>::beginGhostLayerOnlyXYZ(
         stencil::Direction dir, cell_idx_t f ) const
{
   CellInterval ci;
   getGhostRegion(dir,ci);

   return const_iterator( glField_.beginSlice( ci.xMin()  ,  ci.yMin(),   ci.zMin()  , 0,
                                               ci.xMax()+1,  ci.yMax()+1, ci.zMax()+1, 1 ),
                           this);
}

template< typename Functor, uint_t glD >
typename GhostLayerFieldAdaptor<Functor,glD>::const_iterator GhostLayerFieldAdaptor<Functor,glD>::beginSliceBeforeGhostLayer(
         stencil::Direction dir, cell_idx_t width ) const
{
   CellInterval ci;
   getSliceBeforeGhostLayer(dir,ci,width);

   return const_iterator( glField_.beginSlice( ci.xMin()  ,  ci.yMin(),   ci.zMin()  , 0,
                                               ci.xMax()+1,  ci.yMax()+1, ci.zMax()+1, F_SIZE ),
                           this);
}

template< typename Functor, uint_t glD >
typename GhostLayerFieldAdaptor<Functor,glD>::const_iterator GhostLayerFieldAdaptor<Functor,glD>::beginSliceBeforeGhostLayerXYZ(
         stencil::Direction dir, cell_idx_t width, cell_idx_t f) const
{
   CellInterval ci;
   getSliceBeforeGhostLayer(dir,ci,width);

   return const_iterator( glField_.beginSlice( ci.xMin()  ,  ci.yMin(),   ci.zMin()  , 0,
                                               ci.xMax()+1,  ci.yMax()+1, ci.zMax()+1, 1 ),
                          this );
}


template< typename Functor, uint_t glD >
void GhostLayerFieldAdaptor<Functor,glD>::getGhostRegion( stencil::Direction dir, CellInterval & ci ) const
{
   glField_.getGhostRegion(dir,ci, cell_idx_c( glField_.nrOfGhostLayers() - glD ) );
}

template< typename Functor, uint_t glD >
void GhostLayerFieldAdaptor<Functor,glD>::getSliceBeforeGhostLayer( stencil::Direction d, CellInterval & ci, cell_idx_t width ) const
{
   glField_.getSliceBeforeGhostLayer(d, ci, width );
}

template< typename Functor, uint_t glD >
CellInterval GhostLayerFieldAdaptor<Functor,glD>::xyzSizeWithGhostLayer() const
{
   CellInterval ci = glField_.xyzSize();
   for(int i=0; i<3; ++i) {
      ci.min()[i] -= cell_idx_c( nrOfGhostLayers() );
      ci.max()[i] += cell_idx_c( nrOfGhostLayers() );
   }
   return ci;
}

template< typename Functor, uint_t glD >
bool GhostLayerFieldAdaptor<Functor,glD>::coordinatesValid( cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f ) const
{
   cell_idx_t gl = cell_idx_c( nrOfGhostLayers() - glD ) ;
   if ( x < -gl || y < -gl || z < -gl || f < 0 ) return false;

   if ( x >= cell_idx_c( xSize()) + gl ) return false;
   if ( y >= cell_idx_c( ySize()) + gl ) return false;
   if ( z >= cell_idx_c( zSize()) + gl ) return false;
   if ( f >= cell_idx_c( fSize() ) )     return false;
   return true;
}


} // namespace field
} // namespace walberla


