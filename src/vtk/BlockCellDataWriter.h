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
//! \file BlockCellDataWriter.h
//! \ingroup vtk
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "Base64Writer.h"
#include "UtilityFunctions.h"

#include "core/DataTypes.h"
#include "core/debug/Debug.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include <ostream>
#include <string>


namespace walberla {
namespace vtk {



namespace internal {

//**********************************************************************************************************************
/*!
*   \brief Interface that is used by VTKOutput in order to write data to VTK files
*
*   If one wants to write block data to a VTK file, one must implement a class that derives from
*   vtk::internal::BlockCellDataWriter. An instance of this class then can be registered at a VTKOutput object. Every
*   time the output function of this VTKOutput object is triggered, all registered BlockCellDataWriter objects are
*   called.
*
*   Implementations of vtk::internal::BlockCellDataWriter are not supposed to inherit directly from
*   vtk::internal::BlockCellDataWriter but from vtk::BlockCellDataWriter which provides a much nicer interface for the
*   user to implement!
*
*   Every class derived from BlockCellDataWriter must implement two push functions that write one data element to either
*   an output stream or a Base64Writer. Typically, a private "evaluate" function is implemented and called in both push
*   functions:
*   \code
*      void push( std::ostream& os,  const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f )
*      {
*         vtk::toStream( os, evaluate(x,y,z,f) ); // to prevent the data output stream from writing ascii characters
*                                                 // for [unsigned] char data instead of numbers
*      }
*
*      void push( Base64Writer& b64, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f )
*      {
*         b64 << evaluate(x,y,z,f);
*      }
*
*      [DataType] evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f ) { ... }
*   \endcode
*   Both push functions always receive block local cell coordinates corresponding to the two protected members "block_"
*   and "blockStorage_".
*   The data elements which are pushed to the stream/Base64Writer must always be of the same type! Every data type is
*   assigned a specific string in VTK. This string must be returned by the pure virtual member function "typeString()".
*   This function should always be implemented by calling the utility function "vtk::typeToString()":
*   \code
*      std::string typeString() const { return vtk::typeToString< [DataType] >(); }
*   \endcode
*
*   Additionally, every class derived from BlockCellDataWriter can implement two further push functions that are used
*   for resampling the data. The parameters passed to these two functions are as follows:
*   - os/b64: the output stream/base64 writer
*   - x/y/z: block local cell coordinates of the cell that contains the point (globalX,globalY,globalZ)
*   - f: the "f value"/"dimension"/number of components
*   - local[X/Y/Z]Cell: coordinates of the point (globalX,globalY,globalZ) in block local cell coordinates
*                       In contrast to x/y/z (which are discrete values!), these coordinates are continuous/exact.
*   - global[X/Y/Z]: the center of the "resampled" cell in global coordinates (within the global 3D simulation space)
*   - samplingD[x/y/z]: the size/extent of the "resampled" cell in x-/y-/z-direction
*   If you do not provide your own implementation of these two functions, than the default implementation provided by
*   class BlockCellDataWriter is used which performs a simple nearest neighbor interpolation:
*   \code
*      void push( std::ostream& os, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f,
*                                   const real_t, const real_t, const real_t, const real_t, const real_t, const real_t,
*                                   const real_t, const real_t, const real_t ) { push( os, x, y, z, f ); }
*
*      void push( Base64Writer& b64, const cell_idx_t x, const cell_idx_t y,  const cell_idx_t z, const cell_idx_t f,
*                                    const real_t, const real_t, const real_t, const real_t, const real_t, const real_t,
*                                    const real_t, const real_t, const real_t ) { push( b64, x, y, z, f ); }
*   \endcode
*
*   Every time "configure( const IBlock& block, const StructuredBlockStorage& sbs )" is called, the block and its
*   corresponding structured block storage assigned to the BlockCellDataWriter may change. This triggers the call of a
*   protected member function "configure()" which must be implemented and is intended for reacting to these changes.
*
*/
//**********************************************************************************************************************

class BlockCellDataWriter
{
public:
            BlockCellDataWriter( const std::string& id ) : block_( nullptr ), blockStorage_( nullptr ), identifier_( id ) {}
   virtual ~BlockCellDataWriter() = default;

   void configure( const IBlock& block, const StructuredBlockStorage& sbs ) { block_ = &block; blockStorage_ = &sbs; configure(); }

   /// For the documentation of this function, please refer to the documentation/general description of this class.
   virtual void push( std::ostream& os,  const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f ) = 0;
   /// For the documentation of this function, please refer to the documentation/general description of this class.
   virtual void push( Base64Writer& b64, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f ) = 0;

   /// For the documentation of this function, please refer to the documentation/general description of this class.
   virtual inline void push( std::ostream& os,  const cell_idx_t x,      const cell_idx_t y,      const cell_idx_t z,      const cell_idx_t f,
                                                const real_t localXCell, const real_t localYCell, const real_t localZCell,
                                                const real_t globalX,    const real_t globalY,    const real_t globalZ,
                                                const real_t samplingDx, const real_t samplingDy, const real_t samplingDz );
   /// For the documentation of this function, please refer to the documentation/general description of this class.
   virtual inline void push( Base64Writer& b64, const cell_idx_t x,      const cell_idx_t y,      const cell_idx_t z,      const cell_idx_t f,
                                                const real_t localXCell, const real_t localYCell, const real_t localZCell,
                                                const real_t globalX,    const real_t globalY,    const real_t globalZ,
                                                const real_t samplingDx, const real_t samplingDy, const real_t samplingDz );

           uint_t xSize() const { WALBERLA_ASSERT_NOT_NULLPTR( block_ ); WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ ); return blockStorage_->getNumberOfXCells( *block_ ); }
           uint_t ySize() const { WALBERLA_ASSERT_NOT_NULLPTR( block_ ); WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ ); return blockStorage_->getNumberOfYCells( *block_ ); }
           uint_t zSize() const { WALBERLA_ASSERT_NOT_NULLPTR( block_ ); WALBERLA_ASSERT_NOT_NULLPTR( blockStorage_ ); return blockStorage_->getNumberOfZCells( *block_ ); }
   virtual uint_t fSize() const = 0; //!< must return the size of the fourth dimension
                                     /*!< (data fields storing scalar values return "1", data fields storing
                                      *    vector data return the size of the vector) */

   //*******************************************************************************************************************
   /*!
   *   Every data type is assigned a specific string in VTK. This string must be returned by this function, which should
   *   always be implemented by calling the utility function "vtk::typeToString()":
   *   \code
   *      std::string typeString() const { return vtk::typeToString< [DataType] >(); }
   *   \endcode
   */
   //*******************************************************************************************************************
   virtual std::string typeString() const = 0;

   const std::string & identifier() const { return identifier_; };

protected:

   //*******************************************************************************************************************
   /*!
   *   Every time "configure( const IBlock& block, const StructuredBlockStorage& sbs )" is called, the block and its
   *   corresponding structured block storage assigned to the BlockCellDataWriter may change. This triggers the call of
   *   this function which must be implemented and is intended for reacting to these changes.
   */
   //*******************************************************************************************************************
   virtual void configure() = 0;

   void setIdentifier( const std::string & id ) { this->identifier_ = id; }

   const IBlock* block_;
   const StructuredBlockStorage* blockStorage_;

   std::string identifier_;

private:

   BlockCellDataWriter() = default;

}; // class BlockCellDataWriter



inline void BlockCellDataWriter::push( std::ostream& os, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f,
                                       const real_t /*localXCell*/, const real_t /*localYCell*/, const real_t /*localZCell*/,
                                       const real_t /*globalX*/,    const real_t /*globalY*/,    const real_t /*globalZ*/,
                                       const real_t /*samplingDx*/, const real_t /*samplingDy*/, const real_t /*samplingDz*/ )
{
   push( os, x, y, z, f );
}



inline void BlockCellDataWriter::push( Base64Writer& b64, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f,
                                       const real_t /*localXCell*/, const real_t /*localYCell*/, const real_t /*localZCell*/,
                                       const real_t /*globalX*/,    const real_t /*globalY*/,    const real_t /*globalZ*/,
                                       const real_t /*samplingDx*/, const real_t /*samplingDy*/, const real_t /*samplingDz*/ )
{
   push( b64, x, y, z, f );
}

} // namespace internal

using BlockCellDataWriterInterface = internal::BlockCellDataWriter;



//**********************************************************************************************************************
/*!
*   \brief Interface that is used by VTKOutput in order to write data to VTK files
*
*   For a detailed documentation see the documentation of class vtk::internal::BlockCellDataWriter.
*
*   Classes deriving from BlockCellDataWriter must implement a protected member function "evaluate(x,y,z,f)" and
*   "configure()" [pure virtual function inherited from vtk::internal::BlockCellDataWriter]. They can also implement a
*   second evaluate function used during resampling (if this function is not implemented nearest neighbor interpolation
*   will be used for resampling the data).
*/
//**********************************************************************************************************************
template< typename T, uint_t F_SIZE_ARG = 1u >
class BlockCellDataWriter : public BlockCellDataWriterInterface
{

public:

   using value_type = T;

   static const uint_t F_SIZE = F_SIZE_ARG;

            BlockCellDataWriter( const std::string & id ) : BlockCellDataWriterInterface( id ) {}
   ~BlockCellDataWriter() override = default;

   void push( std::ostream & os, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f ) override
   {
      vtk::toStream( os, evaluate( x, y, z, f ) );
   }

   void push( vtk::Base64Writer & b64, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f ) override
   {
      b64 << evaluate( x, y, z, f );
   }

   void push( std::ostream& os,  const cell_idx_t x,      const cell_idx_t y,      const cell_idx_t z,      const cell_idx_t f,
                                 const real_t localXCell, const real_t localYCell, const real_t localZCell,
                                 const real_t globalX,    const real_t globalY,    const real_t globalZ,
                                 const real_t samplingDx, const real_t samplingDy, const real_t samplingDz ) override
   {
      vtk::toStream( os, evaluate( x, y, z, f, localXCell, localYCell, localZCell,
                                   globalX, globalY, globalZ, samplingDx, samplingDy, samplingDz ) );
   }

   void push( Base64Writer& b64, const cell_idx_t x,      const cell_idx_t y,      const cell_idx_t z,      const cell_idx_t f,
                                 const real_t localXCell, const real_t localYCell, const real_t localZCell,
                                 const real_t globalX,    const real_t globalY,    const real_t globalZ,
                                 const real_t samplingDx, const real_t samplingDy, const real_t samplingDz ) override
   {
      b64 << evaluate( x, y, z, f, localXCell, localYCell, localZCell, globalX, globalY, globalZ, samplingDx, samplingDy, samplingDz );
   }

   uint_t fSize() const override { return F_SIZE; }

   std::string typeString() const override { return vtk::typeToString< T >(); }

protected:
   virtual T evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f ) = 0;

   virtual T evaluate( const cell_idx_t x,      const cell_idx_t y,      const cell_idx_t z,      const cell_idx_t f,
                       const real_t localXCell, const real_t localYCell, const real_t localZCell,
                       const real_t globalX,    const real_t globalY,    const real_t globalZ,
                       const real_t samplingDx, const real_t samplingDy, const real_t samplingDz );

   // virtual void configure() = 0; // -> don't forget to implement this member function!

}; // class BlockCellDataWriter



template< typename T, uint_t F_SIZE_ARG >
inline T BlockCellDataWriter< T, F_SIZE_ARG >::evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f,
                                                         const real_t /*localXCell*/, const real_t /*localYCell*/, const real_t /*localZCell*/,
                                                         const real_t /*globalX*/,    const real_t /*globalY*/,    const real_t /*globalZ*/,
                                                         const real_t /*samplingDx*/, const real_t /*samplingDy*/, const real_t /*samplingDz*/ )
{
   return evaluate(x,y,z,f);
}


template< typename T >
class BlockCellDataWriterScalingAdapter : public T
{
public:
   using value_type = typename T::value_type;
   static const uint_t F_SIZE = T::F_SIZE;

   BlockCellDataWriterScalingAdapter( const std::string& id, const T & base, value_type factor ) 
      : T( base ), factor_( factor ) { this->setIdentifier( id ); }

protected:

   virtual value_type evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f )
   {
      return factor_ * T::evaluate( x, y, z, f );
   }

private:
   value_type factor_;
};

template< typename T >
shared_ptr< BlockCellDataWriterScalingAdapter<T> > makeBlockCellDataWriterScalingAdapter( 
   const std::string& id, const shared_ptr<T> & base, typename T::value_type factor )
{
   return walberla::make_shared< BlockCellDataWriterScalingAdapter<T> >( id, *base, factor );
}




} // namespace vtk
} // namespace walberla


