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
//! \file FlagField.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief GhostLayerField that stores flags.
//
//======================================================================================================================

#pragma once

#include "FlagFunctions.h"
#include "FlagUID.h"
#include "GhostLayerField.h"

#include "core/DataTypes.h"
#include "core/cell/CellVector.h"
#include "core/debug/Debug.h"
#include "core/math/Uint.h"

#include <array>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <type_traits>

namespace walberla {
namespace field {

//**********************************************************************************************************************
/*! GhostLayerField that stores bit masks.
 *
 * \ingroup field
 *
 * Two main functionalities:
 *    - Flag access functions: add/remove flags
 *    - Registration/Mapping of strings to bits and vice versa
 *
 *
 * Registration and UID Mechanism:
 *    - Each bit of the flag type is assigned a UID which is basically just a number representing a string
 *      UID's are used here to avoid string comparisons, and to do integer comparisons instead.
 *      The FlagUID class internally holds a uint_t and has a static map that maps the uint_t to a string
 *    - Each registration function exists also in a variant that takes a string. These functions generated the
 *      UID internally.
 *    - All registration data structures are held indirectly by a pointer to allow for a shallow copy of the FlagField
 *
 * Flag / Mask operations:
 *    - Difference Mask/Flag: Flag is a mask where exactly one bit is set. All functions that have "Flag" in their name,
 *      enforce this restriction in debug-mode. In release mode they are equivalent to their "Mask" counterparts.
 *    - Mask operations are function to add, remove and query for bits in a bit-mask.
 *      There are three groups of them:
 *       - There are free functions in FlagFunctions.h that take the mask/flag and the value to modify.
 *         Since they are free functions they cannot check, if the bits they operate one, have been registered.
 *         Intended use case are constructs like: addFlag( fieldIterator.neighbor(d), someFlag )
 *       - The second group takes an iterator, instead of a value. The iterator holds the corresponding flagField
 *         and therefore these functions can do more checks, and should be favored to "free" ones.
 *         Example:  addFlag ( fieldIterator, someFlag)
 *         Note that the iterator is not dereferenced here! When dereferenced the first group of functions would be used
 *       - The third group are the member functions of the FlagField, taking (x,y,z) coordinates. They
 *         of course can also do all checks that the second group can do.
 *
 * See also \ref fieldPage
 *
*/
//**********************************************************************************************************************
template<typename T = uint32_t>
class FlagField : public GhostLayerField<T,1>
{
public:

   //** Type Definitions  **********************************************************************************************
   /*! \name Type Definitions */
   //@{
   using flag_t = T;

   using value_type = typename GhostLayerField<T, 1>::value_type;

   using iterator = typename GhostLayerField<T, 1>::iterator;
   using const_iterator = typename GhostLayerField<T, 1>::const_iterator;

   using reverse_iterator = typename GhostLayerField<T, 1>::reverse_iterator;
   using const_reverse_iterator = typename GhostLayerField<T, 1>::const_reverse_iterator;

   using base_iterator = typename Field<T, 1>::base_iterator;
   using const_base_iterator = typename Field<T, 1>::const_base_iterator;

   using Ptr = typename GhostLayerField<T, 1>::Ptr;
   using ConstPtr = typename GhostLayerField<T, 1>::ConstPtr;
   //@}
   //*******************************************************************************************************************


   //**Construction & Destruction***************************************************************************************
   /*! \name Construction & Destruction */
   //@{

   FlagField( uint_t xSize, uint_t ySize, uint_t zSize, uint_t gl,
              const shared_ptr<FieldAllocator<T> > &alloc = make_shared<StdFieldAlloc<T> >());
   ~FlagField() override;

   inline FlagField<T> * clone()              const;
   inline FlagField<T> * cloneUninitialized() const;
   inline FlagField<T> * cloneShallowCopy()   const;
   //@}
   //*******************************************************************************************************************


   //** Access Functions ***********************************************************************************************
   /*! \name Access Functions */
   //@{
   using idx = cell_idx_t;

   void addMask    (idx x, idx y, idx z, flag_t m) { WALBERLA_ASSERT(isRegistered(m)); field::addMask( this->get(x,y,z), m ); }
   void addFlag    (idx x, idx y, idx z, flag_t f) { WALBERLA_ASSERT(isRegistered(f)); field::addFlag( this->get(x,y,z), f ); }

   void addMask    ( const Cell & cell, flag_t m ) { addMask( cell.x(), cell.y(), cell.z(), m ); }
   void addFlag    ( const Cell & cell, flag_t f ) { addFlag( cell.x(), cell.y(), cell.z(), f ); }

   void removeMask (idx x, idx y, idx z, flag_t m) { WALBERLA_ASSERT(isRegistered(m)); field::removeMask( this->get(x,y,z), m ); }
   void removeFlag (idx x, idx y, idx z, flag_t f) { WALBERLA_ASSERT(isRegistered(f)); field::removeFlag( this->get(x,y,z), f ); }

   void removeMask ( const Cell & cell, flag_t m ) { removeMask( cell.x(), cell.y(), cell.z(), m ); }
   void removeFlag ( const Cell & cell, flag_t f ) { removeFlag( cell.x(), cell.y(), cell.z(), f ); }

   bool isMaskSet       (idx x, idx y, idx z, flag_t m) const;
   bool isFlagSet       (idx x, idx y, idx z, flag_t f) const;
   bool isPartOfMaskSet (idx x, idx y, idx z, flag_t m) const;

   bool isMaskSet       ( const Cell & cell, flag_t m ) const { return isMaskSet      ( cell.x(), cell.y(), cell.z(), m ); }
   bool isFlagSet       ( const Cell & cell, flag_t f ) const { return isFlagSet      ( cell.x(), cell.y(), cell.z(), f ); }
   bool isPartOfMaskSet ( const Cell & cell, flag_t m ) const { return isPartOfMaskSet( cell.x(), cell.y(), cell.z(), m ); }

   inline void getCellsWhereMaskIsSet( flag_t mask, CellVector & out ) const;
   //@}
   //*******************************************************************************************************************



   //** Flag Registration **********************************************************************************************
   /*! \name Flag Registration */
   //@{

   inline flag_t registerFlag( const FlagUID & uid );
   inline flag_t registerFlag( const FlagUID & uid, uint_t bitNr );

   inline flag_t          getFlag( const FlagUID & uid ) const;
   inline const FlagUID & getFlagUID( const flag_t flag ) const;
   inline flag_t          getOrRegisterFlag( const FlagUID & uid );

   template< typename FlagUIDContainer >
   inline flag_t getMask( const FlagUIDContainer & uids ) const;

   inline bool   flagExists( const FlagUID & uid ) const;
   inline bool   flagExists( uint_t bitNr )        const;

   inline void   printRegistered( std::ostream & os ) const;
   inline void   printCell( std::ostream & os, const Cell & cell ) const;

   inline bool   isRegistered ( flag_t mask )  const;

   inline void   getAllRegisteredFlags( std::vector<FlagUID> & out ) const;

   inline const std::map<FlagUID, flag_t> & getMapping() const { WALBERLA_ASSERT_NOT_NULLPTR( data_ ); return data_->uidToFlag; }
   //@}
   //*******************************************************************************************************************


   //** Slicing  *******************************************************************************************************
   /*! \name Slicing */
   //@{
   FlagField<T> * getSlicedField( const CellInterval & interval ) const;
   //@}
   //*******************************************************************************************************************


protected:
   FlagField();


   //** Shallow Copy ***************************************************************************************************
   /*! \name Shallow Copy */
   //@{
   FlagField(const FlagField<T> & other);
   Field<T,1> * cloneShallowCopyInternal()   const override;
   //@}
   //*******************************************************************************************************************

   // All Members are hold inside an extra struct, to enable shallow copies of the field
   struct RegistrationData
   {
      RegistrationData() : usedMask(0), nextFreeBit(0) {}
      RegistrationData( const RegistrationData & o )
         : flagToUID  ( o.flagToUID   ),
           uidToFlag  ( o.uidToFlag   ),
           usedMask   ( o.usedMask    ),
           nextFreeBit( o.nextFreeBit )
      {}

      /// Maps bitNr's to strings
      std::array<FlagUID, sizeof(flag_t)*8> flagToUID;

      /// Maps strings to bit-masks
      std::map<FlagUID, flag_t> uidToFlag;

      /// Mask is one for every registered bit
      flag_t usedMask;

      /// Counter for registerFlag(name)
      /// BitNumbers smaller than nextFreeBit are guaranteed to be
      /// occupied. Bits greater or equal may be occupied, when registerFlac(name,bitNr)
      /// was called.
      uint_t nextFreeBit;
   };
   RegistrationData * data_;



   static_assert( (std::is_same<T,uint8_t >::value ||
                   std::is_same<T,uint16_t>::value ||
                   std::is_same<T,uint32_t>::value ||
                   std::is_same<T,uint64_t>::value),
                  "Only unsigned types of various lengths are allowed as type of FlagFields");


   //** Free Iterator Access functions friend *************************************************************************
   /*! \name Free Iterator Access functions friend */
   //@{
   template <class FT, typename FieldPtrOrIterator > friend void addMask        ( const FieldPtrOrIterator & it, FT mask );
   template <class FT, typename FieldPtrOrIterator > friend void addFlag        ( const FieldPtrOrIterator & it, FT flag );
   template <class FT, typename FieldPtrOrIterator > friend void removeMask     ( const FieldPtrOrIterator & it, FT mask );
   template <class FT, typename FieldPtrOrIterator > friend void removeFlag     ( const FieldPtrOrIterator & it, FT flag );
   template <class FT, typename FieldPtrOrIterator > friend bool isMaskSet      ( const FieldPtrOrIterator & it, FT mask );
   template <class FT, typename FieldPtrOrIterator > friend bool isFlagSet      ( const FieldPtrOrIterator & it, FT flag );
   template <class FT, typename FieldPtrOrIterator > friend bool isPartOfMaskSet( const FieldPtrOrIterator & it, FT mask );

   template <class Stencil, typename FieldPtrOrIterator>
   friend bool isFlagInNeighborhood( const FieldPtrOrIterator & i, typename FieldPtrOrIterator::value_t mask);

   template <class Stencil, typename FieldPtrOrIterator>
   friend typename std::remove_const<typename FieldPtrOrIterator::value_type>::type
      getOredNeighborhood(const FieldPtrOrIterator & i);

   //@}
   //*******************************************************************************************************************

};


} // namespace field
} // namespace walberla

#include "FlagField.impl.h"


//======================================================================================================================
//
//  EXPORTS
//
//======================================================================================================================

namespace walberla {
   // Export flag field class from walberla::field to walberla namespace
   using field::FlagField;

   // Free functions that take iterator as first argument
   using field::addMask;
   using field::removeMask;
   using field::isMaskSet;
   using field::isPartOfMaskSet;

   using field::addFlag;
   using field::removeFlag;
   using field::isFlagSet;

   using field::isFlagInNeighborhood;
   using field::getOredNeighborhood;
}
