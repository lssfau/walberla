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
//! \file FlagField.impl.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Implementation of FlagField
//
//======================================================================================================================

#include <type_traits>

namespace walberla {
namespace field {


   //===================================================================================================================
   //
   //  CONSTRUCTORS
   //
   //===================================================================================================================

   //*******************************************************************************************************************
   /*!\brief Constructor. Creates a field without data.
    *
    * To use this field call GhostLayerField::init()
    *
    *******************************************************************************************************************/
   template<typename T>
   FlagField<T>::FlagField( )
   {
      data_ = new RegistrationData();
   }


   //*******************************************************************************************************************
   /*!\brief Constructor. Initializes field with 0 (no flags set).
    *
    * \param xSize  size of x dimension without ghost layers
    * \param ySize  size of y dimension without ghost layers
    * \param zSize  size of z dimension without ghost layers
    * \param gl     number of ghost layers
    * \param alloc  class that describes how to allocate memory for the field, see FieldAllocator
    *******************************************************************************************************************/
   template<typename T>
   FlagField<T>::FlagField( uint_t xs, uint_t ys, uint_t zs, uint_t gl, const shared_ptr<FieldAllocator<T> > &alloc)
      : GhostLayerField<T,1> (xs,ys,zs,gl,0, fzyx, alloc )
   {
      data_ = new RegistrationData();
   }


   //*******************************************************************************************************************
   /*!\brief Private copy constructor that creates a shallow copy of the data
    *******************************************************************************************************************/
    template<typename T>
    FlagField<T>::FlagField(const FlagField<T> & other)
       : GhostLayerField<T,1>::GhostLayerField(other),
         data_(other.data_)
    {
    }

    //******************************************************************************************************************
    /*!\brief Explicitly declared destructor, needed for handling the shallow copied Registration Data
     ******************************************************************************************************************/
   template<typename T>
   FlagField<T>::~FlagField()
   {
      uint_t refs = Field<T,1>::referenceCount();
      if( refs == 1 ) // last field that uses this data
         delete data_;
   }




   //===================================================================================================================
   //
   //  CLONING AND SLICING
   //
   //===================================================================================================================

   template<typename T>
   Field<T,1> * FlagField<T>::cloneShallowCopyInternal() const
   {
      return new FlagField<T>(*this);
   }

   template<typename T>
   inline FlagField<T> * FlagField<T>::clone() const
   {
      FlagField<T> * ff = dynamic_cast<FlagField<T>* > ( GhostLayerField<T,1>::clone() );
      // make a deep copy of Registration data, reference counting is done by FieldAllocator
      ff->data_ = new RegistrationData( *data_ );
      return ff;
   }

   template<typename T>
   inline FlagField<T> * FlagField<T>::cloneUninitialized() const
   {
      FlagField<T> * ff = dynamic_cast<FlagField<T>* > ( GhostLayerField<T,1>::cloneUninitialized() );
      // make a deep copy of Registration data, reference counting is done by FieldAllocator
      ff->data_ = new RegistrationData();
      return ff;
   }

   template<typename T>
   inline FlagField<T> * FlagField<T>::cloneShallowCopy() const
   {
      return dynamic_cast<FlagField<T>* > ( GhostLayerField<T,1>::cloneShallowCopy() );
   }

   template<typename T>
   FlagField<T> * FlagField<T>::getSlicedField( const CellInterval & ci ) const
   {
      return dynamic_cast<FlagField<T> * >( Field<T,1>::getSlicedField(ci) );
   }


   //===================================================================================================================
   //
   //  ACCESS FUNCTIONS
   //
   //===================================================================================================================


   //*******************************************************************************************************************
   /*!\brief Adds every (inner) cell where at least all bits of mask are set to the CellVector (in inner coordinates)
    *
    * \param mask [in]  bit mask. Test if a cell is added: (content & mask) == true
    * \param cv   [out] cell vector where the cells are added to, in inner coordinates
    *******************************************************************************************************************/
   template<typename T>
   void  FlagField<T>::getCellsWhereMaskIsSet(T mask, CellVector & cv) const
   {
      //check that mask contains only registered bits
      WALBERLA_ASSERT( ! ( mask & (~ data_->usedMask) ));

      for( auto i = this->begin(); i != this->end(); ++i)
         if(*i & mask )
            cv.push_back(Cell(i.x(),i.y(), i.z() ));
   }


   //*******************************************************************************************************************
   /*!\brief Equivalent to field::isMaskSet() with debug checks
   ********************************************************************************************************************/
   template<typename T>
   bool FlagField<T>::isMaskSet (cell_idx_t x, cell_idx_t y, cell_idx_t z, flag_t m) const
   {
      WALBERLA_ASSERT(isRegistered(m));
      return field::isMaskSet ( this->get(x,y,z), m );
   }

   //*******************************************************************************************************************
   /*!\brief Equivalent to field::isFlagSet() with debug checks
   ********************************************************************************************************************/
   template<typename T>
   bool FlagField<T>::isFlagSet (cell_idx_t x, cell_idx_t y, cell_idx_t z, flag_t f) const
   {
      WALBERLA_ASSERT(isRegistered(f));
      return field::isFlagSet( this->get(x,y,z), f );
   }


   //*******************************************************************************************************************
   /*!\brief Equivalent to field::isPartOfMaskSet() with debug checks
   ********************************************************************************************************************/
   template<typename T>
   bool FlagField<T>::isPartOfMaskSet (cell_idx_t x, cell_idx_t y, cell_idx_t z, flag_t m) const
   {
      WALBERLA_ASSERT(isRegistered(m));
      return field::isPartOfMaskSet( this->get(x,y,z), m );
   }


   //===================================================================================================================
   //
   //  FLAG REGISTRATION
   //
   //===================================================================================================================

   //*******************************************************************************************************************
   /*!\brief Registers a flag. Uses the next free bit in the field.
    *
    * \param uid    uid of the flag
    * \return       Mask, that has exactly one bit set to one.
    *******************************************************************************************************************/
   template<typename T>
   T FlagField<T>::registerFlag( const FlagUID & uid )
   {
      if( flagExists(uid) )
         throw std::runtime_error("Flag with name " + uid.getIdentifier() + " was already registered at FlagField.");

      // Skip the already registered bits, which where
      // registered by "void registerFlag(string,value)" variant
      while ( data_->nextFreeBit < sizeof(T)*8  && field::isFlagSet( data_->usedMask, T(T(1) << data_->nextFreeBit)) )
         data_->nextFreeBit++;

      if(data_->nextFreeBit >= sizeof(T)*8 )
         throw std::runtime_error( "Not enough space in flag_type for additional flags." );


      flag_t f = flag_t(T(1) << data_->nextFreeBit);
      data_->uidToFlag[uid] = f;

      data_->flagToUID[data_->nextFreeBit] = uid;
      ++data_->nextFreeBit;
      data_->usedMask = static_cast< flag_t >( data_->usedMask | f );

      return f;
   }



   //*******************************************************************************************************************
   /*!\brief Registers a flag and forces it to a given bitNr.
    *
    * If bitNr is not important use registerFlag(name) which assigns the next free bit to the flag.
    *
    * \param name   string identifier for the flag
    * \param bitNr  The bit nr associated with the flag (NOT the mask!). Has to be smaller than sizeof(flag_type)*8
    *               There exists also a function where this information is not needed, and the next free bit is used.
    * \return       Mask, that has the given bit set to one i.e. 1 << bitNr
    *******************************************************************************************************************/
   template<typename T>
   T FlagField<T>::registerFlag( const FlagUID & uid, uint_t bitNr )
   {
      if( bitNr >= sizeof(T)*8 )
         throw std::runtime_error("flag_t is too small to store the bit with given number.");

      if( flagExists(uid) )
         throw std::runtime_error("Flag with name " + uid.getIdentifier() + " was already registered at FlagField.");

      if ( flagExists(bitNr) )
         throw std::runtime_error("Already registered flag " + data_->flagToUID[bitNr].getIdentifier() +
                                    " at the given position of FlagField.");

      flag_t f = flag_t(T(1) << bitNr);
      data_->uidToFlag[uid] = f;
      data_->flagToUID[bitNr] = uid;
      data_->usedMask = static_cast< flag_t >( data_->usedMask | f );

      return f;
   }



   //*******************************************************************************************************************
   /*!\brief Returns the flag mask for a given UID.
    *
    * \param uid    UID for the flag
    * \return       Mask of this flag.
    *******************************************************************************************************************/
   template<typename T>
   T FlagField<T>::getFlag( const FlagUID & uid ) const
   {
      auto i = data_->uidToFlag.find(uid);
      if(i == data_->uidToFlag.end() )
         throw std::runtime_error("The flag named " + uid.getIdentifier() + " was not registered at the FlagField");

      return i->second;
   }


   //*******************************************************************************************************************
   /*!\brief Returns the flag mask for a given container of UIDs.
    *
    * \param uids   container of UIDs for the flag (has to provide begin() and end())
    * \return       Mask of this all flags.
    *******************************************************************************************************************/
   template<typename T>
   template< typename FlagUIDContainer >
   T FlagField<T>::getMask( const FlagUIDContainer & uids ) const
   {
      T mask(0);
      for( auto flagUID = uids.begin(); flagUID != uids.end(); ++flagUID )
         mask = static_cast< T >( mask | getFlag( *flagUID ) );
      return mask;
   }


   //*******************************************************************************************************************
   /*!\brief Returns the flag UID for a given mask.
   *
   * \param flag   mask for the flag (only one bit may be set)
   * \return       UID of this flag.
   ********************************************************************************************************************/
   template<typename T>
   const FlagUID & FlagField<T>::getFlagUID( const flag_t flag ) const
   {
      WALBERLA_ASSERT( math::uintIsPowerOfTwo( flag ) );

      uint_t bitNr = math::uintMSBPosition( flag ) - uint_t(1);

      WALBERLA_ASSERT_LESS( bitNr, data_->flagToUID.size() );

      if( !flagExists( bitNr ) )
      {
         std::ostringstream oss;
         oss << "The flag mask " << std::hex << flag << " was not registered at the FlagField";
         throw std::runtime_error( oss.str() );
      }

      return data_->flagToUID[ bitNr ];
   }


   //*******************************************************************************************************************
   /*!\brief Returns the flag mask for a given UID if already registered, otherwise registers flag
    *
    * \param uid    UID for the flag
    * \return       Mask of this flag.
    *******************************************************************************************************************/
   template<typename T>
   T FlagField<T>::getOrRegisterFlag( const FlagUID & uid )
   {
      if ( flagExists(uid) )
         return getFlag( uid );
      else
         return registerFlag( uid );
   }



   //*******************************************************************************************************************
   /*!\brief Test if flag was already registered.
    *
    * \param uid    UID for the flag
    * \return       boolean indicating if flag was already registered
    *******************************************************************************************************************/
   template<typename T>
   bool FlagField<T>::flagExists(const FlagUID & uid) const
   {
      return (data_->uidToFlag.find(uid) != data_->uidToFlag.end() );
   }

   //*******************************************************************************************************************
   /*!\brief Test if flag was already registered.
    *
    * \param bitNr  bit number (not mask!) of the flag
    * \return       boolean indicating if flag was already registered
    *******************************************************************************************************************/
   template<typename T>
   bool FlagField<T>::flagExists(uint_t bitNr) const
   {
      WALBERLA_ASSERT_LESS( bitNr, sizeof(T)*8 );
      return (data_->usedMask & (T(1)<<bitNr) ) > T(0);
   }


   //*******************************************************************************************************************
   /*!\brief Prints a list of "bitNr - flagName"
    *
    * \param os   output stream
    *******************************************************************************************************************/
   template<typename T>
   void FlagField<T>::printRegistered(std::ostream & os) const
   {
      for(uint_t i=0; i< sizeof(T)*8; ++i)
      {
         if( ( ( numeric_cast< flag_t >(1) << i ) & data_->usedMask ) != numeric_cast< flag_t >(0) )
            os << i << ":\t" << data_->flagToUID[i] << std::endl;
      }
   }


   //*******************************************************************************************************************
   /*!\brief Prints a list of set flags in the given cell
    *
    * \param os   output stream
    * \param cell all set flags of this cell are printed
    *******************************************************************************************************************/
   template<typename T>
   void FlagField<T>::printCell( std::ostream & os, const Cell & cell ) const
   {
      os << "Flags set in cell " << cell;
      for ( auto i = data_->uidToFlag.begin(); i != data_->uidToFlag.end(); ++i)
      {
          const std::string & name = i->first.getIdentifier();
          const flag_t & flagVal   = i->second;

          if( field::isFlagSet( this->get(cell), flagVal ) )
             os << "\n  - " << name;
      }
      os << "\n";
   }


   //*******************************************************************************************************************
   /*!\brief Appends the string representation of all registered flags at the given vector
    *******************************************************************************************************************/
   template<typename T>
   inline void FlagField<T>::getAllRegisteredFlags(std::vector<FlagUID> & out) const
   {
      for ( auto i = data_->uidToFlag.begin(); i != data_->uidToFlag.end(); ++i) {
         out.push_back( i->first );
      }
   }

   //*******************************************************************************************************************
    /*!\brief True, if all bits that are set in the mask have been registered at this FlagField
     ******************************************************************************************************************/
   template<typename T>
   inline bool FlagField<T>::isRegistered ( T mask) const
   {
      return ! ( mask & (~ data_->usedMask) ) ;
   }





    //==================================================================================================================
    //
    //  FREE ITERATOR FUNCTIONS
    //
    //==================================================================================================================

    // For documentation of the following functions see the equivalent functions in FlagFunctions.h

    template<typename T, typename FieldPtrOrIterator>
    inline void addMask( const FieldPtrOrIterator & it, T mask )
    {
       WALBERLA_ASSERT( dynamic_cast<const FlagField<T> * > ( it.getField() )->isRegistered(mask) );
       addMask( *it, mask );
    }

    template<typename T, typename FieldPtrOrIterator>
    inline void addFlag( const FieldPtrOrIterator & it, T flag )
    {
       WALBERLA_ASSERT( dynamic_cast<const FlagField<T> * > ( it.getField() )->isRegistered(flag) );
       addMask( *it, flag );
    }

    template<class T, typename FieldPtrOrIterator>
    inline void removeMask( const FieldPtrOrIterator & it, T mask )
    {
       WALBERLA_ASSERT( dynamic_cast<const FlagField<T> * > ( it.getField() )->isRegistered(mask) );
       removeMask( *it, mask );
    }

    template<class T, typename FieldPtrOrIterator>
    inline void removeFlag( const FieldPtrOrIterator & it, T flag )
    {
       WALBERLA_ASSERT( dynamic_cast<const FlagField<T> * > ( it.getField() )->isRegistered(flag) );
       removeMask( *it, flag );
    }

    template<class T, typename FieldPtrOrIterator>
    inline bool isMaskSet( const FieldPtrOrIterator & it, T mask )
    {
       WALBERLA_ASSERT( dynamic_cast<const FlagField<T> * > ( it.getField() )->isRegistered(mask) );
       return isMaskSet( *it, mask );
    }

    template<class T, typename FieldPtrOrIterator>
    inline bool isFlagSet( const FieldPtrOrIterator & it, T flag )
    {
       WALBERLA_ASSERT( dynamic_cast<const FlagField<T> * > ( it.getField() )->isRegistered(flag) );
       return isMaskSet( *it, flag );
    }

    template<class T, typename FieldPtrOrIterator>
    inline bool isPartOfMaskSet( const FieldPtrOrIterator & it, T mask )
    {
       WALBERLA_ASSERT( dynamic_cast<const FlagField<T> * > ( it.getField() )->isRegistered(mask) );
       return isPartOfMaskSet( *it, mask );
    }


    //******************************************************************************************************************
    /*!\brief Ores the neighborhood of the specified stencil and checks whether the bits of mask are set
     *
     * \param mask [in]  bit mask. Test if a cell is added: (content & mask) == true
     ******************************************************************************************************************/
    template<class Sten, typename FieldPtrOrIterator>
    inline bool isFlagInNeighborhood(const FieldPtrOrIterator & i, typename FieldPtrOrIterator::value_type mask)
    {
       using T = typename std::remove_const<typename FieldPtrOrIterator::value_type>::type;

       static_assert( (std::is_same< T,uint8_t >::value ||
                       std::is_same< T,uint16_t>::value ||
                       std::is_same< T,uint32_t>::value ||
                       std::is_same< T,uint64_t>::value),
                      "Only unsigned types of various lengths are allowed as type of FlagFields");

       T flag = 0;
       for( auto d = Sten::beginNoCenter(); d != Sten::end(); ++d ) {
          flag = static_cast<T>( flag | i.neighbor(*d) );
       }
       return (flag & mask) > T(0);
    }


    //******************************************************************************************************************
    /*!\brief Ores the neighborhood of the specified stencil and returns mask
     ******************************************************************************************************************/
    template<class Sten, typename FieldPtrOrIterator>
    inline typename std::remove_const<typename FieldPtrOrIterator::value_type>::type
       getOredNeighborhood(const FieldPtrOrIterator & i)
    {
       using RetType = typename std::remove_const<typename FieldPtrOrIterator::value_type>::type;

       RetType flag = 0;
       for( auto d = Sten::beginNoCenter(); d != Sten::end(); ++d ) {
          flag = static_cast<RetType>( flag | i.neighbor(*d) );
       }
       return flag;
    }



} // namespace field
} // namespace walberla




