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
//! \file Printers.impl.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Implementation of field printers
//
//======================================================================================================================



namespace walberla {
namespace field {


   template<typename T, uint_t fs>
   std::ostream &  printSlice( std::ostream & os, const Field<T,fs> & field, int sliceCoord,
                               cell_idx_t sliceValue, cell_idx_t f )
   {
      // If this field is a ghost-layer field call more specific variant of this function
      const GhostLayerField<T,fs> * glField = dynamic_cast< const GhostLayerField<T,fs> * > ( &field );
      if ( glField )
         return printSlice ( os, *glField, sliceCoord, sliceValue, f );


      using std::endl;
      os << endl;

      WALBERLA_ASSERT(sliceCoord >=0 && sliceCoord <3 );

      CellInterval size = field.xyzSize();

      Cell coord = size.min();
      const char * coordNames [3] = { "x", "y", "z" };

      int innerCoord = (sliceCoord + 1) % 3;
      int outerCoord = (sliceCoord + 2) % 3;

      size.min()[uint_c(sliceCoord)] = sliceValue;
      size.max()[uint_c(sliceCoord)] = sliceValue;

      auto sliced = shared_ptr<Field<T,fs> > ( field.getSlicedField( size ) );

      // Headline
      os << coordNames[innerCoord] << " =       ||";
      for (cell_idx_t i = coord[uint_c(innerCoord)] ; i <= size.max()[uint_c(innerCoord)]; ++i )
         os << std::setw(10) << i;
      os<<endl;

      auto precSave = os.precision();
      os.precision(3);

      // Separation line
      os << "================";
      for (cell_idx_t i = coord[uint_c(innerCoord)] ; i <= size.max()[uint_c(innerCoord)]; ++i )
         os << "==========";
      os << endl;

      for ( ; coord[uint_c(outerCoord)] <= size.max()[uint_c(outerCoord)]; ++coord[uint_c(outerCoord)] )
      {
         os << coordNames[uint_c(outerCoord)] << " = " << std::setw(5) << coord[uint_c(outerCoord)] << " ||\t";
         for (coord[uint_c(innerCoord)] = size.min()[uint_c(innerCoord)] ; coord[uint_c(innerCoord)] <= size.max()[uint_c(innerCoord)]; ++ coord[uint_c(innerCoord)] )
            os << std::setw(10) << sliced->get(coord[0],coord[1],coord[2],f);

         os << endl;
      }

      os << endl;

      os.precision(precSave);
      return os;
   }

   template<typename T, uint_t fs>
   std::ostream &  printSlice( std::ostream & os, const GhostLayerField<T,fs> & field, int sliceCoord,
                               cell_idx_t sliceValue, cell_idx_t f )
   {
      // If this field is a FlagField  call more specific variant of this function
      //const FlagField<T> * fField = dynamic_cast< const FlagField<T> * > ( &field );
      //if ( fField )
      //   return printSlice ( os, *fField, sliceCoord, sliceValue, f );

      using std::endl;
      os << endl;

      WALBERLA_ASSERT(sliceCoord >=0 && sliceCoord <3 );
      cell_idx_t glCellIdx = cell_idx_c ( field.nrOfGhostLayers() );

      const char * coordNames [3] = { "x", "y", "z" };

      int innerCoord = (sliceCoord + 1) % 3;
      int outerCoord = (sliceCoord + 2) % 3;

      CellInterval sliceInterval = field.xyzSize();

      sliceInterval.min()[uint_c(sliceCoord)] = sliceValue;
      sliceInterval.max()[uint_c(sliceCoord)] = sliceValue;

      auto sliced = shared_ptr<GhostLayerField<T,fs> > ( field.getSlicedField( sliceInterval ) );
      CellInterval size = field.xyzSizeWithGhostLayer();
      Cell coord = size.min();
      coord[uint_c(sliceCoord)] = 0;

      // Headline
      os << coordNames[uint_c(innerCoord)] << " =       ||    ";
      for (cell_idx_t i = coord[uint_c(innerCoord)] ; i <= size.max()[uint_c(innerCoord)]; ++i )
      {
         if ( i == 0 || i == size.max()[uint_c(innerCoord)]-glCellIdx+1)
            os << std::setw(3) << "|";

         os << std::setw(10) << i;
      }
      os<<endl;

      auto precSave = os.precision();
      os.precision(3);

      // Separation line
      os << "======================";
      for (cell_idx_t i = coord[uint_c(innerCoord)] ; i <= size.max()[uint_c(innerCoord)]; ++i )
         os << "==========";
      os << endl << endl;

      for ( ; coord[uint_c(outerCoord)] <= size.max()[uint_c(outerCoord)]; ++coord[uint_c(outerCoord)] )
      {

         if ( coord[uint_c(outerCoord)] == 0 || coord[uint_c(outerCoord)] == size.max()[uint_c(outerCoord)]-glCellIdx+1)
         {
            os << "----------------------";
            for (cell_idx_t i =  size.min()[uint_c(innerCoord)] ; i <= size.max()[uint_c(innerCoord)]; ++i )
               os << "----------";
            os << endl;
         }

         os << coordNames[uint_c(outerCoord)] << " = " << std::setw(5) << coord[uint_c(outerCoord)] << " ||\t";
         for ( coord[uint_c(innerCoord)] = size.min()[uint_c(innerCoord)] ; coord[uint_c(innerCoord)] <= size.max()[uint_c(innerCoord)]; ++ coord[uint_c(innerCoord)] ) {
            if ( coord[uint_c(innerCoord)] == 0 || coord[uint_c(innerCoord)] == size.max()[uint_c(innerCoord)]-glCellIdx+1 )
               os << std::setw(3) << "|";

            os << std::setw(10) << sliced->get(coord[0],coord[1],coord[2],f);

         }

         os << endl;
      }
      os << endl;
      os.precision(precSave);

      return os;
   }


   template<typename T>
   std::ostream &  printSlice( std::ostream & os, const FlagField<T> & field, int sliceCoord, cell_idx_t sliceValue, cell_idx_t )
   {
      using std::endl;

      WALBERLA_ASSERT(sliceCoord >=0 && sliceCoord <3 );
      cell_idx_t glCellIdx = cell_idx_c ( field.nrOfGhostLayers() );

      const char * coordNames [3] = { "x", "y", "z" };

      int innerCoord = (sliceCoord + 1) % 3;
      int outerCoord = (sliceCoord + 2) % 3;

      CellInterval sliceInterval = field.xyzSize();

      sliceInterval.min()[uint_c(sliceCoord)] = sliceValue;
      sliceInterval.max()[uint_c(sliceCoord)] = sliceValue;

      auto sliced = shared_ptr<FlagField<T> > ( field.getSlicedField( sliceInterval ) );
      CellInterval size = field.xyzSizeWithGhostLayer();
      Cell coord = size.min();
      coord[uint_c(sliceCoord)] = 0;

      // Print legend
      std::vector<field::FlagUID> allFlags;
      sliced->getAllRegisteredFlags( allFlags );
      os << endl << "Registered Flags: " << endl;
      for(auto i = allFlags.begin(); i != allFlags.end(); ++i )
      {
         os << "\t" << i->getIdentifier()[0] << " = " << *i << endl;
      }
      os << endl;


      // Headline
      os << coordNames[uint_c(innerCoord)] << " =       ||    ";
      for (cell_idx_t i = coord[uint_c(innerCoord)] ; i <= size.max()[uint_c(innerCoord)]; ++i )
      {
         if ( i == 0 || i == size.max()[uint_c(innerCoord)]-glCellIdx+1)
            os << std::setw(3) << "|";

         os << std::setw(10) << i;
      }
      os<<endl;

      // Separation line
      os << "======================";
      for (cell_idx_t i = coord[uint_c(innerCoord)] ; i <= size.max()[uint_c(innerCoord)]; ++i )
         os << "==========";
      os << endl << endl;

      for ( ; coord[uint_c(outerCoord)] <= size.max()[uint_c(outerCoord)]; ++coord[uint_c(outerCoord)] )
      {

         if ( coord[uint_c(outerCoord)] == 0 || coord[uint_c(outerCoord)] == size.max()[uint_c(outerCoord)]-glCellIdx+1)
         {
            os << "----------------------";
            for (cell_idx_t i =  size.min()[uint_c(innerCoord)] ; i <= size.max()[uint_c(innerCoord)]; ++i )
               os << "----------";
            os << endl;
         }

         os << coordNames[uint_c(outerCoord)] << " = " << std::setw(5) << coord[uint_c(outerCoord)] << " ||\t";
         for ( coord[uint_c(innerCoord)] = size.min()[uint_c(innerCoord)] ; coord[uint_c(innerCoord)] <= size.max()[uint_c(innerCoord)]; ++ coord[uint_c(innerCoord)] ) {
            if ( coord[uint_c(innerCoord)] == 0 || coord[uint_c(innerCoord)] == size.max()[uint_c(innerCoord)]-glCellIdx+1 )
               os << std::setw(3) << "|";

            std::string out;
            for( auto curFlag = allFlags.begin(); curFlag != allFlags.end(); ++curFlag )
               if ( sliced->isFlagSet(coord[0],coord[1],coord[2], sliced->getFlag(*curFlag) ) )
                  out.append( 1, curFlag->getIdentifier()[0] );

            os << std::setw(10) << out;
         }

         os << endl;
      }

      os << endl;
      return os;
   }


} // namespace field
} // namespace walberla
