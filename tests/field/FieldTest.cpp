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
//! \file FieldTest.cpp
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Tests basic functionality of Field and GhostLayerField
//
//======================================================================================================================

#include "field/Field.h"
#include "field/GhostLayerField.h"
#include "field/Printers.h"
#include "field/SwapableCompare.h"
#include "field/iterators/FieldPointer.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "stencil/D3Q27.h"

#include <iostream>
#include <set>


using namespace walberla;

using std::cout;
using std::endl;


void alignedAllocTest()
{
   void * p;

   p = field::aligned_malloc(5,16);
   WALBERLA_CHECK_EQUAL( ((size_t)p) % 16, 0u  );
   field::aligned_free(p);

   p = field::aligned_malloc(5,64);
   WALBERLA_CHECK_EQUAL( ((size_t)p) % 64, 0u  );
   field::aligned_free(p);
}

void alignedAllocWithOffsetTest()
{
   char * p;
   p = (char*)field::aligned_malloc_with_offset(64 + 2*8,32,8);
   WALBERLA_CHECK_EQUAL( ((size_t)p+8) % 16, 0u  );
   field::aligned_free(p);
}

void simpleCreateAndIterate(field::Layout layout)
{
   const uint_t xs = 3;
   const uint_t ys = 4;
   const uint_t zs = 2;
   const uint_t fs = 2;

   Field<uint_t,fs> field(xs,ys,zs,layout);

   uint_t counter=0;
   for(cell_idx_t f=0; f < cell_idx_c( field.fSize() ); ++f )
      for(cell_idx_t z=0; z < cell_idx_t( field.zSize() ); ++z )
         for(cell_idx_t y=0; y < cell_idx_t( field.ySize() ); ++y )
            for(cell_idx_t x=0; x < cell_idx_t( field.xSize( )); ++x ) {
               uint_t val = uint_c( x*y*z*f + x + y + z + f );
               field(x,y,z,f) = val;
               WALBERLA_CHECK_EQUAL( field(x,y,z,f),     val);
               WALBERLA_CHECK_EQUAL( field.get(x,y,z,f), val);
               counter++;
            }
   WALBERLA_CHECK_EQUAL(counter, xs*ys*zs*fs);


   counter = 0;
   for( auto i = field.begin(); i!= field.end(); ++i )
   {
      //cout << counter << "\t(" <<i.f() << "," << i.z() << "," << i.y() << "," << i.x() <<") Value ="<< *i << endl;
      WALBERLA_CHECK_EQUAL( field(i.x(),i.y(),i.z(),i.f()) , *i );
      counter ++;
   }
   WALBERLA_CHECK_EQUAL(counter, xs*ys*zs*fs);





}

void alignmentTest()
{
   const uint_t xs = 3;
   const uint_t ys = 4;
   const uint_t zs = 2;
   const uint_t fs = 2;

   const uint_t alignment =64;
   typedef field::AllocateAligned<short,alignment> Allocator;
   shared_ptr<Allocator> alloc = make_shared<Allocator>();
   Field<short,fs > field (xs,ys,zs,field::fzyx,alloc);
   for(cell_idx_t f=0; f < cell_idx_c( field.fSize() ); ++f )
      for(cell_idx_t z=0; z < cell_idx_t( field.zSize() ); ++z )
         for(cell_idx_t y=0; y < cell_idx_t( field.ySize() ); ++y )
         {
            void * p = reinterpret_cast<void*>( &(field(0,y,z,f)) );
            // test that each line is aligned
            WALBERLA_CHECK_EQUAL( (size_t)p % alignment, 0u);
         }
}

void ghostLayerFieldAlignmentTest()
{
   const uint_t xs = 3;
   const uint_t ys = 4;
   const uint_t zs = 2;
   const uint_t fs = 2;
   const uint_t gl = 1;

   const uint_t alignment = 4 * sizeof( double );
   typedef field::AllocateAligned<double,alignment> Allocator;
   shared_ptr<Allocator> alloc = make_shared<Allocator>();
   GhostLayerField<double,fs > field (xs,ys,zs,gl,field::fzyx,alloc);
   for(cell_idx_t f=0; f < cell_idx_c( field.fSize() ); ++f )
      for(cell_idx_t z=0; z < cell_idx_t( field.zSize() ); ++z )
         for(cell_idx_t y=0; y < cell_idx_t( field.ySize() ); ++y )
         {
            void * p = reinterpret_cast<void*>( &(field(0,y,z,f)) );
            // test that each line is aligned
            WALBERLA_CHECK_EQUAL( (size_t)p % alignment, 0u);
         }

}


void blockedIterTest(field::Layout layout)
{
   const uint_t xs = 8;
   const uint_t ys = 9;
   const uint_t zs = 10;
   const uint_t fs = 5;

   Field<uint_t,fs> field(xs,ys,zs,layout);
   for(Field<uint_t,fs>::iterator i = field.begin(); i!= field.end(); ++i)
      field(i.x(),i.y(),i.z(),i.f() ) = uint_c( i.x()+i.y()+i.z()+i.f() );


   int counter = 0;
   for(Field<uint_t,fs>::iterator i   = field.beginSlice(5,5,5,1,7,7,7,3);
                                  i  != field.end(); ++i)
   {
      //cout << counter << "\t(" <<i.f() << "," << i.z() << "," << i.y() << "," << i.x() <<") Value ="<< *i << endl;
      //cout << i << endl;
      counter++;
   }
   WALBERLA_CHECK_EQUAL(counter,2*2*2*2);
}


void ghostLayerFieldCreateAndIterate(field::Layout layout)
{
   const uint_t xs = 3;
   const uint_t ys = 4;
   const uint_t zs = 2;
   const uint_t fs = 2;
   const uint_t gl = 2;

   GhostLayerField<cell_idx_t,fs> field(xs,ys,zs,gl,layout);

   uint_t counter=0;
   for(cell_idx_t f=0; f < cell_idx_c( field.fSize() ); ++f )
      for(cell_idx_t z=-cell_idx_c(gl); z < cell_idx_t( field.zSize() + gl ); ++z )
         for(cell_idx_t y=-cell_idx_c(gl); y < cell_idx_t( field.ySize() + gl ); ++y )
            for(cell_idx_t x=-cell_idx_c(gl); x < cell_idx_t( field.xSize() + gl); ++x )
            {
               field(x,y,z,f) = x*y*z*f + x + y + z + f;
               WALBERLA_CHECK_EQUAL(field(x,y,z,f),x*y*z*f + x + y + z + f);
               counter++;
            }
   WALBERLA_CHECK_EQUAL(counter, (xs+2*gl)* (ys+2*gl) * (zs+2*gl) * fs);



   counter = 0;
   for( auto i = field.begin(); i!= field.end(); ++i )
   {
      //cout << counter << "\t(" <<i.f() << "," << i.z() << "," << i.y() << "," << i.x() <<") Value ="<< *i << endl;
      WALBERLA_CHECK_EQUAL( field(i.x(),i.y(),i.z(),i.f()) , *i );
      counter ++;
   }
   WALBERLA_CHECK_EQUAL(counter, xs*ys*zs*fs);


   counter = 0;
   for( auto i = field.beginWithGhostLayer(); i!= field.end(); ++i )
   {
      WALBERLA_CHECK_EQUAL( field(i.x(),i.y(),i.z(),i.f()) , *i );

      Cell cell ( i.x(), i.y(), i.z() );
      WALBERLA_CHECK_EQUAL( cell, i.cell() );

      counter ++;
   }
   WALBERLA_CHECK_EQUAL(counter, (xs+2*gl) * (ys+2*gl) * (zs+2*gl) *fs);

}


void ghostLayerFieldCreateAndIterate2(field::Layout layout)
{
   const uint_t xs = 6;
   const uint_t ys = 3;
   const uint_t zs = 2;
   const uint_t fs = 19;
   const uint_t gl = 1;

   GhostLayerField<cell_idx_t,fs> field(xs,ys,zs,gl,layout);

   for( auto i = field.beginWithGhostLayer(); i != field.end(); ++i) {
      *i = ( i.x() + i.y() + i.z() + i.f() );
      Cell c = i.cell();

      WALBERLA_CHECK_EQUAL( field( c[0], c[1], c[2], i.f() ), *i );
   }

}

void neighborTest(field::Layout layout)
{
   const uint_t xs = 3;
   const uint_t ys = 4;
   const uint_t zs = 2;
   const uint_t fs = 1;
   const uint_t gl = 2;

   uint_t counter = 0;

   GhostLayerField<uint_t,fs> myField(xs,ys,zs,gl,0,layout);
   for( auto i = myField.begin(); i!= myField.end(); ++i)
   {
      i.neighbor(0,0,1,0) = ++counter;
      WALBERLA_CHECK_EQUAL(i.neighbor(0,0,1,0), myField.get(i.x(),i.y(),i.z()+1));

      i.neighbor(0,1,0,0) = ++counter;
      WALBERLA_CHECK_EQUAL(i.neighbor(0,1,0,0), myField.get(i.x(),i.y()+1,i.z()));

      i.neighbor(1,0,0,0) = ++counter;
      WALBERLA_CHECK_EQUAL(i.neighbor(1,0,0,0), myField.get(i.x()+1,i.y(),i.z()));

      i.neighbor(-2,1,-2,0) = ++counter;
      WALBERLA_CHECK_EQUAL(i.neighbor(2,1,-2,0), myField.get(i.x()+2,i.y()+1,i.z()-2));
   }
}


void ghostlayerIterators(field::Layout layout)
{
   using namespace stencil;

   GhostLayerField<std::string,1> field (1,1,1,1,layout);

   for( auto dir = D3Q27::begin(); dir != D3Q27::end(); ++dir )
      for( auto i = field.beginGhostLayerOnly(*dir); i != field.end(); ++i )
         *i = dirToString[*dir];

   WALBERLA_CHECK_EQUAL( field(-1,-1,-1), "BSW" );
   WALBERLA_CHECK_EQUAL( field(-1,-1, 0),  "SW" );
   WALBERLA_CHECK_EQUAL( field(-1,-1, 1), "TSW" );

   WALBERLA_CHECK_EQUAL( field(-1, 0,-1), "BW"  );
   WALBERLA_CHECK_EQUAL( field(-1, 0, 0),  "W"  );
   WALBERLA_CHECK_EQUAL( field(-1, 0, 1), "TW"  );

   WALBERLA_CHECK_EQUAL( field(-1, 1,-1), "BNW" );
   WALBERLA_CHECK_EQUAL( field(-1, 1, 0),  "NW" );
   WALBERLA_CHECK_EQUAL( field(-1, 1, 1), "TNW" );

   WALBERLA_CHECK_EQUAL( field( 0, 0,-1),   "B" );
   WALBERLA_CHECK_EQUAL( field( 0, 0, 0),   "C" );
   WALBERLA_CHECK_EQUAL( field( 0, 0, 1),   "T" );
   WALBERLA_CHECK_EQUAL( field( 0,-1,-1),  "BS" );
   WALBERLA_CHECK_EQUAL( field( 0,-1, 0),   "S" );
   WALBERLA_CHECK_EQUAL( field( 0,-1, 1),  "TS" );
}




void resizeTest(field::Layout layout)
{
   const uint_t gl = 3;
   GhostLayerField<int,3> field ( 3, 4, 5, gl, layout, make_shared<field::StdFieldAlloc<int> >() );
   WALBERLA_CHECK_EQUAL( field.xSize(), 3);
   WALBERLA_CHECK_EQUAL( field.ySize(), 4);
   WALBERLA_CHECK_EQUAL( field.zSize(), 5);

   WALBERLA_CHECK_EQUAL( field.xAllocSize(), 3+2*gl);
   WALBERLA_CHECK_EQUAL( field.yAllocSize(), 4+2*gl);
   WALBERLA_CHECK_EQUAL( field.zAllocSize(), 5+2*gl);

   field.resize(4,2,1);
   WALBERLA_CHECK_EQUAL( field.xSize(), 4);
   WALBERLA_CHECK_EQUAL( field.ySize(), 2);
   WALBERLA_CHECK_EQUAL( field.zSize(), 1);

   WALBERLA_CHECK_EQUAL( field.xAllocSize(), 4+2*gl);
   WALBERLA_CHECK_EQUAL( field.yAllocSize(), 2+2*gl);
   WALBERLA_CHECK_EQUAL( field.zAllocSize(), 1+2*gl);


   Field<int,3> * p = &field;
   WALBERLA_CHECK_EQUAL( p->xSize(), 4);
   WALBERLA_CHECK_EQUAL( p->ySize(), 2);
   WALBERLA_CHECK_EQUAL( p->zSize(), 1);
   WALBERLA_CHECK_EQUAL( p->xAllocSize(), 4+2*gl);
   WALBERLA_CHECK_EQUAL( p->yAllocSize(), 2+2*gl);
   WALBERLA_CHECK_EQUAL( p->zAllocSize(), 1+2*gl);

   p->resize(3,4,5);
   WALBERLA_CHECK_EQUAL( field.xSize(), 3);
   WALBERLA_CHECK_EQUAL( field.ySize(), 4);
   WALBERLA_CHECK_EQUAL( field.zSize(), 5);

   WALBERLA_CHECK_EQUAL( field.xAllocSize(), 3+2*gl);
   WALBERLA_CHECK_EQUAL( field.yAllocSize(), 4+2*gl);
   WALBERLA_CHECK_EQUAL( field.zAllocSize(), 5+2*gl);
}


void swapTest(field::Layout layout)
{
   GhostLayerField<int, 1> field1(2,2,2,1,1,layout);
   GhostLayerField<int, 1> field2(2,2,2,1,2,layout);
   field1.swapDataPointers(field2);

   WALBERLA_CHECK_EQUAL(field1.get(0,0,0), 2);
   WALBERLA_CHECK_EQUAL(field2.get(0,0,0), 1);
}


void sliceTest(field::Layout layout)
{
   const uint_t fSize = 2;

   const int SLICE_VAL = 42;
   const int ELSE_VAL  = 3;

   Field<int, fSize> field ( 3,3,3, ELSE_VAL, layout );
   CellInterval sliceInterval ( 1,1,1, 2,2,1 );

   // Init the cells in the slice area with SLICE_VAL
   int counter = 0;
   for( auto i = field.beginSliceXYZ( sliceInterval,0 ); i!= field.end(); ++i )
   {
      ++counter;
      *i = SLICE_VAL;
   }
   WALBERLA_CHECK_EQUAL( counter, sliceInterval.xSize()*sliceInterval.ySize()*sliceInterval.zSize() );

   // Iterate over complete field and test for correct values
   counter = 0;
   for( auto i = field.begin(); i != field.end(); ++i)
   {
      ++counter;
      if( i.f() == 0 && sliceInterval.contains( i.cell() ) )
      {
         WALBERLA_CHECK_EQUAL( *i, SLICE_VAL );
      }
      else
      {
         WALBERLA_CHECK_EQUAL( *i, ELSE_VAL);
      }
   }
   WALBERLA_CHECK_EQUAL ( counter, field.xSize() * field.ySize() * field.zSize() * fSize );



   auto sliced = shared_ptr<Field<int,fSize> > ( field.getSlicedField(sliceInterval) );

   WALBERLA_CHECK_EQUAL ( sliced->xSize(), sliceInterval.xSize() );
   WALBERLA_CHECK_EQUAL ( sliced->ySize(), sliceInterval.ySize() );
   WALBERLA_CHECK_EQUAL ( sliced->zSize(), sliceInterval.zSize() );
   WALBERLA_CHECK_EQUAL ( sliced->fSize(), fSize );

   counter = 0;
   for( cell_idx_t z = 0; z < cell_idx_c( sliceInterval.zSize() ); ++z )
      for( cell_idx_t y = 0; y < cell_idx_c( sliceInterval.ySize() ); ++y )
         for( cell_idx_t x = 0; x < cell_idx_c( sliceInterval.xSize() ); ++x )
            for( cell_idx_t f = 0; f < cell_idx_c( fSize ); ++f )
            {
               ++counter;
               if( f == 0 )
               {
                  WALBERLA_CHECK_EQUAL( sliced->get(x,y,z,f), SLICE_VAL );
               }
               else
               {
                  WALBERLA_CHECK_EQUAL( sliced->get(x,y,z,f), ELSE_VAL );
               }
            }
   WALBERLA_CHECK_EQUAL( counter, sliceInterval.xSize()*sliceInterval.ySize()*sliceInterval.zSize() * fSize );


   counter = 0;
   for( auto i = sliced->begin(); i != sliced->end(); ++i )
   {
      ++counter;
      if( i.f() == 0 )
      {
         WALBERLA_CHECK_EQUAL( *i, SLICE_VAL );
      }
      else
      {
         WALBERLA_CHECK_EQUAL( *i, ELSE_VAL );
      }
   }

   WALBERLA_CHECK_EQUAL( counter, sliceInterval.xSize()*sliceInterval.ySize()*sliceInterval.zSize() * fSize );
}



void printerTest()
{
   GhostLayerField<double, 1> field (5,5,5,2, 42.424242);
   field::printSlice( cout, field, 0, 1 ) << endl;
}


void reverseIteratorTest(field::Layout layout)
{
   const uint_t xs = 7;
   const uint_t ys = 6;
   const uint_t zs = 5;
   const uint_t fs = 2;
   const uint_t gl = 3;


   GhostLayerField<cell_idx_t,fs>          field(xs,ys,zs,gl,layout);
   const GhostLayerField<cell_idx_t,fs> &  constField = field;

   // Without ghost layers

   for( auto i = field.begin(); i != field.end(); ++i ) {
      *i = 2*i.x() + 3*i.y() + 7*i.z() + 11*i.f();
   }

   for (auto i = constField.rbegin(); i!= constField.rend(); ++i ) {
      WALBERLA_CHECK_EQUAL ( *i, 2*i.x() + 3*i.y() + 7*i.z() + 11*i.f() );
   }


   // With ghost layers
   for( auto i = field.beginWithGhostLayer(); i != field.end(); ++i ) {
      *i = 2*i.x() + 3*i.y() + 7*i.z() + 11*i.f();
   }

   for (auto i = constField.rbeginWithGhostLayer(); i!= constField.rend(); ++i ) {
      WALBERLA_CHECK_EQUAL ( *i, 2*i.x() + 3*i.y() + 7*i.z() + 11*i.f() );
   }

}

void isIteratorConsecutiveTest ( field::Layout layout )
{
   const uint_t xs = 7;
   const uint_t ys = 6;
   const uint_t zs = 5;
   const uint_t fs = 2;
   const uint_t gl = 3;

   GhostLayerField<int,fs> field(xs,ys,zs,gl,0,layout, make_shared<field::StdFieldAlloc<int> >());
   int * lastPointer = & ( *field.beginWithGhostLayer() );
   for( auto i = ++field.beginWithGhostLayer(); i != field.end(); ++i )
   {
      int * curElement = &(*i);
      WALBERLA_CHECK_EQUAL( curElement - lastPointer , 1 );
      lastPointer = curElement;
   }

}

void swapableCompareTest ( )
{
   typedef Field<unsigned char, 1> MyField;
   std::set<MyField*, field::SwapableCompare<MyField*> >  fieldSet;

   fieldSet.insert( new MyField(2,4,1, field::zyxf) );
   fieldSet.insert( new MyField(1,2,5, field::zyxf) );
   fieldSet.insert( new MyField(1,2,2, field::zyxf) );
   fieldSet.insert( new MyField(1,2,3, field::zyxf) );

   fieldSet.insert( new MyField(2,4,1, field::fzyx) );
   fieldSet.insert( new MyField(1,2,5, field::fzyx) );
   fieldSet.insert( new MyField(1,2,2, field::fzyx) );
   fieldSet.insert( new MyField(1,2,3, field::fzyx) );

   for ( auto i = fieldSet.begin(); i != fieldSet.end(); ++i )
   {
      std::cout << "( " << (*i)->xSize()
                << ","  << (*i)->ySize()
                << ","  << (*i)->zSize()
                << ","  << (*i)->layout() << ")" << std::endl;

   }

   for ( auto i = fieldSet.begin(); i != fieldSet.end(); ++i )
   {
      delete *i;
   }
}


real_t functionTakingAConstIterator( const GhostLayerField<real_t,1>::const_base_iterator & it )
{
   return *it;
}


void iteratorToConstConversionTest()
{
   GhostLayerField<real_t,1> glField ( 1,1,1,1 );
   auto it = glField.begin();
   *it = 1.0; // make sure that iterator is not const
   real_t v = functionTakingAConstIterator( it ); // test conversion to const iterator
   WALBERLA_CHECK_FLOAT_EQUAL( v, 1.0 );
}


void sizeTest()
{
   const uint_t xs = 7;
   const uint_t ys = 6;
   const uint_t zs = 5;
   const uint_t fs = 2;
   const uint_t gl = 1;

   GhostLayerField<cell_idx_t,fs>  field(xs,ys,zs,gl,field::fzyx);
   WALBERLA_CHECK_EQUAL( field.xSize(), xs );
   WALBERLA_CHECK_EQUAL( field.ySize(), ys );
   WALBERLA_CHECK_EQUAL( field.zSize(), zs );

   WALBERLA_CHECK_EQUAL( field.xSizeWithGhostLayer(), xs + 2 * gl );
   WALBERLA_CHECK_EQUAL( field.ySizeWithGhostLayer(), ys + 2 * gl );
   WALBERLA_CHECK_EQUAL( field.zSizeWithGhostLayer(), zs + 2 * gl );
}



void fieldPointerTest()
{
   Field<real_t,1> field ( 2,2,1,1 );

   field( 0,0,0 ) = 0;
   field( 0,1,0 ) = 1;
   field( 1,0,0 ) = 2;
   field( 1,1,0 ) = 3;

   Field<real_t,1>::Ptr fieldPtr ( field, 0,0,0 );

   WALBERLA_CHECK_FLOAT_EQUAL(*fieldPtr, 0.0 );
   WALBERLA_CHECK_FLOAT_EQUAL(fieldPtr.neighbor( stencil::E ), 2.0 );
   WALBERLA_CHECK_FLOAT_EQUAL(fieldPtr.neighbor( stencil::N ), 1.0 );
   WALBERLA_CHECK_FLOAT_EQUAL(fieldPtr.neighbor( stencil::NE ), 3.0 );
}


template<uint_t fSize>
void flattenTest()
{
   Field<Vector3<uint_t>, fSize> field ( 2,2,1 );

   for( cell_idx_t x = 0; x < cell_idx_c(field.xSize()); ++x )
      for( cell_idx_t y = 0; y < cell_idx_c(field.ySize()); ++y )
         for( cell_idx_t z = 0; z < cell_idx_c(field.zSize()); ++z )
            for( cell_idx_t f = 0; f < cell_idx_c(field.fSize()); ++f )
               for( uint_t g = 0; g < 3; ++g )
               {
                  uint_t val = uint_t(&(field( x,y,z,f )[g]));
                  field( x,y,z,f )[g] = val;
               }

   shared_ptr<Field<uint_t, 3*fSize>> flattened(field.flattenedShallowCopy());

   Field<uint_t, 3*fSize> cmp ( 2,2,1 );
   WALBERLA_CHECK_EQUAL(cmp.xSize(), flattened->xSize());
   WALBERLA_CHECK_EQUAL(cmp.ySize(), flattened->ySize());
   WALBERLA_CHECK_EQUAL(cmp.zSize(), flattened->zSize());
   WALBERLA_CHECK_EQUAL(cmp.fSize(), flattened->fSize());
   WALBERLA_CHECK_EQUAL(cmp.xAllocSize(), flattened->xAllocSize());
   WALBERLA_CHECK_EQUAL(cmp.yAllocSize(), flattened->yAllocSize());
   WALBERLA_CHECK_EQUAL(cmp.zAllocSize(), flattened->zAllocSize());
   WALBERLA_CHECK_EQUAL(cmp.fAllocSize(), flattened->fAllocSize());
   WALBERLA_CHECK_EQUAL(cmp.allocSize(), flattened->allocSize());
   WALBERLA_CHECK_EQUAL(cmp.xStride(), flattened->xStride());
   WALBERLA_CHECK_EQUAL(cmp.yStride(), flattened->yStride());
   WALBERLA_CHECK_EQUAL(cmp.zStride(), flattened->zStride());
   WALBERLA_CHECK_EQUAL(cmp.fStride(), flattened->fStride());
   WALBERLA_CHECK_EQUAL(cmp.xOff(), flattened->xOff());
   WALBERLA_CHECK_EQUAL(cmp.yOff(), flattened->yOff());
   WALBERLA_CHECK_EQUAL(cmp.zOff(), flattened->zOff());

   for( cell_idx_t x = 0; x < cell_idx_c(field.xSize()); ++x )
      for( cell_idx_t y = 0; y < cell_idx_c(field.ySize()); ++y )
         for( cell_idx_t z = 0; z < cell_idx_c(field.zSize()); ++z )
            for( cell_idx_t f = 0; f < cell_idx_c(field.fSize()); ++f )
               for( uint_t g = 0; g < 3; ++g )
               {
                  WALBERLA_CHECK_EQUAL(field(x,y,z,f)[g], flattened->get(x,y,z,3*f+cell_idx_c(g)));
               }
}


template<uint_t fSize>
void ghostFlattenTest()
{
   GhostLayerField<Vector3<uint_t>, fSize> field ( 2,2,1, 1 );

   for( cell_idx_t x = -cell_idx_c(field.nrOfGhostLayers()); x < cell_idx_c(field.xSize()+field.nrOfGhostLayers()); ++x )
      for( cell_idx_t y = -cell_idx_c(field.nrOfGhostLayers()); y < cell_idx_c(field.ySize()+field.nrOfGhostLayers()); ++y )
         for( cell_idx_t z = -cell_idx_c(field.nrOfGhostLayers()); z < cell_idx_c(field.zSize()+field.nrOfGhostLayers()); ++z )
            for( cell_idx_t f = 0; f < cell_idx_c(field.fSize()); ++f )
               for( uint_t g = 0; g < 3; ++g )
               {
                  uint_t val = uint_t(&(field( x,y,z,f )[g]));
                  field( x,y,z,f )[g] = val;
               }

   shared_ptr<GhostLayerField<uint_t, 3*fSize>> flattened(field.flattenedShallowCopy());

   GhostLayerField<uint_t, 3*fSize> cmp ( 2,2,1, 1 );
   WALBERLA_CHECK_EQUAL(cmp.xSize(), flattened->xSize());
   WALBERLA_CHECK_EQUAL(cmp.ySize(), flattened->ySize());
   WALBERLA_CHECK_EQUAL(cmp.zSize(), flattened->zSize());
   WALBERLA_CHECK_EQUAL(cmp.fSize(), flattened->fSize());
   WALBERLA_CHECK_EQUAL(cmp.xAllocSize(), flattened->xAllocSize());
   WALBERLA_CHECK_EQUAL(cmp.yAllocSize(), flattened->yAllocSize());
   WALBERLA_CHECK_EQUAL(cmp.zAllocSize(), flattened->zAllocSize());
   WALBERLA_CHECK_EQUAL(cmp.fAllocSize(), flattened->fAllocSize());
   WALBERLA_CHECK_EQUAL(cmp.allocSize(), flattened->allocSize());
   WALBERLA_CHECK_EQUAL(cmp.xStride(), flattened->xStride());
   WALBERLA_CHECK_EQUAL(cmp.yStride(), flattened->yStride());
   WALBERLA_CHECK_EQUAL(cmp.zStride(), flattened->zStride());
   WALBERLA_CHECK_EQUAL(cmp.fStride(), flattened->fStride());
   WALBERLA_CHECK_EQUAL(cmp.xOff(), flattened->xOff());
   WALBERLA_CHECK_EQUAL(cmp.yOff(), flattened->yOff());
   WALBERLA_CHECK_EQUAL(cmp.zOff(), flattened->zOff());
   WALBERLA_CHECK_EQUAL(cmp.nrOfGhostLayers(), flattened->nrOfGhostLayers());

   for( cell_idx_t x = -cell_idx_c(field.nrOfGhostLayers()); x < cell_idx_c(field.xSize()+field.nrOfGhostLayers()); ++x )
      for( cell_idx_t y = -cell_idx_c(field.nrOfGhostLayers()); y < cell_idx_c(field.ySize()+field.nrOfGhostLayers()); ++y )
         for( cell_idx_t z = -cell_idx_c(field.nrOfGhostLayers()); z < cell_idx_c(field.zSize()+field.nrOfGhostLayers()); ++z )
            for( cell_idx_t f = 0; f < cell_idx_c(field.fSize()); ++f )
               for( uint_t g = 0; g < 3; ++g )
               {
                  WALBERLA_CHECK_EQUAL(field(x,y,z,f)[g], flattened->get(x,y,z,3*f+cell_idx_c(g)));
               }
}


int main( int argc, char**argv )
{
   walberla::Environment walberlaEnv( argc, argv );
   using field::fzyx;
   using field::zyxf;

   debug::enterTestMode();
   alignedAllocTest();
   alignedAllocWithOffsetTest();

   alignmentTest();
   ghostLayerFieldAlignmentTest();
   sizeTest();
   iteratorToConstConversionTest();
   fieldPointerTest();


   simpleCreateAndIterate(fzyx);
   simpleCreateAndIterate(zyxf);

   blockedIterTest(fzyx);
   blockedIterTest(zyxf);

   ghostLayerFieldCreateAndIterate(fzyx);
   ghostLayerFieldCreateAndIterate(zyxf);

   ghostLayerFieldCreateAndIterate2(fzyx);
   ghostLayerFieldCreateAndIterate2(zyxf);

   neighborTest(fzyx);
   neighborTest(zyxf);

   ghostlayerIterators(fzyx);
   ghostlayerIterators(zyxf);

   resizeTest(fzyx);
   resizeTest(zyxf);

   swapTest(fzyx);
   swapTest(zyxf);

   sliceTest(fzyx);
   sliceTest(zyxf);

   reverseIteratorTest(fzyx);
   reverseIteratorTest(zyxf);

   isIteratorConsecutiveTest( fzyx );
   isIteratorConsecutiveTest( zyxf );

   flattenTest<1>();
   flattenTest<3>();
   ghostFlattenTest<1>();
   ghostFlattenTest<3>();


   //swapableCompareTest();
   //printerTest();

   return 0;
}


