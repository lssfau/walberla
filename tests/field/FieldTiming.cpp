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
//! \file FieldTiming.cpp
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Timing of various iteration variants for fields
//
//======================================================================================================================

#include "field/Field.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Random.h"
#include "core/timing/TimingPool.h"

#include <fstream>
#include <iostream>
#include <sstream>


namespace walberla {

using namespace field;

using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::stringstream;

const uint_t fSize = 19;

typedef Field<double,fSize > DoubleField;


double sumHandWritten(const DoubleField & field)
{
   double sum = 0;

   WALBERLA_CHECK( field.layout()== field::fzyx );

   cell_idx_t sf = cell_idx_c( field.fSize() );
   cell_idx_t sz = cell_idx_c( field.zSize() );
   cell_idx_t sy = cell_idx_c( field.ySize() );

   for(cell_idx_t f=0; f < sf; ++f )
      for(cell_idx_t z=0; z < sz; ++z )
         for(cell_idx_t y=0; y < sy; ++y )
         {
            const double* begin = &field.get(0,y,z,f);
            const double* end   = &field.get(cell_idx_c(field.xSize())-1,y,z,f)+1;
            for(; begin!=end; ++begin) {
               sum += *begin;
            }
         }

   return sum;
}

double sumHandWrittenPlusIf(const DoubleField & field)
{
   double sum = 0;

   cell_idx_t sf = cell_idx_c( field.fSize() );
   cell_idx_t sz = cell_idx_c( field.zSize() );
   cell_idx_t sy = cell_idx_c( field.ySize() );

   for(cell_idx_t f=0; f < sf; ++f )
      for(cell_idx_t z=0; z < sz; ++z )
         for(cell_idx_t y=0; y < sy; ++y )
         {
            const double* begin = &field.get(0,y,z,f);
            const double* end   = &field.get(cell_idx_c(field.xSize())-1,y,z,f)+1;
            for(; begin!=end; ++begin) {
               sum += *begin;
               if(begin == nullptr) // Artificial if condition to simulate iterators
                  cout << "Should not happen" << endl;
            }
         }

   return sum;
}



double sumIterator(const DoubleField & field)
{
   //WALBERLA_CHECK_EQUAL(field.layout() == fzyx );
   double sum = 0;

   for(DoubleField::const_iterator i = field.begin(); i != field.end(); ++i)
      sum += *i;

   return sum;
}

double sumIteratorCachedEnd(const DoubleField & field)
{
   //WALBERLA_CHECK_EQUAL(field.layout() == fzyx );
   double sum = 0;

   const DoubleField::const_iterator& myEnd = field.end();
   for(DoubleField::const_iterator i = field.begin(); i != myEnd; ++i)
      sum += *i;

   return sum;
}

double sumIteratorFast(const DoubleField & field)
{
   double sum = 0;
   for( DoubleField::const_iterator i = field.begin(); i != field.end(); i.incrOuter() )
      for (; i.testInner(); i.incrInner() )
         sum += *i;


   return sum;
}

double sumGetFzyx(const DoubleField & field)
{
   //WALBERLA_CHECK_EQUAL(field.layout() == fzyx );
   double sum = 0;

   cell_idx_t sf = cell_idx_c( field.fSize() );
   cell_idx_t sx = cell_idx_c( field.xSize() );
   cell_idx_t sy = cell_idx_c( field.ySize() );
   cell_idx_t sz = cell_idx_c( field.zSize() );

   for(cell_idx_t f=0; f < sf; ++f )
      for(cell_idx_t z=0; z < sz; ++z )
         for(cell_idx_t y=0; y < sy; ++y )
            for(cell_idx_t x=0; x < sx; ++x )
               sum += field.get(x,y,z,f);

   return sum;
}

double sumGetZyxf(const DoubleField & field)
{
   //WALBERLA_CHECK_EQUAL(field.layout() == fzyx );
   double sum = 0;

   cell_idx_t sf = cell_idx_c( field.fSize() );
   cell_idx_t sx = cell_idx_c( field.xSize() );
   cell_idx_t sy = cell_idx_c( field.ySize() );
   cell_idx_t sz = cell_idx_c( field.zSize() );

   for(cell_idx_t z=0; z < sz; ++z )
      for(cell_idx_t y=0; y < sy; ++y )
         for(cell_idx_t x=0; x < sx; ++x )
            for(cell_idx_t f=0; f < sf; ++f )
               sum += field.get(x,y,z,f);

   return sum;
}



double xyzfTest ( const DoubleField & field )
{
   double sum = 0;
   double sum2 = 0;
   const DoubleField::const_iterator& myEnd = field.end();
   for(DoubleField::const_iterator i = field.begin(); i != myEnd; ++i)
   {
      sum += *i;
      sum2 += i.x() + i.y() + i.z() + i.f();
   }

   if ( isIdentical(sum2,42.424242) )
      cout << "Test only for making sure that sum2 is not optimized away";

   return sum;
}


double initFieldRandom(DoubleField & field)
{
   double sum =0;
   for(Field<double,fSize>::iterator i = field.begin(); i != field.end(); ++i) {
      *i = math::realRandom( 1.0, 2.0 );
      sum += *i;
   }
   return sum;
}

using Func = std::function<double (const DoubleField &)>;
void timeFunction( WcTimer & timer,  Func f, const DoubleField & field, double sum, double epsilon, int nrExecutions=30 )
{
   for(int i=0 ; i < nrExecutions; ++i)
   {
      timer.start();
      double val = f(field);
      WALBERLA_CHECK( floatIsEqual(sum, val, epsilon ) );
      if(val < -10000)
         cout << "Shouldn't happen - only there that function call cannot be optimized away" << endl;

      timer.end();
   }
}


int main(int argc, char ** argv)
{
   debug::enterTestMode();

#ifdef NDEBUG
   uint_t xs = 90;
   uint_t ys = 70;
   uint_t zs = 50;
#else
   uint_t xs = 20;
   uint_t ys = 10;
   uint_t zs = 10;
#endif

   if(argc > 1)
   {
      if( argc != 4) {
         cerr << "Usage: " << argv[0] << " xSize ySize zSize" << endl;
         return 1;
      }
      xs = uint_c( atoi(argv[1]) );
      ys = uint_c( atoi(argv[2]) );
      zs = uint_c( atoi(argv[3]) );
   }

   Field<double,fSize > field (xs,ys,zs,zyxf);
   if ( field.layout() == fzyx )
      cout << "Test for layout fzyx"<< endl;
   else
      cout << "Test for layout zyxf"<< endl;

   const double expectedValue = double_c( xs ) * double_c( ys ) * double_c( zs ) * double_c( fSize ) * 1.5;
   const double epsilon       = expectedValue * 10e-12;

   double sum = initFieldRandom(field);

   WcTimingPool tp;
   timeFunction ( tp["Get Zyxf"],  sumGetZyxf,    field, sum, epsilon );
   timeFunction ( tp["Get Fzyx"],  sumGetFzyx,    field, sum, epsilon );
   //timeFunction ( tp["Handwritten"],       sumHandWritten,       field, sum, epsilon );
   timeFunction ( tp["Iterator"],          sumIterator,          field, sum, epsilon );
   timeFunction ( tp["Iterator Fast"],     sumIteratorFast,      field, sum, epsilon );

   // Print timing results
   cout << tp;

   stringstream fileName;
   fileName << "timingResults_" << xs << "_" << ys << "_" << zs << ".m";
   ofstream of ( fileName.str().c_str() );
   of << "size = [ " << xs << "," << ys << "," << zs <<"];" << endl;
   tp.printMatlab( of );

   return EXIT_SUCCESS;
}
}

int main(int argc, char** argv){
   return walberla::main(argc, argv);
}
