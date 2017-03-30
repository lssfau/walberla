//=================================================================================================
/*!
 *  \file   StencilTest.cpp
 *  \brief  
 *  \author Martin Bauer
 */
//=================================================================================================


#include "stencil/D3Q27.h"
#include "stencil/D3Q19.h"
#include "stencil/D2Q4.h"
#include "stencil/D2Q5.h"

#include "core/debug/TestSubsystem.h"

#include <iostream>

using namespace walberla;
using namespace stencil;

using std::cout;
using std::endl;


template <typename S>
void stencilIteration()
{
   uint_t counter=0;
   for( typename S::iterator i = S::begin(); i !=  S::end(); ++i) {
      counter++;
   }

   // the "+" is required, because WALBERLA_CHECK__EQUAL takes a reference, and Size is a constant
   WALBERLA_CHECK_EQUAL(counter,+S::Size);
}


void directionTest()
{
   for( D3Q27::iterator i = D3Q27::begin(); i != D3Q27::end(); ++i) {
      WALBERLA_CHECK_EQUAL(i.cx(), - cx[i.mirrorX()]);
      WALBERLA_CHECK_EQUAL(i.cy(), - cy[i.mirrorY()]);
      WALBERLA_CHECK_EQUAL(i.cz(), - cz[i.mirrorZ()]);

      WALBERLA_CHECK_EQUAL(-cx[i.inverseDir()], i.cx());
      WALBERLA_CHECK_EQUAL(-cy[i.inverseDir()], i.cy());
      WALBERLA_CHECK_EQUAL(-cz[i.inverseDir()], i.cz());
   }
}


int main()
{
   debug::enterTestMode();

   stencilIteration<D3Q27>();
   stencilIteration<D3Q19>();
   stencilIteration<D2Q4>();
   stencilIteration<D2Q5>();


   directionTest();
}


