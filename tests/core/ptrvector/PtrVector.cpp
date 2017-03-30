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
//! \file PtrVector.cpp
//! \ingroup field
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/ptrvector/PtrVector.h"
#include "core/ptrvector/policies/NoDelete.h"

#include "core/debug/TestSubsystem.h"
   
using namespace walberla;

struct A{
    int c;
    int val;
    static int refCount;

    explicit A(int v = 0) : c(0), val(v) { ++refCount; }
    virtual ~A(){ --refCount; }
};

struct B : public A {
    B() : A(1) {}
};

struct C : public A {
    C() : A(2) {}
};

int A::refCount = 0;

int main( int /*argc*/, char** /*argv*/ )
{
    debug::enterTestMode();
    {
        auto a = new A();
        auto b = new B();
        auto c = new C();

        WALBERLA_CHECK_EQUAL(A::refCount, 3);

        PtrVector<A, PtrDelete> vec;
        vec.pushBack(a);
        vec.pushBack(b);
        vec.pushBack(c);

        WALBERLA_CHECK_EQUAL(A::refCount, 3);

        for (auto it = vec.begin(); it!=vec.end(); ++it){
            ++(it->c);
        }
        for (auto it = vec.begin<B>(); it!=vec.end<B>(); ++it){
            ++(it->c);
        }

        WALBERLA_CHECK_EQUAL(a->c, 1);
        WALBERLA_CHECK_EQUAL(b->c, 2);

        vec.erase(++vec.begin());
        int sum = 0;
        for (auto it = vec.begin(); it!=vec.end(); ++it){
            sum += it->val;
        }

        WALBERLA_CHECK_EQUAL(sum, 2);
    }
    WALBERLA_CHECK_EQUAL(A::refCount, 0);

    auto a = new A();
    auto b = new B();
    WALBERLA_CHECK_EQUAL(A::refCount, 2);
    {
        PtrVector<A, NoDelete> vec;
        vec.pushBack(a);
        vec.pushBack(b);
        WALBERLA_CHECK_EQUAL(A::refCount, 2);
    }
    WALBERLA_CHECK_EQUAL(A::refCount, 2);

    delete a;
    delete b;
}
