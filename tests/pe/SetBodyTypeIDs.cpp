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
//! \file SetBodyTypeIDs.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"

#include "pe/rigidbody/SetBodyTypeIDs.h"

#include "pe/rigidbody/Box.h"
#include "pe/rigidbody/Capsule.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Union.h"

namespace walberla {
namespace pe {

class Base{
public:
    Base(id_t const typeID) : typeID_(typeID) {}
    virtual ~Base() {}
    id_t getTypeID() const { return typeID_;}
private:
    id_t typeID_;
};

class A : public Base{
    template <class T>
    friend struct SetBodyTypeIDs;
private:
    static id_t staticTypeID_;
    static void setStaticTypeID(id_t typeID) {staticTypeID_ = typeID;}
public:
    A() : Base(getStaticTypeID()) { }

    static id_t getStaticTypeID() {return staticTypeID_;}
};

class B : public Base{
    template <class T>
    friend struct SetBodyTypeIDs;
private:
    static id_t staticTypeID_;
    static void setStaticTypeID(id_t typeID) {staticTypeID_ = typeID;}
public:
    B() : Base(getStaticTypeID()) { }

    static id_t getStaticTypeID() {return staticTypeID_;}
};

class C : public Base{
   template <class T>
   friend struct SetBodyTypeIDs;
private:
    static id_t staticTypeID_;
    static void setStaticTypeID(id_t typeID) {staticTypeID_ = typeID;}
public:
    C() : Base(getStaticTypeID()) { }

    static id_t getStaticTypeID() {return staticTypeID_;}
};

id_t A::staticTypeID_ = 100;
id_t B::staticTypeID_ = 100;
id_t C::staticTypeID_ = 100;

}
}

using namespace walberla::pe;

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   typedef boost::tuple<A, B, C> BodyTuple ;
   SetBodyTypeIDs<BodyTuple>::execute();

   WALBERLA_CHECK_UNEQUAL(A::getStaticTypeID(), 100);
   WALBERLA_CHECK_UNEQUAL(B::getStaticTypeID(), 100);
   WALBERLA_CHECK_UNEQUAL(C::getStaticTypeID(), 100);

   WALBERLA_CHECK_UNEQUAL(A::getStaticTypeID(), B::getStaticTypeID());
   WALBERLA_CHECK_UNEQUAL(A::getStaticTypeID(), C::getStaticTypeID());
   WALBERLA_CHECK_UNEQUAL(B::getStaticTypeID(), C::getStaticTypeID());

   typedef boost::tuple<Plane, Sphere, Box, Capsule, Union< boost::tuple<Sphere, Box> > > BodyTuple2 ;
   SetBodyTypeIDs<BodyTuple2>::execute();
   WALBERLA_CHECK_UNEQUAL(Plane::getStaticTypeID(), Sphere::getStaticTypeID());
   WALBERLA_CHECK_UNEQUAL(Plane::getStaticTypeID(), Box::getStaticTypeID());
   WALBERLA_CHECK_UNEQUAL(Plane::getStaticTypeID(), Capsule::getStaticTypeID());

   return EXIT_SUCCESS;
}
