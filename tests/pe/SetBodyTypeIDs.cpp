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

#include "pe/Materials.h"

#include "pe/rigidbody/SetBodyTypeIDs.h"

#include "pe/rigidbody/Box.h"
#include "pe/rigidbody/Capsule.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/Union.h"

#include "pe/utility/BodyCast.h"

namespace walberla {
namespace pe {

class Base{
public:
    Base(id_t const typeID) : typeID_(typeID) {}
    virtual ~Base() = default;
    id_t getTypeID() const { return typeID_;}
private:
    id_t typeID_;
};

class A : public Base{
    template <class T, int N>
    friend struct SetBodyTypeIDs;
private:
    static id_t staticTypeID_;
    static void setStaticTypeID(id_t typeID) {staticTypeID_ = typeID;}
public:
    A() : Base(getStaticTypeID()) { }

    static id_t getStaticTypeID() {return staticTypeID_;}
};

class B : public Base{
    template <class T, int N>
    friend struct SetBodyTypeIDs;
private:
    static id_t staticTypeID_;
    static void setStaticTypeID(id_t typeID) {staticTypeID_ = typeID;}
public:
    B() : Base(getStaticTypeID()) { }

    static id_t getStaticTypeID() {return staticTypeID_;}
};

class C : public Base{
   template <class T, int N>
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

struct SingleTypeFunctor
{
   std::vector<walberla::id_t> ids_;

   void operator()( const BodyID /*bd*/) { }
   void operator()( const SphereID /*bd*/) { ids_.push_back(Sphere::getStaticTypeID()); }
   void operator()( const BoxID /*bd*/) { ids_.push_back(Box::getStaticTypeID()); }
   void operator()( const CapsuleID /*bd*/) { ids_.push_back(Capsule::getStaticTypeID()); }
};

struct DoubleTypeFunctor
{
   std::vector<walberla::id_t> ids_;

   void operator()( const BodyID /*bd1*/, const BodyID /*bd2*/) { }
   void operator()( const BoxID /*bd1*/, const CapsuleID /*bd2*/) { ids_.push_back(Box::getStaticTypeID()); ids_.push_back(Capsule::getStaticTypeID()); }
   void operator()( const CapsuleID /*bd1*/, const CapsuleID /*bd2*/) { ids_.push_back(5); }
};

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   typedef std::tuple<A, B, C> BodyTuple ;
   SetBodyTypeIDs<BodyTuple>::execute();

   WALBERLA_CHECK_UNEQUAL(A::getStaticTypeID(), 100);
   WALBERLA_CHECK_UNEQUAL(B::getStaticTypeID(), 100);
   WALBERLA_CHECK_UNEQUAL(C::getStaticTypeID(), 100);

   WALBERLA_CHECK_UNEQUAL(A::getStaticTypeID(), B::getStaticTypeID());
   WALBERLA_CHECK_UNEQUAL(A::getStaticTypeID(), C::getStaticTypeID());
   WALBERLA_CHECK_UNEQUAL(B::getStaticTypeID(), C::getStaticTypeID());

   typedef std::tuple<Plane, Sphere, Box, Capsule, Union<Sphere, Box> > BodyTuple2 ;
   SetBodyTypeIDs<BodyTuple2>::execute();
   WALBERLA_CHECK_UNEQUAL(Plane::getStaticTypeID(), Sphere::getStaticTypeID());
   WALBERLA_CHECK_UNEQUAL(Plane::getStaticTypeID(), Box::getStaticTypeID());
   WALBERLA_CHECK_UNEQUAL(Plane::getStaticTypeID(), Capsule::getStaticTypeID());

   SingleTypeFunctor singleFunc;
   Box     bx (0, 0, Vec3(0), Quat(), Vec3(1), Material::find("iron"), false, false, false);
   Capsule cap(0, 0, Vec3(0), Quat(), 1, 1, Material::find("iron"), false, false, false);

   SingleCast<BodyTuple2, SingleTypeFunctor, void>::execute(Box::getStaticTypeID(), singleFunc);
   SingleCast<BodyTuple2, SingleTypeFunctor, void>::execute(Capsule::getStaticTypeID(), singleFunc);
   WALBERLA_CHECK_EQUAL( singleFunc.ids_.size(), 2);
   WALBERLA_CHECK_EQUAL( singleFunc.ids_[0], Box::getStaticTypeID());
   WALBERLA_CHECK_EQUAL( singleFunc.ids_[1], Capsule::getStaticTypeID());

   singleFunc.ids_.clear();
   SingleCast<BodyTuple2, SingleTypeFunctor, void>::execute(&bx, singleFunc);
   SingleCast<BodyTuple2, SingleTypeFunctor, void>::execute(&cap, singleFunc);
   WALBERLA_CHECK_EQUAL( singleFunc.ids_.size(), 2);
   WALBERLA_CHECK_EQUAL( singleFunc.ids_[0], Box::getStaticTypeID());
   WALBERLA_CHECK_EQUAL( singleFunc.ids_[1], Capsule::getStaticTypeID());

   DoubleTypeFunctor doubleFunc;
   DoubleCast<BodyTuple2, BodyTuple2, DoubleTypeFunctor, void>::execute(&bx, &cap, doubleFunc);
   DoubleCast<BodyTuple2, BodyTuple2, DoubleTypeFunctor, void>::execute(&cap, &cap, doubleFunc);
   WALBERLA_CHECK_EQUAL( doubleFunc.ids_.size(), 3);
   WALBERLA_CHECK_EQUAL( doubleFunc.ids_[0], Box::getStaticTypeID());
   WALBERLA_CHECK_EQUAL( doubleFunc.ids_[1], Capsule::getStaticTypeID());;
   WALBERLA_CHECK_EQUAL( doubleFunc.ids_[2], 5);;

   return EXIT_SUCCESS;
}
