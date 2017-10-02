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
//! \file Types.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/AABBFwd.h"
#include "stencil/D3Q27.h"

#include <boost/array.hpp>

#include <vector>

namespace walberla{
namespace math {

template<typename T> class Vector2;
template<typename T> class Vector3;
template<typename T> class Matrix2;
template<typename T> class Matrix3;
template<typename T> class Quaternion;
template<typename T> class MatrixMxN;

}
namespace pe{

typedef math::Vector2<real_t> Vec2;
typedef math::Vector3<real_t> Vec3;
typedef math::Matrix2<real_t> Mat2;
typedef math::Matrix3<real_t> Mat3;
typedef math::Quaternion<real_t>  Quat;
typedef math::MatrixMxN<real_t> MatN;

using math::AABB;

//=================================================================================================
//
//  walberla::pe NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

class Attachable;
class BallJoint;
class BodyManager;
class BodyStorage;
class Box;
class Capsule;
class Contact;
class Cylinder;
class CylindricalBoundary;
class FixedJoint;
class ForceGenerator;
class GeomPrimitive;
class Gravity;
class HingeJoint;
class Joint;
class Link;
class Material;
class MPISystem;
class Node;
class Plane;
class Process;
class RigidBody;
class Section;
class SliderJoint;
class Sphere;
class Spring;
class Squirmer;
class TriangleMesh;
template <typename BodyTypeTuple>
class Union;

typedef RigidBody             BodyType;            //!< Type of the rigid bodies.
typedef RigidBody*            BodyID;              //!< Handle for a rigid body.
typedef const RigidBody*      ConstBodyID;         //!< Handle for a constant rigid body.

typedef GeomPrimitive*  GeomID;
typedef const GeomPrimitive*  ConstGeomID;

typedef Box                   BoxType;             //!< Type of the box geometric primitive.
typedef Box*                  BoxID;               //!< Handle for a box primitive.
typedef const Box*            ConstBoxID;          //!< Handle for a constant box primitive.

typedef Capsule               CapsuleType;         //!< Type of the capsule geometric primitive.
typedef Capsule*              CapsuleID;           //!< Handle for a capsule primitive.
typedef const Capsule*        ConstCapsuleID;      //!< Handle for a constant capsule primitive.

typedef Cylinder              CylinderType;        //!< Type of the cylinder geometric primitive.
typedef Cylinder*             CylinderID;          //!< Handle for a cylinder primitive.
typedef const Cylinder*       ConstCylinderID;     //!< Handle for a constant cylinder primitive.

typedef CylindricalBoundary        CylindricalBoundaryType;        //!< Type of the cylindrical boundary geometric primitive.
typedef CylindricalBoundary*       CylindricalBoundaryID;          //!< Handle for a cylindrical boundary primitive.
typedef const CylindricalBoundary* ConstCylindricalBoundaryID;     //!< Handle for a constant cylindrical boundary primitive.

typedef Plane                 PlaneType;           //!< Type of the plane geometric primitive.
typedef Plane*                PlaneID;             //!< Handle for a plane primitive.
typedef const Plane*          ConstPlaneID;        //!< Handle for a constant plane primitive.

typedef Sphere                SphereType;          //!< Type of the sphere geometric primitive.
typedef Sphere*               SphereID;            //!< Handle for a sphere primitive.
typedef const Sphere*         ConstSphereID;       //!< Handle for a constant sphere primitive.

typedef Squirmer              SquirmerType;        //!< Type of the squirmer geometric primitive.
typedef Squirmer*             SquirmerID;          //!< Handle for a squirmer primitive.
typedef const Squirmer*       ConstSquirmerID;     //!< Handle for a constant squirmer primitive.

typedef TriangleMesh          MeshType;             //!< Type of the triangle mesh geometric primitive.
typedef TriangleMesh*         MeshID;               //!< Handle for a triangle mesh primitive.
typedef const TriangleMesh*   ConstMeshID;          //!< Handle for a constant triangle mesh primitive.

typedef Attachable            AttachableType;      //!< Type of the attachables.
typedef Attachable*           AttachableID;        //!< Handle for an attachable.
typedef const Attachable*     ConstAttachableID;   //!< Handle for a constant attachable.

typedef Gravity               GravityType;         //!< Type of the gravity force generators.
typedef Gravity*              GravityID;           //!< Handle for a gravity force generator.
typedef const Gravity*        ConstGravityID;      //!< Handle for a constant gravity force generator.

typedef Spring                SpringType;          //!< Type of the spring force generators.
typedef Spring*               SpringID;            //!< Handle for a spring force generator.
typedef const Spring*         ConstSpringID;       //!< Handle for a constant spring force generator.

typedef Contact               ContactType;         //!< Type of the contacts.
typedef Contact*              ContactID;           //!< Handle for a contact.
typedef const Contact*        ConstContactID;      //!< Handle for a constant contact.

typedef Joint                 JointType;           //!< Type of the joints.
typedef Joint*                JointID;             //!< Handle for a joint.
typedef const Joint*          ConstJointID;        //!< Handle for a constant joint.

typedef SliderJoint           SliderJointType;     //!< Type of the slider joint.
typedef SliderJoint*          SldierJointID;       //!< Handle for a slider joint.
typedef const SliderJoint*    ConstSliderJointID;  //!< Handle for a constant slider joint.

typedef HingeJoint            HingeJointType;      //!< Type of the hinge joint.
typedef HingeJoint*           HingeJointID;        //!< Handle for a hinge joint.
typedef const HingeJoint*     ConstHingeJointID;   //!< Handle for a constant hinge joint.

typedef BallJoint             BallJointType;       //!< Type of the ball joint.
typedef BallJoint*            BallJointID;         //!< Handle for a ball joint.
typedef const BallJoint*      ConstBallJointID;    //!< Handle for a constant ball joint.

typedef Node                  NodeType;
typedef Node*                 NodeID;       //!< Handle to a BodyTrait instance.
typedef const Node*           ConstNodeID;  //!< Handle to a constant BodyTrait instance.

typedef Process               ProcessType;         //!< Type of the remote processes.
typedef Process*              ProcessID;           //!< Handle for a remote process.
typedef const Process*        ConstProcessID;      //!< Handle for a constant remote process.


typedef BodyManager*          ManagerID;           //!< Handle for a BodyManager.
typedef const BodyManager*    ConstManagerID;      //!< Handle for a constant BodyManager.

typedef std::vector<Material>  Materials;          //!< Vector for materials.


//*************************************************************************************************
/*!\brief Unique material ID.
 *
 * Every registered material has a unique MaterialID that can be used wherever the material is
 * required. The \b pe engine provides a couple of predefined materials (see the \ref materials
 * module). However, it is possible to define a custom material via the createMaterial() function:

   \code
   // Creates the material "myMaterial" with the following material properties:
   //  - material density               : 2.54
   //  - coefficient of restitution     : 0.8
   //  - coefficient of static friction : 0.1
   //  - coefficient of dynamic friction: 0.05
   MaterialID myMaterial = createMaterial( "myMaterial", real_c(2.54), real_c(0.8), real_c(0.1), real_c(0.05) );
   \endcode
 */
typedef Materials::size_type  MaterialID;
//*************************************************************************************************

///Output data type of coarse collision detection module
typedef std::vector<std::pair<BodyID, BodyID> > PossibleContacts;

///Output data type of fine collision detection module
typedef std::vector< Contact > Contacts;

///Class enum for BodyStorages
/// \see Storage
struct StorageType
{
   enum StorageTypeEnum
   {
      LOCAL    = 0,
      SHADOW   = 1
   };
};

struct StorageSelect
{
   enum StorageSelectEnum
   {
      LOCAL    = 1 << 0,
      SHADOW   = 1 << 1,
      GLOBAL   = 1 << 2
   };
};

///Container for local and shadow body storage
typedef boost::array<BodyStorage, 2> Storage;

}  // namespace pe
}  // namespace walberla
