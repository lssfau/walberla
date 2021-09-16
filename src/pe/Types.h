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

#include <array>
#include <memory>
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

using Vec2 = math::Vector2<real_t>;
using Vec3 = math::Vector3<real_t>;
using Mat2 = math::Matrix2<real_t>;
using Mat3 = math::Matrix3<real_t>;
using Quat = math::Quaternion<real_t>;
using MatN = math::MatrixMxN<real_t>;

using math::AABB;

//=================================================================================================
//
//  walberla::pe NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

class BodyManager;
class BodyStorage;
class Box;
class Capsule;
class Contact;
class Cylinder;
class Ellipsoid;
class CylindricalBoundary;
class GeomPrimitive;
class Material;
class MPISystem;
class Plane;
class RigidBody;
class Sphere;
class Spring;
class Squirmer;
class TriangleMesh;
template <typename... BodyTypes>
class Union;

using BodyType = RigidBody;            //!< Type of the rigid bodies.
using BodyID = RigidBody *;              //!< Handle for a rigid body.
using ConstBodyID = const RigidBody *;         //!< Handle for a constant rigid body.
using   BodyPtr             = std::unique_ptr<RigidBody>;

using GeomID = GeomPrimitive *;
using ConstGeomID = const GeomPrimitive *;

using BoxType = Box;             //!< Type of the box geometric primitive.
using BoxID = Box *;               //!< Handle for a box primitive.
using ConstBoxID = const Box *;          //!< Handle for a constant box primitive.
using   BoxPtr              = std::unique_ptr<Box>;

using CapsuleType = Capsule;         //!< Type of the capsule geometric primitive.
using CapsuleID = Capsule *;           //!< Handle for a capsule primitive.
using ConstCapsuleID = const Capsule *;      //!< Handle for a constant capsule primitive.
using   CapsulePtr          = std::unique_ptr<Capsule>;

using CylinderType = Cylinder;        //!< Type of the cylinder geometric primitive.
using CylinderID = Cylinder *;          //!< Handle for a cylinder primitive.
using ConstCylinderID = const Cylinder *;     //!< Handle for a constant cylinder primitive.
using   CylinderPtr         = std::unique_ptr<Cylinder>;

using CylindricalBoundaryType = CylindricalBoundary;        //!< Type of the cylindrical boundary geometric primitive.
using CylindricalBoundaryID = CylindricalBoundary *;          //!< Handle for a cylindrical boundary primitive.
using ConstCylindricalBoundaryID = const CylindricalBoundary *;     //!< Handle for a constant cylindrical boundary primitive.
using   CylindricalBoundaryPtr   = std::unique_ptr<CylindricalBoundary>;

using EllipsoidType = Ellipsoid;       //!< Type of the ellipsoid geometric primitive.
using EllipsoidID = Ellipsoid *;         //!< Handle for a ellipsoid primitive.
using ConstEllipsoidID = const Ellipsoid *;    //!< Handle for a constant ellipsoid primitive.
using   EllipsoidPtr        = std::unique_ptr<Ellipsoid>;

using PlaneType = Plane;           //!< Type of the plane geometric primitive.
using PlaneID = Plane *;             //!< Handle for a plane primitive.
using ConstPlaneID = const Plane *;        //!< Handle for a constant plane primitive.
using   PlanePtr            = std::unique_ptr<Plane>;

using SphereType = Sphere;          //!< Type of the sphere geometric primitive.
using SphereID = Sphere *;            //!< Handle for a sphere primitive.
using ConstSphereID = const Sphere *;       //!< Handle for a constant sphere primitive.
using   SpherePtr           = std::unique_ptr<Sphere>;

using SquirmerType = Squirmer;        //!< Type of the squirmer geometric primitive.
using SquirmerID = Squirmer *;          //!< Handle for a squirmer primitive.
using ConstSquirmerID = const Squirmer *;     //!< Handle for a constant squirmer primitive.
using   SquirmerPtr         = std::unique_ptr<Squirmer>;

using MeshType = TriangleMesh;             //!< Type of the triangle mesh geometric primitive.
using MeshID = TriangleMesh *;               //!< Handle for a triangle mesh primitive.
using ConstMeshID = const TriangleMesh *;          //!< Handle for a constant triangle mesh primitive.
using   TriangleMeshPtr     = std::unique_ptr<TriangleMesh>;

using ContactType = Contact;         //!< Type of the contacts.
using ContactID = Contact *;           //!< Handle for a contact.
using ConstContactID = const Contact *;      //!< Handle for a constant contact.

using ManagerID = BodyManager *;           //!< Handle for a BodyManager.
using ConstManagerID = const BodyManager *;      //!< Handle for a constant BodyManager.

using Materials = std::vector<Material>;          //!< Vector for materials.


//*************************************************************************************************
/*!\brief Unique material ID.
 *
 * Every registered material has a unique MaterialID that can be used wherever the material is
 * required. The \b pe engine provides a couple of predefined materials (see the \ref Materials
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
using MaterialID = Materials::size_type;
//*************************************************************************************************

///Output data type of coarse collision detection module
using PossibleContacts = std::vector<std::pair<BodyID, BodyID>>;

///Output data type of fine collision detection module
using Contacts = std::vector<Contact>;

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
using Storage = std::array<BodyStorage, 2>;

}  // namespace pe
}  // namespace walberla
