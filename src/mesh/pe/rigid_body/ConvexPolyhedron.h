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
//! \file ConvexPolyhedron.h
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <array>
#include "mesh_common/TriangleMeshes.h"
#include "mesh/pe/Types.h"

#include "core/DataTypes.h"
#include "core/math/Constants.h"
#include "core/math/Matrix3.h"
#include "core/math/Vector3.h"
#include "core/DataTypes.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/math/Limits.h"
#include "core/math/Utility.h"

#include "pe/rigidbody/GeomPrimitive.h"
#include "pe/rigidbody/SetBodyTypeIDs.h"
#include "pe/Types.h"

namespace walberla {
namespace mesh {
namespace pe {

using namespace walberla::pe;

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/**
 * \ingroup pe
 * \brief pe body representing a ConvexPolyhedron.
 */
class ConvexPolyhedron : public GeomPrimitive
{
public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   ConvexPolyhedron( id_t sid, id_t uid, const Vec3& gpos, const Quat& q,
                     const TriangleMesh & mesh, MaterialID material,
                     const bool global, const bool communicating, const bool infiniteMass );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ConvexPolyhedron() override;
   //@}
   //**********************************************************************************************
   //**********************************************************************************************

public:
   void init( const Vec3& gpos, const Quat& q,
              const bool global, const bool communicating, const bool infiniteMass );

   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   real_t getVolume() const override;
   real_t getSurfaceArea() const;
   const TriangleMesh & getMesh() const { return mesh_; }
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static inline id_t getStaticTypeID();
   //@}
   //**********************************************************************************************

   //**Output functions****************************************************************************
   /*!\name Output functions */
   //@{
   void print( std::ostream& os, const char* tab ) const override;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline Vec3 support( const Vec3& d ) const override;
   inline Vec3 supportContactThreshold( const Vec3& d ) const override;
   //@}
   //**********************************************************************************************

protected:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   bool containsRelPointImpl ( real_t px, real_t py, real_t pz ) const override;
   bool isSurfaceRelPointImpl( real_t px, real_t py, real_t pz ) const override;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void calcBoundingBox() override;  // Calculation of the axis-aligned bounding box
   virtual TriangleMesh::VertexHandle supportVertex( const TriangleMesh::Normal & d, const TriangleMesh::VertexHandle startVertex ) const;
   //@}
   //**********************************************************************************************

   TriangleMesh mesh_;
   real_t boundingSphereRadius_;
   std::array<TriangleMesh::VertexHandle, 8> octandVertices_;
private:
   static id_t staticTypeID_;  //< type id of ConvexPolyhedron, will be set by SetBodyTypeIDs
   static void setStaticTypeID(id_t typeID) {staticTypeID_ = typeID;}

   //** friend declaration
   /// needed to be able to set static type ids with setStaticTypeID
   template <class T, int N>
   friend struct walberla::pe::SetBodyTypeIDs;
};
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Returns unique type id of this type
 *
 * \return geometry specific type id
 */
inline id_t ConvexPolyhedron::getStaticTypeID()
{
   return staticTypeID_;
}
//*************************************************************************************************



//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name ConvexPolyhedron operators */
//@{
std::ostream& operator<<( std::ostream& os, const ConvexPolyhedron& s );
std::ostream& operator<<( std::ostream& os, ConstConvexPolyhedronID s );
//@}
//*************************************************************************************************


} // namespace pe
} // namespace mesh
} // namespace walberla