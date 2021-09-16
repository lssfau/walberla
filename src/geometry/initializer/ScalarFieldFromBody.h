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
//! \file ScalarFieldFromBody.h
//! \ingroup geometry
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#pragma once

#include "geometry/initializer/Initializer.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "core/math/Parser.h"


namespace walberla {
namespace geometry {
namespace initializer {



   //*******************************************************************************************************************
   /*! Initializes a scalar field from a geometric body
   *
   * Currently supported are Sphere, Ellipsoid and Box (= AABB)
   *
   * Examples:
   * \verbatim
       <InitializerUID> {
            initialFill : drop;
            someArbitraryId {
               shape sphere;
               add;                  // add to the current value in the field
               value 1.0;
               midpoint < 3,4,5>;
               radius 4;
               id 0;                 // If given a vector of scalar fields, which one of them to operate on.
                                     // Should be zero (the default value) if given a scalar field directly.
            }
            object2 {
               shape box;
               add;
               value -3.14;
               min <1,2,3>;
               max <3,4,5>;
            }
            object3_ellipse {
               set;                  // overwrite the current value in the field
               value 0;
               shape ellipsoid;
               midpoint < 3,4,2>;
               axis1    <1,0,0>;
               axis2    <0,1,0>;
               radii    <1,1,4>;
            }
       }
     \endverbatim
   *
   * \ingroup geometry
   *
   */
   //*******************************************************************************************************************
   template <typename Field_T>
   class ScalarFieldFromBody : public Initializer
   {
   public:
      using Value_T = typename Field_T::value_type;
      
      /*************************************************************************************************************//**
      * Constructor
      *
      * \param scalarFieldID    the scalar field to initialize,
                                or a vector of scalar fields to initialize based on the id specified with the body
      *
      *****************************************************************************************************************/
      ScalarFieldFromBody( StructuredBlockStorage & structuredBlockStorage, BlockDataID scalarFieldID )
         : ScalarFieldFromBody(structuredBlockStorage, std::vector<BlockDataID>(1, scalarFieldID))
      {}
      ScalarFieldFromBody( StructuredBlockStorage & structuredBlockStorage, std::vector<BlockDataID> scalarFieldID );



      /*************************************************************************************************************//**
      * Interface implementation for Initializer - sets a body on a scalar field with options from configuration file
      *
      *****************************************************************************************************************/
      void init( BlockStorage & blockStorage, const Config::BlockHandle & blockHandle ) override;



      /*************************************************************************************************************//**
      * Sets a body on the scalar field
      *
      * \param body       The body object - has to implement either overlapFraction(...), or contains(...)
      *                   see BodyOverlapFunctions for detailed body concept
      * \param value      The value to set on the matched cells in the field.
      * \param addOrSet   If true, the value is added to scalar field
      *                   If false, the value is set on the scalar field.
      * \param id         If operating on a vector of fields, which field to treat. Zero otherwise.
      *
      *  Supported bodies are Sphere, Ellipsoid, AABB.
      *  To add a new supported body implement concept defined in BodyOverlapFunctions.h, and
      *  add an explicit template instantiation in ScalarFieldFromBody.cpp for the new body.
      *
      *****************************************************************************************************************/
      template<typename Body>
      void init( const Body & body, Value_T value, bool addOrSet, std::vector<BlockDataID>::size_type id = 0 );
      /*************************************************************************************************************//**
      * Sets a body on the scalar field
      *
      * \param body       The body object - has to implement either overlapFraction(...), or contains(...)
      *                   see BodyOverlapFunctions for detailed body concept
      * \param parser     A function parser which will have the variables x,y,z bound before it is evaluated
      * \param addOrSet   If true, the value is added to scalar field
      *                   If false, the value is set on the scalar field.
      * \param id         If operating on a vector of fields, which field to treat. Zero otherwise.
      *
      *  Supported bodies are Sphere, Ellipsoid, AABB.
      *  To add a new supported body implement concept defined in BodyOverlapFunctions.h, and
      *  add an explicit template instantiation in ScalarFieldFromBody.cpp for the new body.
      *
      *****************************************************************************************************************/
      template<typename Body>
      void init( const Body & body, math::FunctionParser & parser, bool addOrSet, std::vector<BlockDataID>::size_type id = 0 );


   protected:

      StructuredBlockStorage & structuredBlockStorage_;
      std::vector<BlockDataID> scalarFieldID_;

      std::string              addKeyword_;
      std::string              setKeyword_;
   };




} // namespace initializer
} // namespace geometry
} // namespace walberla

#include "ScalarFieldFromBody.impl.h"

