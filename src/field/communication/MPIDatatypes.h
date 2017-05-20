
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
//! \file MPIDatatypes.h
//! \ingroup field
//! \author Christian Godenschwager
//! \brief Functions to generate MPI data types for fields
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/mpi/MPIWrapper.h"

#include "field/Field.h"
#include "field/GhostLayerField.h"

#include "stencil/Directions.h"

namespace walberla {
namespace field {
namespace communication {

/*
 * All function in this file work for types "Field_T" and "GhostLayerField_T" which implement
 * the following functions:
 *
 * Trait Field_T:
 *     xOff(), yOff(), zOff()
 *     xAllocSize(), yAllocSize(), zAllocSize(), fAllocSize()
 *     xSize(), ySize(), zSize(), fSize()
 * Trait GhostLayerField_T,  same as Field_T plus additional:
 *     getGhostRegion()
 *     getSliceBeforeGhostLayer()
 *
 * For the semantics of these functions see documentation of Field / GhostLayerField
 */



//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate the contents of a whole Field
*
*  Only the inner cells of a Field are communicated, when the returned datatype is used in MPI communication.
*  
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field        The Field to be communicated
*  \returns            The MPI datatype
*/
//======================================================================================================================
template<typename Field_T>
MPI_Datatype mpiDatatype( const Field_T & field );



//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate an interval of cells of a Field
*
*  Only the cells in the described interval of a Field are communicated, when the returned datatype is used in MPI
*  communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field        The Field to be communicated
*  \param xBeg         Beginning of the interval in x direction
*  \param yBeg         Beginning of the interval in y direction
*  \param zBeg         Beginning of the interval in z direction
*  \param fBeg         Beginning of the interval in f direction
*  \param xEnd         End of the interval in x direction (included in the interval)
*  \param yEnd         End of the interval in y direction (included in the interval)
*  \param zEnd         End of the interval in z direction (included in the interval)
*  \param fEnd         End of the interval in f direction (included in the interval)
*  \returns            The MPI datatype
*/
//======================================================================================================================
template<typename Field_T>
MPI_Datatype mpiDatatypeSlice( const Field_T & field,
                               const cell_idx_t xBeg, const cell_idx_t yBeg, const cell_idx_t zBeg, const cell_idx_t fBeg,
                               const cell_idx_t xEnd, const cell_idx_t yEnd, const cell_idx_t zEnd, const cell_idx_t fEnd );



//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate an interval of cells of a Field
*
*  Only the cells in the described interval of a Field are communicated, when the returned datatype is used in MPI
*  communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field        The Field to be communicated
*  \param interval     Interval of communicated cells
*  \param f            f component that is communicated
*  \returns            The MPI datatype
*/
//======================================================================================================================
template<typename Field_T>
MPI_Datatype mpiDatatypeSliceXYZ( const Field_T & field, const CellInterval & interval, cell_idx_t f = 0 );



//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate an interval of cells of a Field
*
*  Only the cells in the described interval of a Field are communicated, when the returned datatype is used in MPI
*  communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field        The Field to be communicated
*  \param interval     Interval of communicated cells
*  \param fBeg         Beginning of the interval in f direction
*  \param fEnd         End of the interval in f direction (included in the interval)
*  \returns            The MPI datatype
*/
//======================================================================================================================
template<typename Field_T>
MPI_Datatype mpiDatatypeSliceXYZ( const Field_T & field, const CellInterval & interval, const cell_idx_t fBeg, const cell_idx_t fEnd );



//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate an interval of cells of a Field
*
*  Only the cells in the described interval of a Field are communicated, when the returned datatype is used in MPI
*  communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field        The Field to be communicated
*  \param interval     Interval of communicated cells
*  \param fs           Set of f components to be communicated
*  \returns            The MPI datatype
*/
//======================================================================================================================
template<typename Field_T>
MPI_Datatype mpiDatatypeSliceXYZ( const Field_T & field, const CellInterval & interval, const std::set<cell_idx_t> & fs );



//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate the contents of a whole GhostLayerField including its ghost layers
*
*  All cells of a Field are communicated, when the returned datatype is used in MPI communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field        The GhostLayerField to be communicated
*  \returns            The MPI datatype
*/
//======================================================================================================================
template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeWithGhostLayer( const GhostLayerField_T & field );



//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate the contents of a whole GhostLayerField including some of its ghost layers
*
*  All inner cells and the specified ghost layers of a Field are communicated, when the returned datatype is used in MPI
*  communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field          The GhostLayerField to be communicated
*  \param numGhostLayers The number of ghost layers to be communicated (from inner to outer layers)
*  \returns              The MPI datatype
*/
//======================================================================================================================
template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeWithGhostLayer( const GhostLayerField_T & field, const uint_t numGhostLayers );



//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate parts of the ghost layers of a GhostLayerField
*
*  All ghost layers in the specified direction of a Field are communicated, when the returned datatype is used in MPI
*  communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field          The GhostLayerField to be communicated
*  \param dir            As described at GhostLayerField::getGhostRegion
*  \param fullSlice      As described at GhostLayerField::getGhostRegion
*  \returns              The MPI datatype
*/
//======================================================================================================================
template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeGhostLayerOnly( const GhostLayerField_T & field, const stencil::Direction dir, const bool fullSlice = false );



//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate parts of the ghost layers of a GhostLayerField
*
*  The specified ghost layers of a Field are communicated, when the returned datatype is used in MPI communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field          The GhostLayerField to be communicated
*  \param thickness      As described at GhostLayerField::getGhostRegion
*  \param dir            As described at GhostLayerField::getGhostRegion
*  \param fullSlice      As described at GhostLayerField::getGhostRegion
*  \returns              The MPI datatype
*/
//======================================================================================================================
template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeGhostLayerOnly( const GhostLayerField_T & field, const uint_t thickness, const stencil::Direction dir, const bool fullSlice = false );



//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate parts of the ghost layers of a GhostLayerField
*
*  The specified ghost layers of a Field are communicated, when the returned datatype is used in MPI communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field          The GhostLayerField to be communicated
*  \param dir            As described at GhostLayerField::getGhostRegion
*  \param fullSlice      As described at GhostLayerField::getGhostRegion
*  \param f              Component which is communicated
*  \returns              The MPI datatype
*/
//======================================================================================================================
template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeGhostLayerOnlyXYZ( const GhostLayerField_T & field, const stencil::Direction dir, const bool fullSlice = false, const cell_idx_t f = 0 );


//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate parts of the ghost layers of a GhostLayerField
*
*  The specified ghost layers of a Field are communicated, when the returned datatype is used in MPI communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field          The GhostLayerField to be communicated
*  \param dir            As described at GhostLayerField::getGhostRegion
*  \param fullSlice      As described at GhostLayerField::getGhostRegion
*  \param fBeg           Beginning of the interval of communicated components
*  \param fEnd           End of the interval of communicated components (included in the interval)
*  \returns              The MPI datatype
*/
//======================================================================================================================
template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeGhostLayerOnlyXYZ( const GhostLayerField_T & field, const stencil::Direction dir, const bool fullSlice, const cell_idx_t fBeg, const cell_idx_t fEnd );


//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate parts of the ghost layers of a GhostLayerField
*
*  The specified ghost layers of a Field are communicated, when the returned datatype is used in MPI communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field          The GhostLayerField to be communicated
*  \param dir            As described at GhostLayerField::getGhostRegion
*  \param fullSlice      As described at GhostLayerField::getGhostRegion
*  \param fs             Communicated components
*  \returns              The MPI datatype
*/
//======================================================================================================================
template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeGhostLayerOnlyXYZ( const GhostLayerField_T & field, const stencil::Direction dir, const bool fullSlice, const std::set<cell_idx_t> & fs );


//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate inner parts if a GhostLayerField near the ghost layers
*
*  The specified inner layers of a Field are communicated, when the returned datatype is used in MPI communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field          The GhostLayerField to be communicated
*  \param dir            As described at GhostLayerField::getSliceBeforeGhostLayer
*  \param thickness      As described at GhostLayerField::getSliceBeforeGhostLayer
*  \param fullSlice      As described at GhostLayerField::getSliceBeforeGhostLayer
*  \returns              The MPI datatype
*/
//======================================================================================================================
template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeSliceBeforeGhostlayer( const GhostLayerField_T & field, const stencil::Direction dir, const uint_t thickness = 1, const bool fullSlice = false );


//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate inner parts if a GhostLayerField near the ghost layers
*
*  The specified inner layers of a Field are communicated, when the returned datatype is used in MPI communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field          The GhostLayerField to be communicated
*  \param dir            As described at GhostLayerField::getSliceBeforeGhostLayer
*  \param thickness      As described at GhostLayerField::getSliceBeforeGhostLayer
*  \param fullSlice      As described at GhostLayerField::getSliceBeforeGhostLayer
*  \param f              Component which is communicated
*  \returns              The MPI datatype
*/
//======================================================================================================================
template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeSliceBeforeGhostlayerXYZ( const GhostLayerField_T & field, const stencil::Direction dir, const uint_t thickness = 1, const cell_idx_t f = 0, const bool fullSlice = false );


//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate inner parts if a GhostLayerField near the ghost layers
*
*  The specified inner layers of a Field are communicated, when the returned datatype is used in MPI communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field          The GhostLayerField to be communicated
*  \param dir            As described at GhostLayerField::getSliceBeforeGhostLayer
*  \param thickness      As described at GhostLayerField::getSliceBeforeGhostLayer
*  \param fullSlice      As described at GhostLayerField::getSliceBeforeGhostLayer
*  \param fBeg           Beginning of the interval of communicated components
*  \param fEnd           End of the interval of communicated components (included in the interval)
*  \returns              The MPI datatype
*/
//======================================================================================================================
template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeSliceBeforeGhostlayerXYZ( const GhostLayerField_T & field, const stencil::Direction dir, const uint_t thickness, const cell_idx_t fBeg, const cell_idx_t fEnd, const bool fullSlice );


//======================================================================================================================
/*!
*  \brief Creates a MPI datatype to communicate inner parts if a GhostLayerField near the ghost layers
*
*  The specified inner layers of a Field are communicated, when the returned datatype is used in MPI communication.
*
*  The returned MPI_Datatype still has to be committed before used in communication and should be freed if it is not
*  needed any longer.
*
*  \param field          The GhostLayerField to be communicated
*  \param dir            As described at GhostLayerField::getSliceBeforeGhostLayer
*  \param thickness      As described at GhostLayerField::getSliceBeforeGhostLayer
*  \param fullSlice      As described at GhostLayerField::getSliceBeforeGhostLayer
*  \param fs             Communicated components
*  \returns              The MPI datatype
*/
//======================================================================================================================
template<typename GhostLayerField_T>
MPI_Datatype mpiDatatypeSliceBeforeGhostlayerXYZ( const GhostLayerField_T & field, const stencil::Direction dir, const uint_t thickness, const std::set<cell_idx_t> & fs, const bool fullSlice );


} // namespace communication
} // namespace field
} // namespace walberla

#include "MPIDatatypes.impl.h"
