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
//! \file VTKWriter.h
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "field/Field.h"
#include "core/VectorTrait.h"
#include "core/debug/Debug.h"
#include "vtk/BlockCellDataWriter.h"
#include <vtk/VTKOutput.h>


namespace walberla {
namespace field {


//**********************************************************************************************************************
/*! Write out a field in VTK representation (for Paraview )
*
*  Template parameters:
*      Field_T:    - can be a Field, GhostLayerField or FieldAdaptor (or anything that implements this concept )
*                  - Field_T::value_type has to be a native datatype, or a Vector3<NativeDataType>
*      OutputType: - Before writing out the data, it is converted to this type
*                  - Typical use case: instead of writing double fields, write out float, so the output files are smaller
*                  - This type defaults to Field_T::value_type or in case of vector fields,
*                    the type of the Vector3 (i.e. Field_T::value_type::value_type )
*                  - Examples: Field_T = GhostLayerField<double,3>
*                                  -> By default OutputType would be double,
*                                     however can be set to float to get smaller output
*                              Field_T = GhostLayerField<Vector3<double>,1>
*                                  -> Default would be double, not Vector3<double>!
*
* \param fieldId              BlockDataID for the field that should be written out
* \param blocks               StructuredBlockStorage where the field is stored
* \param identifier           name of the field in the paraview output (has to be unique for all outputs)
* \param writeFrequency       write output only every  writeFrequency'th timestep (call of returned functor)
* \param ghostLayers          number of ghost layers to include in the output
* \param forcePVTU            force unstructured output ( currently required if not all blocks have the field allocated )
* \param baseFolder           name of the base folder
* \param executionFolder      name for the folder that stores all files corresponding to one step
* \param continuousNumbering  important when writeFrequency > 1,  if true the simulation steps are numbered continuously
*                             otherwise they are numbered after the timestep
* \param binary               switch for binary output
* \param litteEndian          byte order
* \param simultaneousIOOps    number of simultaneous IO operations, 0 means all processes write concurrently
*                             limiting the number of simultaneous IO operations makes sense for huge number of processes
*                             in order to not overload the file system
* \param requiredStates       see selectable concept
* \param incompatibleStates   see selectable concept
*
* \return  Returns a functor. Every writeFrequency'th call of this functor, the field is written out.
*          This functor usually is added to timeloop.
*/
//**********************************************************************************************************************

// template default argument emulated using overload:  OutputType = internal::DefaultOutputType<Field_T>::type
template< typename Field_T, typename OutputType  >
inline vtk::VTKOutput::Write createVTKOutput  ( const ConstBlockDataID & fieldId, const StructuredBlockStorage & blocks,
                                                const std::string & identifier,
                                                const uint_t writeFrequency = 1,
                                                const uint_t ghostLayers = 0,
                                                const bool forcePVTU = false,
                                                const std::string & baseFolder = std::string("vtk_out"),
                                                const std::string & executionFolder = std::string("simulation_step"),
                                                const bool continuousNumbering = false,
                                                const bool binary = true,
                                                const bool littleEndian = true,
                                                const int simultaneousIOOperations = 0,
                                                const Set<SUID>& requiredStates     = Set<SUID>::emptySet(),
                                                const Set<SUID>& incompatibleStates = Set<SUID>::emptySet(),
                                                bool useMPIIO = true, const uint_t initialExecutionCount = 0 );




//======================================================================================================================
//
//  VTKWriter
//
//======================================================================================================================

template< typename Field_T, typename OutputType = typename VectorTrait<typename Field_T::value_type >::OutputType >
class VTKWriter : public vtk::BlockCellDataWriter< OutputType,
                                                   VectorTrait<typename Field_T::value_type>::F_SIZE * Field_T::F_SIZE >
{
public:
   using OutputTrait = VectorTrait<typename Field_T::value_type>;

   using base_t = vtk::BlockCellDataWriter<OutputType, OutputTrait::F_SIZE * Field_T::F_SIZE>;

   VTKWriter( const ConstBlockDataID bdid, const std::string& id ) :
      base_t( id ), bdid_( bdid ), field_( nullptr ) {}

protected:

   void configure() override {
      WALBERLA_ASSERT_NOT_NULLPTR( this->block_ );
      field_ = this->block_->template getData< Field_T >( bdid_ );
   }

   OutputType evaluate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f ) override
   {
      WALBERLA_ASSERT_NOT_NULLPTR( field_ );

      if ( Field_T::F_SIZE == 1 )
      {
         // common case that can be handled faster
         return numeric_cast<OutputType>( OutputTrait::get( field_->get(x,y,z,0), uint_c(f) )  );
      }
      else
      {
         const cell_idx_t fField = f / cell_idx_c( OutputTrait::F_SIZE );
         const cell_idx_t fType  = f % cell_idx_c( OutputTrait::F_SIZE );

         return numeric_cast<OutputType>( OutputTrait::get( field_->get(x,y,z,fField), uint_c(fType) )  );
      }
   }

   const ConstBlockDataID bdid_;
   const Field_T * field_;

};




//======================================================================================================================
//
//  Implementation convenience functions - for documentation see above
//
//======================================================================================================================


template< typename Field_T, typename OutputType >
inline vtk::VTKOutput::Write createVTKOutput  ( const ConstBlockDataID & fieldId, const StructuredBlockStorage & blocks,
                                                const std::string & identifier,
                                                const uint_t writeFrequency,
                                                const uint_t ghostLayers,
                                                const bool forcePVTU,
                                                const std::string & baseFolder,
                                                const std::string & executionFolder,
                                                const bool continuousNumbering,
                                                const bool binary,
                                                const bool littleEndian,
                                                const int simultaneousIOOperations,
                                                const Set<SUID>& requiredStates,
                                                const Set<SUID>& incompatibleStates,
                                                bool useMPIIO, const uint_t initialExecutionCount )
{
   auto vtkOutput = vtk::createVTKOutput_BlockData( blocks, identifier, writeFrequency, ghostLayers, forcePVTU, baseFolder, executionFolder,
                                                    continuousNumbering, binary, littleEndian, useMPIIO, initialExecutionCount );

   vtkOutput->addCellDataWriter( make_shared< VTKWriter< Field_T,OutputType > >( fieldId, identifier ) );

   return writeFiles( vtkOutput, true, simultaneousIOOperations, requiredStates, incompatibleStates );
}


// mimic template default argument OutputType = internal::DefaultOutputType<Field_T>::type
template< typename Field_T  >
inline vtk::VTKOutput::Write createVTKOutput  ( const ConstBlockDataID & fieldId, const StructuredBlockStorage & blocks,
                                                const std::string & identifier,
                                                const uint_t writeFrequency = 1,
                                                const uint_t ghostLayers = 0,
                                                const bool forcePVTU = false,
                                                const std::string & baseFolder = std::string("vtk_out"),
                                                const std::string & executionFolder = std::string("simulation_step"),
                                                const bool continuousNumbering = false,
                                                const bool binary = true,
                                                const bool littleEndian = true,
                                                const int simultaneousIOOperations = 0,
                                                const Set<SUID>& requiredStates     = Set<SUID>::emptySet(),
                                                const Set<SUID>& incompatibleStates = Set<SUID>::emptySet(),
                                                bool useMPIIO = true, const uint_t initialExecutionCount = 0 )
{
   using OutputType = typename VectorTrait<typename Field_T::value_type>::OutputType;

   return createVTKOutput<Field_T, OutputType> (
               fieldId, blocks, identifier, writeFrequency, ghostLayers, forcePVTU, baseFolder,executionFolder,
               continuousNumbering, binary, littleEndian, simultaneousIOOperations, requiredStates, incompatibleStates,
               useMPIIO, initialExecutionCount
            );
}


template< typename Field_T, typename OutputType >
inline vtk::VTKOutput::Write createScalingVTKOutput( const ConstBlockDataID & fieldId, const StructuredBlockStorage & blocks,
                                                     const std::string & identifier,
                                                     const uint_t writeFrequency,
                                                     const typename VectorTrait<typename Field_T::value_type >::OutputType factor,
                                                     const uint_t ghostLayers = 0,
                                                     const bool forcePVTU = false,
                                                     const std::string & baseFolder = std::string( "vtk_out" ),
                                                     const std::string & executionFolder = std::string( "simulation_step" ),
                                                     const bool continuousNumbering = false,
                                                     const bool binary = true,
                                                     const bool littleEndian = true,
                                                     const int simultaneousIOOperations = 0,
                                                     const Set<SUID>& requiredStates = Set<SUID>::emptySet(),
                                                     const Set<SUID>& incompatibleStates = Set<SUID>::emptySet(),
                                                     bool useMPIIO = true, const uint_t initialExecutionCount = 0 )
{
   auto vtkOutput = vtk::createVTKOutput_BlockData( blocks, identifier, writeFrequency, ghostLayers, forcePVTU, baseFolder, executionFolder,
                                                    continuousNumbering, binary, littleEndian, useMPIIO, initialExecutionCount );

   auto writer = make_shared< VTKWriter< Field_T,OutputType > >( fieldId, identifier );
   vtkOutput->addCellDataWriter( vtk::makeBlockCellDataWriterScalingAdapter(identifier, writer, factor) );

   return writeFiles( vtkOutput, true, simultaneousIOOperations, requiredStates, incompatibleStates );
}


// mimic template default argument OutputType = internal::DefaultOutputType<Field_T>::type
template< typename Field_T  >
inline vtk::VTKOutput::Write createScalingVTKOutput  ( const ConstBlockDataID & fieldId, const StructuredBlockStorage & blocks,
                                                       const std::string & identifier,
                                                       const uint_t writeFrequency,
                                                       const typename VectorTrait<typename Field_T::value_type >::OutputType factor,
                                                       const uint_t ghostLayers = 0,
                                                       const bool forcePVTU = false,
                                                       const std::string & baseFolder = std::string("vtk_out"),
                                                       const std::string & executionFolder = std::string("simulation_step"),
                                                       const bool continuousNumbering = false,
                                                       const bool binary = true,
                                                       const bool littleEndian = true,
                                                       const int simultaneousIOOperations = 0,
                                                       const Set<SUID>& requiredStates     = Set<SUID>::emptySet(),
                                                       const Set<SUID>& incompatibleStates = Set<SUID>::emptySet(),
                                                       bool useMPIIO = true, const uint_t initialExecutionCount = 0 )
{
   using OutputType = typename VectorTrait<typename Field_T::value_type>::OutputType;

   return createScalingVTKOutput<Field_T, OutputType> (
               fieldId, blocks, identifier, writeFrequency, factor, ghostLayers, forcePVTU, baseFolder,executionFolder,
               continuousNumbering, binary, littleEndian, simultaneousIOOperations, requiredStates, incompatibleStates,
               useMPIIO, initialExecutionCount
            );
}



} // namespace field
} // namespace walberla
