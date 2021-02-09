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
//! \file UniformMPIDatatypeInfo.h
//! \ingroup communication
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/mpi/Datatype.h"
#include "domain_decomposition/IBlock.h"
#include "stencil/Directions.h"

namespace walberla {
namespace communication {


   //*******************************************************************************************************************
   /*!
   * Interface for direction MPI datatype based communication using blockforest::UniformDirectScheme
   *
   * \ingroup communication
   *
   * For direct bufferless MPI communication the MPI data type concept has to be used.
   * With this interface a MPI data type can be associated with block data, such that this block data can be
   * communicated using e.g. the blockforest::communication::UniformDirectScheme.
   * In addition to providing the MPI data types for sending and receiving in each direction, a pointer
   * to the data has to be provided.
   */
   //*******************************************************************************************************************
   class UniformMPIDatatypeInfo
   {
   public:

      //**Construction & Destruction************************************************************************************
      /*! \name Construction & Destruction */
      //@{
               UniformMPIDatatypeInfo() = default;
      virtual ~UniformMPIDatatypeInfo() = default;
      //@}
      //****************************************************************************************************************

      /*************************************************************************************************************//**
      * Return the MPI data type that should be used for sending to neighbor in  specified direction
      *****************************************************************************************************************/
      virtual shared_ptr<mpi::Datatype> getSendDatatype ( IBlock * block, const stencil::Direction dir ) = 0;

      /*************************************************************************************************************//**
      * Return the MPI data type that should be used for receiving from neighbor in specified direction
      *****************************************************************************************************************/
      virtual shared_ptr<mpi::Datatype> getRecvDatatype ( IBlock * block, const stencil::Direction dir ) = 0;




      /*************************************************************************************************************//**
      * Return pointer to data that should be send to neighbor in specified direction
      *****************************************************************************************************************/
      virtual void * getSendPointer( IBlock * block, const stencil::Direction dir ) = 0;


      /*************************************************************************************************************//**
      * Return pointer to memory where received data is written to
      *****************************************************************************************************************/
      virtual void * getRecvPointer( IBlock * block, const stencil::Direction dir ) = 0;



      /*************************************************************************************************************//**
      * Return how many data items of the datatype should be communicated per block and direction
      * Due to custom aggregated MPI datatypes this is usually 1
      *****************************************************************************************************************/
      virtual int getNumberOfItemsToCommunicate( IBlock * , const stencil::Direction ) { return 1; }
   };




} // namespace communication
} // namespace walberla


