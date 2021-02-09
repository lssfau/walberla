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
//! \file DumpBlockStructureLevel.h
//! \ingroup vtk
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "BlockCellDataWriter.h"


namespace walberla {
namespace vtk {



class DumpBlockStructureLevel : public ::walberla::vtk::BlockCellDataWriter< uint8_t > {

public:

   DumpBlockStructureLevel( const std::string& id ) : ::walberla::vtk::BlockCellDataWriter<uint8_t>( id ), level_( uint8_c(0) ) {}

protected:

   void configure() override { WALBERLA_ASSERT_NOT_NULLPTR( this->block_ ); WALBERLA_ASSERT_NOT_NULLPTR( this->blockStorage_ ); level_ = uint8_c( this->blockStorage_->getLevel( *(this->block_) ) ); }

   uint8_t evaluate( const cell_idx_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) override { return level_; }

   uint8_t level_;

}; // DumpBlockStructureLevel



} // namespace vtk
} // namespace walberla
