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
//! \file ScalarFieldFromCellInterval.impl.h
//! \ingroup geometry
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//! \brief Implementations for ScalarFieldFromCellInterval
//
//======================================================================================================================



namespace walberla {
namespace geometry {
namespace initializer {



template< typename Field_T >
ScalarFieldFromCellInterval<Field_T>::ScalarFieldFromCellInterval( StructuredBlockStorage & blocks, std::vector<BlockDataID> fieldId )
   : structuredBlockStorage_( blocks ), scalarFieldID_ ( fieldId )
{}

template < typename Field_T>
void ScalarFieldFromCellInterval<Field_T>::init ( BlockStorage & /*blockStorage*/, const Config::BlockHandle & block )
{
	return init(block);
}

template< typename Field_T >
void ScalarFieldFromCellInterval<Field_T>::init( const Config::BlockHandle & blockHandle )
{
   if( !blockHandle )
      return;
   
   CellInterval globalCellInterval;
   if ( blockHandle.isDefined("CellInterval") )
      globalCellInterval = blockHandle.getParameter<CellInterval>("CellInterval");
   else
   {
      const Vector3<cell_idx_t> min = blockHandle.getParameter< Vector3<cell_idx_t> >( "min");
      const Vector3<cell_idx_t> max = blockHandle.getParameter< Vector3<cell_idx_t> >( "max");

      globalCellInterval.min() = Cell( min[0], min[1], min[2] );
      globalCellInterval.max() = Cell( max[0], max[1], max[2] );
   }
   
   std::vector<BlockDataID>::size_type id = blockHandle.getParameter<std::vector<BlockDataID>::size_type>("id");
   std::string expression                 = blockHandle.getParameter<std::string>                        ("value");
   
   try
   {
      Value_T value = string_to_num<Value_T>(expression);
      init(globalCellInterval, value, id);
   }
   
   catch(std::invalid_argument&)
   {
      math::FunctionParser p;
      p.parse(expression);
      init(globalCellInterval, p, id);
   }
}


template< typename Field_T >
void ScalarFieldFromCellInterval<Field_T>::init( const CellInterval & globalCellInterval, Value_T value, std::vector<BlockDataID>::size_type id )
{
   for( auto blockIt = structuredBlockStorage_.begin(); blockIt != structuredBlockStorage_.end(); ++blockIt )
   {
      auto f = blockIt->template getData<Field_T>(scalarFieldID_[id]);
      
      CellInterval localCellInterval;
      structuredBlockStorage_.transformGlobalToBlockLocalCellInterval( localCellInterval, *blockIt, globalCellInterval );
      localCellInterval.intersect( f->xyzSizeWithGhostLayer() );

      for (auto cell: localCellInterval)
         f->get( cell ) = value;
   }
}

template< typename Field_T >
void ScalarFieldFromCellInterval<Field_T>::init( const CellInterval & globalCellInterval, math::FunctionParser & parser, std::vector<BlockDataID>::size_type id )
{
   for( auto blockIt = structuredBlockStorage_.begin(); blockIt != structuredBlockStorage_.end(); ++blockIt )
   {
      auto f = blockIt->template getData<Field_T>(scalarFieldID_[id]);
      
      CellInterval localCellInterval;
      structuredBlockStorage_.transformGlobalToBlockLocalCellInterval( localCellInterval, *blockIt, globalCellInterval );
      localCellInterval.intersect( f->xyzSizeWithGhostLayer() );

      for (auto cell: localCellInterval)
      {
         std::map<std::string,Value_T> params;
         Cell global_cell;
         structuredBlockStorage_.transformBlockLocalToGlobalCell(global_cell, *blockIt, cell);
         params["x"] = global_cell.x();
         params["y"] = global_cell.y();
         params["z"] = global_cell.z();
         f->get( cell ) = parser.evaluate(params);
      }
   }
}



} // namespace initializer
} // namespace geometry
} // namespace walberla
