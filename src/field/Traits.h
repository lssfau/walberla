
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
//! \file Traits.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once


namespace walberla {
namespace field {




   template<typename Pair>
   struct ToField
   {
       typedef Field< typename Pair::first, Pair::second::value > type;
   };

   template<typename TupleList>
   struct ToFieldList
   {
      typedef mpl::transform< TupleList, ToField<mpl::_1> >::type type;
   };

   template<typename Pair>
   struct ToGhostLayerField
   {
      typedef GhostLayerField< typename Pair::first, Pair::second::value > type;
   };

   template<typename TupleList>
   struct ToGhostLayerFieldList
   {
      typedef mpl::transform< TupleList, ToGhostLayerField<mpl::_1> >::type type;
   };


} // namespace field
} // namespace walberla



