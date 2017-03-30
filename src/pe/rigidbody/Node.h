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
//! \file Node.h
//! \ingroup pe
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <pe/Types.h>
#include <core/DataTypes.h>

namespace walberla {
namespace pe {

//=================================================================================================
//
//  SPECIALIZATION FOR THE UNION FIND ALGORITHM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Node for the 'Union Find' algorithm.
 * \ingroup batch_generation
 *
 * This specialization of the BodyTrait class template adapts rigid bodies to the 'Union Find'
 * batch generation algorithm. In this algorithm, rigid bodies act as nodes in a graph and
 * contacts between the rigid bodies act as edges.
 */
class Node
{
public:

   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline Node();
   //@}
   //**********************************************************************************************

   //**Node functions******************************************************************************
   /*!\name Node functions */
   //@{
   inline ConstNodeID getRoot()   const;
   inline void        resetNode() const;
   //@}
   //**********************************************************************************************

public:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   mutable ConstNodeID root_;  //!< The root of the contact graph containing the rigid body.
   mutable size_t rank_;       //!< The current rank of the rigid body in the contact graph.
   mutable size_t batch_;      //!< Index of the batch containing all contacts in the contact graph.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default constructor for the BodyTrait<UnionFind> specialization.
 */
inline Node::Node()
   : rank_ ( 0 )
   , batch_( 0 )
{
   root_ = this ;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculating the current root of the contact graph.
 *
 * \return The current root of the contact graph.
 *
 * This function returns the current root of the contact graph and compresses the graph for
 * a fast access to the root node.
 */
inline ConstNodeID Node::getRoot() const
{
   ConstNodeID root( this );

   // Traverse the contact graph to find the root node
   while( root->root_ != root ) {
      root = root->root_;
   }

   // Retraverse the graph for graph compression
   ConstNodeID tmp( this );
   while( tmp->root_ != tmp ) {
      ConstNodeID last = tmp;
      tmp = tmp->root_;
      last->root_ = root;
   }

   return root;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the data members of the node.
 *
 * \return void
 */
inline void Node::resetNode() const
{
   root_  = this;
   rank_  = 0;
   batch_ = 0;
}
//*************************************************************************************************

} // namespace pe

} // namespace walberla
