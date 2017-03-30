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
//! \file HilbertCurveConstruction.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "Types.h"


namespace walberla {
namespace blockforest {



// see: "Dynamic Octree Load Balancing Using Space-Filling Curves", Campbell, Devine, Gervasio, Teresco, 2003
static const uint_t hilbertOrder[24][8] = { { 0, 1, 3, 2, 6, 7, 5, 4 },
                                            { 0, 4, 6, 2, 3, 7, 5, 1 },
                                            { 0, 1, 5, 4, 6, 7, 3, 2 },
                                            { 5, 1, 0, 4, 6, 2, 3, 7 },
                                            { 3, 7, 6, 2, 0, 4, 5, 1 },
                                            { 6, 7, 3, 2, 0, 1, 5, 4 },
                                            { 5, 1, 3, 7, 6, 2, 0, 4 },
                                            { 0, 4, 5, 1, 3, 7, 6, 2 },
                                            { 5, 4, 0, 1, 3, 2, 6, 7 },
                                            { 5, 4, 6, 7, 3, 2, 0, 1 },
                                            { 0, 2, 3, 1, 5, 7, 6, 4 },
                                            { 6, 4, 0, 2, 3, 1, 5, 7 },
                                            { 5, 7, 3, 1, 0, 2, 6, 4 },
                                            { 3, 7, 5, 1, 0, 4, 6, 2 },
                                            { 6, 4, 5, 7, 3, 1, 0, 2 },
                                            { 0, 2, 6, 4, 5, 7, 3, 1 },
                                            { 6, 2, 0, 4, 5, 1, 3, 7 },
                                            { 6, 2, 3, 7, 5, 1, 0, 4 },
                                            { 3, 2, 0, 1, 5, 4, 6, 7 },
                                            { 6, 7, 5, 4, 0, 1, 3, 2 },
                                            { 5, 7, 6, 4, 0, 2, 3, 1 },
                                            { 3, 2, 6, 7, 5, 4, 0, 1 },
                                            { 3, 1, 0, 2, 6, 4, 5, 7 },
                                            { 3, 1, 5, 7, 6, 4, 0, 2 } };

static const uint_t hilbertOrientation[24][8] = { {  1,  2,  0,  3,  4,  0,  5,  6 },
                                                  {  0,  7,  1,  8,  5,  1,  4,  9 },
                                                  { 15,  0,  2, 22, 20,  2, 19, 23 },
                                                  { 20,  6,  3, 23, 15,  3, 16, 22 },
                                                  { 22, 13,  4, 12, 11,  4,  1, 20 },
                                                  { 11, 19,  5, 20, 22,  5,  0, 12 },
                                                  {  9,  3,  6,  2, 21,  6, 17,  0 },
                                                  { 10,  1,  7, 11, 12,  7, 13, 14 },
                                                  { 12,  9,  8, 14, 10,  8, 18, 11 },
                                                  {  6,  8,  9,  7, 17,  9, 21,  1 },
                                                  {  7, 15, 10, 16, 13, 10, 12, 17 },
                                                  {  5, 14, 11,  9,  0, 11, 22,  8 },
                                                  {  8, 20, 12, 19, 18, 12, 10,  5 },
                                                  { 18,  4, 13,  5,  8, 13,  7, 19 },
                                                  { 17, 11, 14,  1,  6, 14, 23,  7 },
                                                  {  2, 10, 15, 18, 19, 15, 20, 21 },
                                                  { 19, 17, 16, 21,  2, 16,  3, 18 },
                                                  { 14, 16, 17, 15, 23, 17,  6, 10 },
                                                  { 13, 21, 18, 17,  7, 18,  8, 16 },
                                                  { 16,  5, 19,  4,  3, 19,  2, 13 },
                                                  {  3, 12, 20, 13, 16, 20, 15,  4 },
                                                  { 23, 18, 21, 10, 14, 21,  9, 15 },
                                                  {  4, 23, 22,  6,  1, 22, 11,  3 },
                                                  { 21, 22, 23,  0,  9, 23, 14,  2 } };



} // namespace blockforest
} // namespace walberla


