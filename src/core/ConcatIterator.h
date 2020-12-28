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
//! \file ConcatIterator.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//======================================================================================================================

#pragma once

#include <iterator>

namespace walberla {

template< typename T >
class ConcatIterator
{
public:
    using iterator_category = std::input_iterator_tag;
    using value_type = typename T::value_type;
    using difference_type = typename T::difference_type;
    using pointer = typename T::pointer;
    using reference = typename T::reference;

    /// default constructed iterator points to end()
    ConcatIterator( )
        : ended_( true )
    {}

    ConcatIterator( const T& beginFirst,
                    const T& endFirst,
                    const T& beginSecond,
                    const T& endSecond )
        : it_(beginFirst)
        , itEndFirst_(endFirst)
        , itBeginSecond_(beginSecond)
        , itEndSecond_(endSecond)
        , first_(true)
        , ended_(false)
    {
        checkStateAndAdapt();
    }

    ConcatIterator& operator++()    { ++it_; checkStateAndAdapt(); return *this; }            // prefix ++X
    ConcatIterator  operator++(int) { ConcatIterator it( *this ); operator++(); return it; }; // postfix X++

    bool operator==( const ConcatIterator& rhs ) const
    {
        if (ended_ || rhs.ended_)
        {
           if (ended_ == rhs.ended_)
           {
              return true;
           }
           else
           {
              return false;
           }
        }

        return it_ == rhs.it_;
    }
    bool operator!=( const ConcatIterator & rhs ) const { return !(*this == rhs); }

    typename T::reference  operator*()  { return it_.operator*(); }
    typename T::pointer    operator->() { return it_.operator->(); }

private:

    /// checks if iterator has to jump to second range
    void checkStateAndAdapt()
    {
        if( first_ && it_ == itEndFirst_ )
        {
           it_ = itBeginSecond_;
           first_ = false;
        }

        if( !first_ && it_ == itEndSecond_ )
        {
            ended_ = true;
        }
    }

    /// current position
    T it_;
    /// end of first range
    const T itEndFirst_;
    /// begin of second range
    const T itBeginSecond_;
    /// end of second range
    const T itEndSecond_;

    ///located in first range?
    bool first_;

    /// end reached?
    bool ended_;
};

} // namespace walberla
