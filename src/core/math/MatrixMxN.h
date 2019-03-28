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
//! \file MatrixMxN.h
//! \ingroup core
//! \author Klaus Iglberger
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <stdexcept>
#include <core/math/Utility.h>
#include <core/math/MathTrait.h>
#include <core/debug/Debug.h>
#include <core/DataTypes.h>
#include <core/Macros.h>
#include <core/math/Shims.h>

#include <type_traits>


namespace walberla {
namespace math {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Efficient implementation of a \f$ M \times N \f$ matrix.
 *
 * The MatrixMxN class is the representation of a dynamic \f$ M \times N \f$ matrix with a total
 * of \f$ M \cdot N \f$ dynamically allocated elements. These elements can be directly accessed
 * with the 1D subscript operator or with the 2D function operator. The matrix is stored in a
 * row-wise fashion:

               \f[\left(\begin{array}{*{5}{c}}
               0           & 1             & 2             & \cdots & N-1         \\
               N           & N+1           & N+2           & \cdots & 2 \cdot N-1 \\
               \vdots      & \vdots        & \vdots        & \ddots & \vdots      \\
               M \cdot N-N & M \cdot N-N+1 & M \cdot N-N+2 & \cdots & M \cdot N-1 \\
               \end{array}\right)\f]

 * MatrixMxN can be used with any non-cv-qualified element type. The arithmetic operators for
 * matrix/matrix, matrix/vector and matrix/element operations with the same element type work
 * for any element type as long as the element type supports the arithmetic operation. Arithmetic
 * operations between matrices, vectors and elements of different element types are only supported
 * for all data types supported by the MathTrait class template (for details see the MathTrait
 * class description).

   \code
   MatrixMxN< double > a, b, c;
   MatrixMxN< float  > d;
   MatrixMxN< std::complex<double> > e, f, g;
   MatrixMxN< std::complex<float>  > h;

   ...         // Appropriate resizing

   c = a + b;  // OK: Same element type, supported
   c = a + d;  // OK: Different element types, supported by the MathTrait class template

   g = e + f;  // OK: Same element type, supported
   g = e + h;  // Error: Different element types, not supported by the MathTrait class template
   \endcode
 */
template< typename Type >  // Data type of the matrix
class MatrixMxN
{
   //**Compile time checks*************************************************************************
   /*! \cond internal */
   static_assert(!std::is_const<Type>::value, "only non const Types are allowed!");
   static_assert(!std::is_volatile<Type>::value, "only non volatile types are allowed!");
   /*! \endcond */
   //**********************************************************************************************

public:
   //**Type definitions****************************************************************************
   typedef MatrixMxN<Type>   This;           //!< Type of this MatrixMxN instance.
   typedef This              ResultType;     //!< Result type for expression template evaluations.
   typedef Type              ElementType;    //!< Type of the matrix elements.
   typedef const MatrixMxN&  CompositeType;  //!< Data type for composite expression templates.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
                           explicit inline MatrixMxN();
                           explicit inline MatrixMxN( size_t m, size_t n );
                           explicit inline MatrixMxN( size_t m, size_t n, Type init );
                                    inline MatrixMxN( const MatrixMxN& m );

   template< typename Other, size_t M, size_t N >
   inline MatrixMxN( const Other (&rhs)[M][N] );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~MatrixMxN();
   //@}
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   template< typename Other, size_t M, size_t N >
   inline MatrixMxN& operator=( const Other (&rhs)[M][N] );

                           inline MatrixMxN&  operator= ( Type set );
                           inline MatrixMxN&  operator= ( const MatrixMxN& set );
                           inline Type&       operator[]( size_t index );
                           inline const Type& operator[]( size_t index )       const;
                           inline Type&       operator()( size_t i, size_t j );
                           inline const Type& operator()( size_t i, size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
                              inline size_t          rows()               const;
                              inline size_t          columns()            const;
                              inline size_t          capacity()           const;
                              inline size_t          nonZeros()           const;
                              inline size_t          nonZeros( size_t i ) const;
                              inline void            reset();
                              inline void            clear();
                                     void            resize ( size_t m, size_t n, bool preserve=true );
                              inline void            extend ( size_t m, size_t n, bool preserve=true );
                              inline void            reserve( size_t elements );
                              inline MatrixMxN&      transpose();
                              inline bool            isDiagonal()         const;
                              inline bool            isSymmetric()        const;
   template< typename Other > inline MatrixMxN&      scale( Other scalar );
                              inline void            swap( MatrixMxN& m ) /* throw() */;
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool isAliased ( const Other* alias ) const;
   //@}
   //**********************************************************************************************

private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;             //!< The current number of rows of the matrix.
   size_t n_;             //!< The current number of columns of the matrix.
   size_t capacity_;      //!< The maximum capacity of the matrix.
   Type* WALBERLA_RESTRICT v_;  //!< The dynamically allocated matrix elements.
                          /*!< Access to the matrix elements is gained via the subscript or
                               function call operator. The order of the elements is
                               \f[\left(\begin{array}{*{5}{c}}
                               0            & 1             & 2             & \cdots & N-1         \\
                               N            & N+1           & N+2           & \cdots & 2 \cdot N-1 \\
                               \vdots       & \vdots        & \vdots        & \ddots & \vdots      \\
                               M \cdot N-N  & M \cdot N-N+1 & M \cdot N-N+2 & \cdots & M \cdot N-1 \\
                               \end{array}\right)\f] */
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor for MatrixMxN.
 */
template< typename Type >  // Data type of the matrix
inline MatrixMxN<Type>::MatrixMxN()
   : m_       ( 0 )  // The current number of rows of the matrix
   , n_       ( 0 )  // The current number of columns of the matrix
   , capacity_( 0 )  // The maximum capacity of the matrix
   , v_       ( 0 )  // The matrix elements
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a matrix of size \f$ m \times n \f$. No element initialization is performed!
 *
 * \param m The number of rows of the matrix.
 * \param n The number of columns of the matrix.
 *
 * \b Note: This constructor is only responsible to allocate the required dynamic memory. No
 *          element initialization is performed!
 */
template< typename Type >  // Data type of the matrix
inline MatrixMxN<Type>::MatrixMxN( size_t m, size_t n )
   : m_       ( m )                    // The current number of rows of the matrix
   , n_       ( n )                    // The current number of columns of the matrix
   , capacity_( m*n )                  // The maximum capacity of the matrix
   , v_       ( new Type[capacity_] )  // The matrix elements
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a homogenous initialization of all \f$ m \times n \f$ matrix elements.
 *
 * \param m The number of rows of the matrix.
 * \param n The number of columns of the matrix.
 * \param init The initial value of the matrix elements.
 *
 * All matrix elements are initialized with the specified value.
 */
template< typename Type >  // Data type of the matrix
inline MatrixMxN<Type>::MatrixMxN( size_t m, size_t n, Type init )
   : m_       ( m )                    // The current number of rows of the matrix
   , n_       ( n )                    // The current number of columns of the matrix
   , capacity_( m*n )                  // The maximum capacity of the matrix
   , v_       ( new Type[capacity_] )  // The matrix elements
{
   for( size_t i=0; i<capacity_; ++i )
      v_[i] = init;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for MatrixMxN.
 *
 * \param m Matrix to be copied.
 *
 * The copy constructor is explicitly defined due to the required dynamic memory management
 * and in order to enable/facilitate NRV optimization.
 */
template< typename Type >  // Data type of the matrix
inline MatrixMxN<Type>::MatrixMxN( const MatrixMxN& m )
   : m_       ( m.m_  )                // The current number of rows of the matrix
   , n_       ( m.n_  )                // The current number of columns of the matrix
   , capacity_( m_*n_ )                // The maximum capacity of the matrix
   , v_       ( new Type[capacity_] )  // The matrix elements
{
   for( size_t i=0; i<capacity_; ++i )
      v_[i] = m.v_[i];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array initialization of all matrix elements.
 *
 * \param rhs \f$ M \times N \f$ dimensional array for the initialization.
 * \return Reference to the assigned matrix.
 *
 * This constructor offers the option to directly initialize the elements of the matrix:

   \code
   const real init[3][3] = { { 1, 2, 3 },
                             { 4, 5 },
                             { 7, 8, 9 } };
   MatrixMxN<real> A = init;
   \endcode

 * The matrix is sized accoring to the size of the array and initialized with the given values.
 * Missing values are initialized with zero (as e.g. the value 6 in the example).
 */
template< typename Type >  // Data type of the matrix
template< typename Other   // Data type of the initialization array
        , size_t M         // Number of rows of the initialization array
        , size_t N >       // Number of columns of the initialization array
inline MatrixMxN<Type>::MatrixMxN( const Other (&rhs)[M][N] )
   : m_       ( M )              // The current number of rows of the matrix
   , n_       ( N )              // The current number of columns of the matrix
   , capacity_( M*N )            // The maximum capacity of the matrix
   , v_       ( new Type[M*N] )  // The matrix elements
{
   for( size_t i=0; i<M; ++i )
      for( size_t j=0; j<N; ++j )
         v_[i*N+j] = rhs[i][j];
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for MatrixMxN.
 */
template< typename Type >  // Data type of the matrix
inline MatrixMxN<Type>::~MatrixMxN()
{
   delete [] v_;
}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Array assignment to all matrix elements.
 *
 * \param rhs \f$ M \times N \f$ dimensional array for the assignment.
 * \return Reference to the assigned matrix.
 *
 * This assignment operator offers the option to directly set all elements of the matrix:

   \code
   const real init[3][3] = { { 1, 2, 3 },
                             { 4, 5 },
                             { 7, 8, 9 } };
   MatrixMxN<real> A;
   A = init;
   \endcode

 * The matrix is resized accoring to the size of the array and initialized with the given values.
 * Missing values are initialized with zero (as e.g. the value 6 in the example).
 */
template< typename Type >  // Data type of the matrix
template< typename Other   // Data type of the initialization array
        , size_t M         // Number of rows of the initialization array
        , size_t N >       // Number of columns of the initialization array
inline MatrixMxN<Type>& MatrixMxN<Type>::operator=( const Other (&rhs)[M][N] )
{
   resize( M, N, false );

   for( size_t i=0; i<M; ++i )
      for( size_t j=0; j<N; ++j )
         v_[i*N+j] = rhs[i][j];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Homogenous assignment to all matrix elements.
 *
 * \param rhs Scalar value to be assigned to all matrix elements.
 * \return Reference to the assigned matrix.
 */
template< typename Type >  // Data type of the matrix
inline MatrixMxN<Type>& MatrixMxN<Type>::operator=( Type rhs )
{
   const size_t sqrsize( m_*n_ );
   for( size_t i=0; i<sqrsize; ++i )
      v_[i] = rhs;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for MatrixMxN.
 *
 * \param rhs Matrix to be copied.
 * \return Reference to the assigned matrix.
 *
 * The matrix is resized according to the given \f$ M \times N \f$ matrix and initialized as a
 * copy of this matrix.
 */
template< typename Type >  // Data type of the matrix
inline MatrixMxN<Type>& MatrixMxN<Type>::operator=( const MatrixMxN& rhs )
{
   if( &rhs == this ) return *this;

   resize( rhs.m_, rhs.n_, false );

   const size_t sqrsize( m_*n_ );
   for( size_t i=0; i<sqrsize; ++i )
      v_[i] = rhs.v_[i];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 1D-access to the matrix elements.
 *
 * \param index Access index. The index has to be in the range \f$[0..M \cdot N-1]\f$.
 * \return Reference to the accessed value.
 */
template< typename Type >  // Data type of the matrix
inline Type& MatrixMxN<Type>::operator[]( size_t index )
{
   WALBERLA_ASSERT( index < m_*n_, "Invalid matrix access index" );
   return v_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 1D-access to the matrix elements.
 *
 * \param index Access index. The index has to be in the range \f$[0..M \cdot N-1]\f$.
 * \return Reference to the accessed value.
 */
template< typename Type >  // Data type of the matrix
inline const Type& MatrixMxN<Type>::operator[]( size_t index ) const
{
   WALBERLA_ASSERT( index < m_*n_, "Invalid matrix access index" );
   return v_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 2D-access to the matrix elements.
 *
 * \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
 * \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
 * \return Reference to the accessed value.
 */
template< typename Type >  // Data type of the matrix
inline Type& MatrixMxN<Type>::operator()( size_t i, size_t j )
{
   WALBERLA_ASSERT( i<m_ && j<n_, "Invalid matrix access index" );
   return v_[i*n_+j];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 2D-access to the matrix elements.
 *
 * \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
 * \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
 * \return Reference to the accessed value.
 */
template< typename Type >  // Data type of the matrix
inline const Type& MatrixMxN<Type>::operator()( size_t i, size_t j ) const
{
   WALBERLA_ASSERT( i<m_ && j<n_, "Invalid matrix access index" );
   return v_[i*n_+j];
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current number of rows of the matrix.
 *
 * \return The number of rows of the matrix.
 */
template< typename Type >  // Data type of the matrix
inline size_t MatrixMxN<Type>::rows() const
{
   return m_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of columns of the matrix.
 *
 * \return The number of columns of the matrix.
 */
template< typename Type >  // Data type of the matrix
inline size_t MatrixMxN<Type>::columns() const
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the matrix.
 *
 * \return The capacity of the matrix.
 */
template< typename Type >  // Data type of the matrix
inline size_t MatrixMxN<Type>::capacity() const
{
   return capacity_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the total number of non-zero elements in the matrix
 *
 * \return The number of non-zero elements in the sparse matrix.
 */
template< typename Type >  // Data type of the matrix
inline size_t MatrixMxN<Type>::nonZeros() const
{
   size_t nonzeros( 0 );

   const size_t sqrsize( m_*n_ );
   for( size_t i=0; i<sqrsize; ++i )
      if( v_[i] != Type() )
         ++nonzeros;

   return nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the specified row.
 *
 * \param i The index of the row.
 * \return The number of non-zero elements of row \a i.
 */
template< typename Type >  // Data type of the matrix
inline size_t MatrixMxN<Type>::nonZeros( size_t i ) const
{
   WALBERLA_ASSERT( i < rows(), "Invalid row access index" );

   const size_t end( (i+1)*n_ );
   size_t nonzeros( 0 );

   for( size_t j=i*n_; j<end; ++j )
      if( v_[j] != Type() )
         ++nonzeros;

   return nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
 *
 * \return void
 */
template< typename Type >  // Data type of the matrix
inline void MatrixMxN<Type>::reset()
{
   const size_t sqrsize( m_*n_ );
   for( size_t i=0; i<sqrsize; ++i )
      math::reset(v_[i]);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the \f$ M \times N \f$ matrix.
 *
 * \return void
 *
 * After the clear() function, the size of the matrix is 0.
 */
template< typename Type >  // Data type of the matrix
inline void MatrixMxN<Type>::clear()
{
   m_ = 0;
   n_ = 0;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the matrix.
 *
 * \param m The new number of rows of the matrix.
 * \param n The new number of columns of the matrix.
 * \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
 * \return void
 *
 * This function resizes the matrix using the given size to \f$ m \times n \f$. During this
 * operation, new dynamic memory may be allocated in case the capacity of the matrix is too
 * small. Therefore this function potentially changes all matrix elements. In order to preserve
 * the old matrix values, the \a preserve flag can be set to \a true. However, new matrix
 * elements are not initialized!\n
 * The following example illustrates the resize operation of a \f$ 2 \times 4 \f$ matrix to a
 * \f$ 4 \times 2 \f$ matrix. The new, uninitialized elements are marked with \a x:

                              \f[
                              \left(\begin{array}{*{4}{c}}
                              1 & 2 & 3 & 4 \\
                              5 & 6 & 7 & 8 \\
                              \end{array}\right)

                              \Longrightarrow

                              \left(\begin{array}{*{2}{c}}
                              1 & 2 \\
                              5 & 6 \\
                              x & x \\
                              x & x \\
                              \end{array}\right)
                              \f]
 */
template< typename Type >  // Data type of the matrix
void MatrixMxN<Type>::resize( size_t m, size_t n, bool preserve )
{
   if( m == m_ && n == n_ ) return;

   if( preserve )
   {
      Type* WALBERLA_RESTRICT v = new Type[m*n];
      const size_t min_m( min( m, m_ ) );
      const size_t min_n( min( n, n_ ) );

      for( size_t i=0; i<min_m; ++i )
         for( size_t j=0; j<min_n; ++j )
            v[i*n+j] = v_[i*n_+j];

      std::swap( v_, v );
      delete [] v;
      capacity_ = m*n;
   }
   else if( m*n > capacity_ ) {
      Type* WALBERLA_RESTRICT v = new Type[m*n];
      std::swap( v_, v );
      delete [] v;
      capacity_ = m*n;
   }

   m_ = m;
   n_ = n;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Extending the size of the matrix.
 *
 * \param m Number of additional rows.
 * \param n Number of additional columns.
 * \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
 * \return void
 *
 * This function increases the matrix size by \a m rows and \a n columns. During this operation,
 * new dynamic memory may be allocated in case the capacity of the matrix is too small. Therefore
 * this function potentially changes all matrix elements. In order to preserve the old matrix
 * values, the \a preserve flag can be set to \a true. However, new matrix elements are not
 * initialized!
 */
template< typename Type >  // Data type of the matrix
inline void MatrixMxN<Type>::extend( size_t m, size_t n, bool preserve )
{
   resize( m_+m, n_+n, preserve );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of the matrix.
 *
 * \param elements The new minimum capacity of the sparse matrix.
 * \return void
 *
 * This function increases the capacity of the sparse matrix to at least \a elements elements.
 * The current values of the matrix elements are preserved.
 */
template< typename Type >  // Data type of the matrix
inline void MatrixMxN<Type>::reserve( size_t elements )
{
   if( elements > capacity_ )
   {
      // Allocating a new array
      Type* WALBERLA_RESTRICT tmp = new Type[elements];

      // Replacing the old array
      std::copy( v_, v_+m_*n_, tmp );
      std::swap( tmp, v_ );
      delete [] tmp;

      capacity_ = elements;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Transposing the matrix.
 *
 * \return Reference to the transposed matrix.
 */
template< typename Type >  // Data type of the matrix
inline MatrixMxN<Type>& MatrixMxN<Type>::transpose()
{
   MatrixMxN tmp( trans(*this) );
   swap( tmp );
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the matrix is diagonal.
 *
 * \return \a true if the matrix is diagonal, \a false if not.
 *
 * This function tests whether the matrix is diagonal, i.e. if the non-diagonal elements are
 * default elements. In case of integral or floating point data types, a diagonal matrix has
 * the form

                        \f[\left(\begin{array}{*{5}{c}}
                        aa     & 0      & 0      & \cdots & 0  \\
                        0      & bb     & 0      & \cdots & 0  \\
                        0      & 0      & cc     & \cdots & 0  \\
                        \vdots & \vdots & \vdots & \ddots & 0  \\
                        0      & 0      & 0      & 0      & mn \\
                        \end{array}\right)\f]
 */
template< typename Type >  // Data type of the matrix
inline bool MatrixMxN<Type>::isDiagonal() const
{
   const size_t iend( m_-1 );

   for( size_t i=0; i<iend; ++i ) {
      for( size_t j=i+1; j<n_; ++j ) {
         if( !isDefault( v_[i*6+j] ) || !isDefault( v_[j*6+i] ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the matrix is symmetric.
 *
 * \return \a true if the matrix is symmetric, \a false if not.
 */
template< typename Type >  // Data type of the matrix
inline bool MatrixMxN<Type>::isSymmetric() const
{
   const size_t iend( m_-1 );

   for( size_t i=0; i<iend; ++i ) {
      for( size_t j=i+1; j<n_; ++j ) {
         if( !equal( v_[i*n_+j], v_[j*n_+i] ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the matrix by the scalar value \a scalar (\f$ A=B*s \f$).
 *
 * \param scalar The scalar value for the matrix scaling.
 * \return Reference to the matrix.
 */
template< typename Type >   // Data type of the matrix
template< typename Other >  // Data type of the scalar value
inline MatrixMxN<Type>& MatrixMxN<Type>::scale( Other scalar )
{
   const size_t sqrsize( m_*n_ );
   for( size_t i=0; i<sqrsize; ++i )
      v_[i] *= scalar;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two matrices.
 *
 * \param m The matrix to be swapped.
 * \return void
 * \exception no-throw guarantee.
 */
template< typename Type >  // Data type of the matrix
inline void MatrixMxN<Type>::swap( MatrixMxN& m ) /* throw() */
{
   std::swap( m_, m.m_ );
   std::swap( n_, m.n_ );
   std::swap( capacity_, m.capacity_ );
   std::swap( v_, m.v_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the matrix is aliased with the given address \a alias.
 *
 * \param alias The alias to be checked.
 * \return \a true in case the alias corresponds to this matrix, \a false if not.
 */
template< typename Type >   // Data type of the matrix
template< typename Other >  // Data type of the foreign expression
inline bool MatrixMxN<Type>::isAliased( const Other* alias ) const
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name MatrixMxN operators */
//@{
template< typename Type >
inline bool isnan( const MatrixMxN<Type>& m );

template< typename Type >
inline void reset( MatrixMxN<Type>& m );

template< typename Type >
inline void clear( MatrixMxN<Type>& m );

template< typename Type >
inline bool isDefault( const MatrixMxN<Type>& m );

template< typename Type >
inline const MatrixMxN<Type> inv( const MatrixMxN<Type>& m );

template< typename Type >
inline void swap( MatrixMxN<Type>& a, MatrixMxN<Type>& b ) /* throw() */;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given matrix for not-a-number elements.
 * \ingroup dense_matrix_MxN
 *
 * \param m The matrix to be checked for not-a-number elements.
 * \return \a true if at least one element of the matrix is not-a-number, \a false otherwise.
 */
template< typename Type >  // Data type of the matrix
inline bool isnan( const MatrixMxN<Type>& m )
{
   for( size_t i=0; i<m.rows(); ++i ) {
      for( size_t j=0; j<m.columns(); ++j )
         if( isnan( m(i,j) ) ) return true;
   }
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given dense matrix.
 * \ingroup dense_matrix_MxN
 *
 * \param m The dense matrix to be resetted.
 * \return void
 */
template< typename Type >  // Data type of the matrix
inline void reset( MatrixMxN<Type>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given dense matrix.
 * \ingroup dense_matrix_MxN
 *
 * \param m The dense matrix to be cleared.
 * \return void
 */
template< typename Type >  // Data type of the matrix
inline void clear( MatrixMxN<Type>& m )
{
   m.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given dense matrix is in default state.
 * \ingroup dense_matrix_MxN
 *
 * \param m The dense matrix to be tested for its default state.
 * \return \a true in case the given matrix is component-wise zero, \a false otherwise.
 */
template< typename Type >  // Data type of the matrix
inline bool isDefault( const MatrixMxN<Type>& m )
{
   const size_t sqrsize( m.rows()*m.columns() );
   for( size_t i=0; i<sqrsize; ++i )
      if( !isDefault( m[i] ) ) return false;
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inverting the given dense matrix.
 * \ingroup dense_matrix_MxN
 *
 * \param m The dense matrix to be inverted.
 * \return The inverse of the matrix.
 *
 * This function returns the inverse of the given dense matrix. It has the same effect as
 * calling the getInverse() member function of the matrix.
 */
template< typename Type >  // Data type of the matrix
inline const MatrixMxN<Type> inv( const MatrixMxN<Type>& m )
{
   return m.getInverse();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two matrices.
 * \ingroup dense_matrix_MxN
 *
 * \param a The first matrix to be swapped.
 * \param b The second matrix to be swapped.
 * \return void
 * \exception no-throw guarantee.
 */
template< typename Type >  // Data type of the matrices
inline void swap( MatrixMxN<Type>& a, MatrixMxN<Type>& b ) /* throw() */
{
   a.swap( b );
}
//*************************************************************************************************

} // namespace math
}
