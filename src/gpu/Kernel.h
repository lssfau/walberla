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
//! \file Kernel.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Abort.h"
#include "core/debug/Debug.h"
#include "core/FunctionTraits.h"

#include "gpu/GPUWrapper.h"
#include "ErrorChecking.h"

#include <type_traits>
#include <vector>


namespace walberla {
namespace gpu
{


   //*******************************************************************************************************************
   /*! Wrapper class around a GPU kernel, to call kernels also from code not compiled with the device compiler
   *
   * Example:
   * \code
         // Declaration of kernel, implementation has to be in a file compiled with nvcc
         void kernel_func ( double * inputData, int size );

         auto kernel = make_kernel( kernel_func );
         kernel.addParam<double*> ( argument1 );
         kernel.addParam<int>     ( 20 );
         kernel.configure( dim3( 3,3,3), dim3( 4,4,4) );
         kernel();
         // this code is equivalent to:
         kernel_func<<< dim3( 3,3,3), dim3( 4,4,4) >> ( argument1, 20 );
   * \endcode
   *
   * Why use this strange wrapper class instead of the nice kernel call syntax "<<<griddim, blockdim >>>" ??
   *     - This syntax is nice but has to be compiled with the device compiler
   *     - The wrapper allows to compile the kernel call with the host compiler
   *
   * Drawbacks of this class compared to kernel call syntax:
   * Type checking of parameters can only be done at runtime (is done only in Debug mode!).
   * Consider the following example:
   * \code
         // Declaration of kernel, implementation has to be in a file compiled with nvcc
         void kernel_func ( double * inputData, int size );

         auto kernel = make_kernel( kernel_func );
         kernel.addParam<float*>       ( argument1 );
         kernel.addParam<unsigned int> ( 40 );
         kernel.configure( dim3( 3,3,3), dim3( 4,4,4) );
         kernel();
         // this code is equivalent to:
         kernel_func<<< dim3( 3,3,3), dim3( 4,4,4) >> ( argument1, 20 );
   * \endcode
   * The parameter types of the kernel and the parameters added at the gpu::Kernel class do not match.
   * This is only detected when the code is run and was compiled in DEBUG mode!
   *
   *
   * Advantages of this class compared to kernel call syntax: Integrates nicely with waLBerlas field indexing and
   * accessor concepts:
   * \code
         void kernel_func( gpu::SimpleFieldAccessor<double> f );

         auto myKernel = gpu::make_kernel( &kernel_double );
         myKernel.addFieldIndexingParam( gpu::SimpleFieldIndexing<double>::xyz( gpuField ) );
         myKernel();
   * \endcode
   * When using at least one FieldIndexingParameter configure() does not have to be called, since the thread and grid
   * setup is done by the indexing scheme. If two FieldIndexingParameters are passed, the two indexing schemes have to
   * be consistent.
   */
   //*******************************************************************************************************************
   template<typename FuncPtr>
   class Kernel
   {
   public:
      Kernel( FuncPtr funcPtr );

      template<typename T>  void addParam( const T & param );
      template<typename T>  void addFieldIndexingParam( const T & indexing );


      void configure( dim3 gridDim, dim3 blockDim, std::size_t sharedMemSize = 0 );
      void operator() ( gpuStream_t stream = nullptr ) const;


   protected:
      //** Members        **********************************************************************************************
      /*! \name Members  */
      //@{
      FuncPtr funcPtr_;

      bool configured_{ false };
      dim3 gridDim_;
      dim3 blockDim_;
      std::size_t sharedMemSize_{ 0 };

      std::vector< std::vector<char> > params_;
      //@}
      //****************************************************************************************************************


      //** Type checking of parameters **********************************************************************************
      /*! \name Type checking of parameters  */
      //@{
      typedef typename std::remove_pointer<FuncPtr>::type FuncType;

      #define CHECK_PARAMETER_FUNC( Number ) \
      template<typename T> \
      bool checkParameter##Number( typename std::enable_if< (FunctionTraits<FuncType>::arity > Number ), T >::type *  = 0 ) { \
         typedef typename FunctionTraits<FuncType>::template argument<Number>::type ArgType; \
         return std::is_same< T, ArgType >::value; \
      } \
      template<typename T> \
      bool checkParameter##Number( typename std::enable_if< (FunctionTraits<FuncType>::arity <= Number ),T >::type *  = 0 ) { \
         return false; \
      }

      CHECK_PARAMETER_FUNC(0)
      CHECK_PARAMETER_FUNC(1)
      CHECK_PARAMETER_FUNC(2)
      CHECK_PARAMETER_FUNC(3)
      CHECK_PARAMETER_FUNC(4)
      CHECK_PARAMETER_FUNC(5)
      CHECK_PARAMETER_FUNC(6)
      CHECK_PARAMETER_FUNC(7)
      CHECK_PARAMETER_FUNC(8)
      CHECK_PARAMETER_FUNC(9)
      CHECK_PARAMETER_FUNC(10)
      CHECK_PARAMETER_FUNC(11)
      CHECK_PARAMETER_FUNC(12)
      CHECK_PARAMETER_FUNC(13)
      CHECK_PARAMETER_FUNC(14)
      CHECK_PARAMETER_FUNC(15)
      CHECK_PARAMETER_FUNC(16)
      CHECK_PARAMETER_FUNC(17)
      CHECK_PARAMETER_FUNC(18)
      CHECK_PARAMETER_FUNC(19)

      #undef CHECK_PARAMETER_FUNC

      template<typename T> bool checkParameter( uint_t n );
      //@}
      //****************************************************************************************************************
   };


   template<typename FuncPtr>
   Kernel<FuncPtr> make_kernel( FuncPtr funcPtr ) {
      return Kernel<FuncPtr> ( funcPtr );
   }







   //===================================================================================================================
   //
   //  Implementation
   //
   //===================================================================================================================

   template<typename FP>
   Kernel<FP>::Kernel( FP funcPtr )
      : funcPtr_ ( funcPtr ) {}

   template<typename FP>
   template<typename T>
   void Kernel<FP>::addParam( const T & param )
   {
      std::vector<char> paramInfo;
      paramInfo.resize( sizeof(T) );
      std::memcpy ( paramInfo.data(), &param, sizeof(T) );

      WALBERLA_ASSERT( checkParameter<T>( params_.size() ),
                       "gpu::Kernel type mismatch of parameter " << params_.size()  )

      params_.push_back( paramInfo );
   }


   template<typename FP>
   template<typename Indexing>
   void Kernel<FP>::addFieldIndexingParam( const Indexing & indexing )
   {
      configure( indexing.gridDim(), indexing.blockDim() );
      addParam( indexing.gpuAccess() );
   }

   template<typename FP>
   void Kernel<FP>::configure( dim3 gridDim, dim3 blockDim, std::size_t sharedMemSize )
   {
      if ( ! configured_ )
      {
         gridDim_ = gridDim;
         blockDim_ = blockDim;
         sharedMemSize_ = sharedMemSize;
         configured_ = true;
      }
      else
      {
         if ( gridDim.x  != gridDim_.x  || gridDim.y != gridDim_.y   || gridDim.z != gridDim_.z ||
              blockDim.x != blockDim_.x || blockDim.y != blockDim_.y || blockDim.z != blockDim_.z  )
         {
            WALBERLA_ABORT( "Error when configuring gpu::Kernel: Inconsistent setup. " )
         }
      }
   }

   template<typename FP>
   void Kernel<FP>::operator() ( gpuStream_t stream ) const
   {
      // check for correct number of parameter calls
      if ( params_.size() != FunctionTraits<FuncType>::arity ) {
         WALBERLA_ABORT( "Error when calling gpu::Kernel - Wrong number of arguments. " <<
                         "Expected " << FunctionTraits<FuncType>::arity << ", received " << params_.size() )
      }

      // register all parameters
      std::vector<void*> args;
      for( auto paramIt = params_.begin(); paramIt != params_.end(); ++paramIt )  {
         args.push_back( const_cast<char*>(paramIt->data()) );
      }

      // .. and launch the kernel
      static_assert( sizeof(void *) == sizeof(void (*)(void)),
                     "object pointer and function pointer sizes must be equal" );
      WALBERLA_DEVICE_SECTION()
      {
         WALBERLA_GPU_CHECK(gpuLaunchKernel((void*) funcPtr_, gridDim_, blockDim_, args.data(), sharedMemSize_, stream))
      }
   }


   template<typename FP>
   template<typename T>
   bool Kernel<FP>::checkParameter( uint_t n )
   {
      switch (n) {
         case 0: return checkParameter0<T>();
         case 1: return checkParameter1<T>();
         case 2: return checkParameter2<T>();
         case 3: return checkParameter3<T>();
         case 4: return checkParameter4<T>();
         case 5: return checkParameter5<T>();
         case 6: return checkParameter6<T>();
         case 7: return checkParameter7<T>();
         case 8: return checkParameter8<T>();
         case 9: return checkParameter9<T>();
         case 10: return checkParameter10<T>();
         case 11: return checkParameter11<T>();
         case 12: return checkParameter12<T>();
         case 13: return checkParameter13<T>();
         case 14: return checkParameter14<T>();
         case 15: return checkParameter15<T>();
         case 16: return checkParameter16<T>();
         case 17: return checkParameter17<T>();
         case 18: return checkParameter18<T>();
         case 19: return checkParameter19<T>();
         default:
            WALBERLA_ABORT("Too many parameters passed to kernel")
      }
      return false;
   }




} // namespace gpu
} // namespace walberla
