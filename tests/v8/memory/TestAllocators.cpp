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
//! \file TestAllocators.cpp
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
//======================================================================================================================

#include "core/all.h"

#include "walberla/v8/Memory.hpp"
#include "walberla/v8/Testing.hpp"


#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)
#include "TestAllocators.h"
#endif

namespace test_allocators
{
using namespace walberla;
using namespace walberla::v8;

template<typename T>
void checkAlignment(T* p, size_t alignment){
    std::uintptr_t addr = reinterpret_cast<std::uintptr_t>(p);
    testing::assert_equal((((addr & -addr) % alignment)), static_cast<std::uintptr_t>(0));
}

// ***************************************************************************************************
// Tests
// ***************************************************************************************************
template<memory::MemTag mem_tag, typename T>
struct testAlignment{    
    void operator()(){
        using MT = memory::MemoryTraits< mem_tag, T >;
        using AllocatorType = MT::AllocatorType;
        const size_t alignment = 2048;
        const size_t n = 10;
        AllocatorType allocator(0, alignment);
        T* p = allocator.allocate(n);
        checkAlignment<T>(p, alignment);
        allocator.deallocate(p, n);
    }
};

template<memory::MemTag mem_tag, typename T>
struct testAlignmentMT{    
    void operator()(){    
        using MT = memory::MemoryTraits< mem_tag, T >;
        using AllocatorType = MT::AllocatorType;
        const size_t alignment = MT::alignment();
        const size_t n = 10;
        AllocatorType allocator{ MT::getAllocator(0) };
        T* p = allocator.allocate(n);
        checkAlignment<T>(p, alignment);
        allocator.deallocate(p, n);
    }
};

template<memory::MemTag mem_tag, typename T>
struct testOffsetAlignment{   
    void operator()(){  
        using MT = memory::MemoryTraits< mem_tag, T >;
        using AllocatorType = MT::AllocatorType;
        const size_t alignment = 2048;
        const size_t n = 100;
        const size_t offset = 7 * sizeof(T);
        AllocatorType allocator(offset, alignment);
        T* p = allocator.allocate(n);
        checkAlignment<T>(p + 7, alignment);
        allocator.deallocate(p, n);
    }
};

template<memory::MemTag mem_tag, typename T>
struct testInterchange{   
    void operator()(){  
        using MT = memory::MemoryTraits< mem_tag, T >;
        using AllocatorType = MT::AllocatorType;
        const size_t n = 100;
        AllocatorType a1(4 * sizeof(T), 256);
        AllocatorType a2(9 * sizeof(T), 64);

        testing::assert_equal(a1==a2, true);
        testing::assert_equal(a1!=a2, false);

        T* p = a1.allocate(n);     // Allocate with a1
        a2.deallocate(p, n);       // Deallocate with a2
    }
};

template<memory::MemTag mem_tag, typename T>
struct testBadAlloc{ 
    void operator()(){  
        using MT = memory::MemoryTraits< mem_tag, T >;
        using A = MT::AllocatorType;

        A a(100 * sizeof(T), alignof(T));

        testing::throws< std::bad_alloc >([&](){
            //  overflow
            [[maybe_unused]] T* p = a.allocate(size_t(-1));
        });

        testing::throws< std::bad_alloc >([&](){
            // offset larger than data
            [[maybe_unused]] T* p = a.allocate(40);
        });
    }
};


template<memory::MemTag mem_tag, typename T>
struct testDefaultAndCopy{
    void operator()(){  
        using MT = memory::MemoryTraits< mem_tag, T >;
        using AllocatorType = MT::AllocatorType;

        static_assert(std::is_default_constructible_v<AllocatorType>);
        static_assert(std::is_copy_constructible_v<AllocatorType>);

        // Default constructor
        AllocatorType a1;
        T* p = a1.allocate(10);
        checkAlignment<T>(p, alignof(T));
        a1.deallocate(p, 10);

        // Copy constrcutor
        AllocatorType a2(a1);
        T* p2 = a2.allocate(10);
        checkAlignment<T>(p2, alignof(T));
        a1.deallocate(p2, 10);
    }
};

template<memory::MemTag mem_tag, typename T1, typename T2>
struct testRebindAndConvert{
    void operator()(){  
        static_assert(!std::is_same_v< T1, T2 >);

        using MT1 = memory::MemoryTraits< mem_tag, T1 >;
        using A = MT1::AllocatorType;

        using MT2 = memory::MemoryTraits< mem_tag, T2 >;
        using B_expected = MT2::AllocatorType;

        // Rebind A to a different value_type T2
        using B = A::template rebind<T2>::other;

        static_assert(std::is_same_v< B, B_expected >);
        static_assert(std::is_constructible_v<A, B>);
        static_assert(std::is_constructible_v<B, A>);

        B b(8, 64);
        A a(b);

        testing::assert_equal(B(a)==b, true);
        testing::assert_equal(A(b)==a, true);
    }
};

template<memory::MemTag mem_tag, typename T>
struct testSTLVector{
    static_assert(mem_tag::isHostAccessible, "testSTLVector requires memory to be host accessible");
    
    // Checks if the values on device match the values set on host
    void testUnified([[maybe_unused]] T* ptr, [[maybe_unused]] size_t size){
        #if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)

        bool* flag = nullptr;
        WALBERLA_GPU_CHECK(gpuMallocManaged(&flag, sizeof(bool)));

        checkVector(ptr, size, flag);
        testing::assert_true(*flag);

        #endif
    }

    void operator()(){
        using MT = memory::MemoryTraits< mem_tag, T >;
        using AllocatorType = MT::AllocatorType;

        AllocatorType allocator1(6 * sizeof(T), 128);
        std::vector<T, AllocatorType> vec1(60, allocator1);
        vec1.resize(20); 
        vec1.shrink_to_fit();
        vec1.clear();
        vec1.reserve(200);

        for(size_t i = 0; i < 200; ++i) vec1.push_back(static_cast<T>(i));

        checkAlignment<T>(vec1.data() + 6, 128);

        AllocatorType allocator2(10 * sizeof(T), 256);
        std::vector<T, AllocatorType> vec2(20, allocator2);
        
        vec2 = vec1;

        if constexpr(std::is_same_v<mem_tag,memory::memtag::unified>){
            testUnified(vec1.data(), vec1.size());
            testUnified(vec2.data(), vec2.size());
        }

        testing::assert_equal(vec1==vec2, true);
    }
};

template<memory::MemTag mem_tag, typename T>
struct testSTLMap{
    void operator()(){  
        using PairType = std::pair<const T, T>;
        using MT = memory::MemoryTraits< mem_tag, PairType >;
        using AllocatorType = MT::AllocatorType;

        AllocatorType allocator1(0, 256);
        std::map<T, T, std::less<T>, AllocatorType> m1(std::less<T>(), allocator1);

        for (int i = 0; i < 50; ++i)
            m1.emplace(i, static_cast<T>(i));

        for (int i = 0; i < 50; ++i)
            testing::assert_equal(static_cast<int>(m1.at(i)), i);
        
        AllocatorType allocator2(0, 128);
        std::map<T, T, std::less<T>, AllocatorType> m2(std::less<T>(), allocator2);
        
        m2 = m1;
        testing::assert_equal(m1==m2, true);
    }
};

// ***************************************************************************************************
// Test all
// ***************************************************************************************************

template <template<memory::MemTag, typename> class Test>
void runForAllMemTags()
{
    Test<memory::memtag::host, double>()();
    Test<memory::memtag::unified, double>()();
    Test<memory::memtag::pinned, double>()();

#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)
    Test<memory::memtag::device, double>()();
#endif
}

void testRebindAndConvertAll(){
    testRebindAndConvert<memory::memtag::host, int, double>()();
    testRebindAndConvert<memory::memtag::unified, int, double>()();
    testRebindAndConvert<memory::memtag::pinned, int, double>()();

#if defined(WALBERLA_BUILD_WITH_CUDA) || defined(WALBERLA_BUILD_WITH_HIP)
    testRebindAndConvert<memory::memtag::device, int, double>()();
#endif
}

void testSTLVectorAll(){
    testSTLVector<memory::memtag::host, double>()();
    testSTLVector<memory::memtag::unified, double>()();
    testSTLVector<memory::memtag::pinned, double>()();
}

void testSTLMapAll(){
    testSTLMap<memory::memtag::host, double>()();
    testSTLMap<memory::memtag::unified, double>()();
    testSTLMap<memory::memtag::pinned, double>()();
}

} // namespace test_allocators


// ***************************************************************************************************
// ***************************************************************************************************

int main(int argc, char** argv){
   walberla::mpi::Environment env{ argc, argv };

   using namespace test_allocators;

   return walberla::v8::testing::TestsRunner(
             {
                { "testAlignment",           &runForAllMemTags< testAlignment > },
                { "testAlignmentMT",         &runForAllMemTags< testAlignmentMT > },
                { "testOffsetAlignment",     &runForAllMemTags< testOffsetAlignment > },
                { "testInterchange",         &runForAllMemTags< testInterchange > },
                { "testBadAlloc",            &runForAllMemTags< testBadAlloc > },
                { "testDefaultAndCopy",      &runForAllMemTags< testDefaultAndCopy > },
                { "testRebindAndConvert",    &testRebindAndConvertAll },
                { "testSTLVector",           &testSTLVectorAll },
                { "testSTLMap",              &testSTLMapAll },
             })
      .run(argc, argv);
}
