#pragma once

#include "core/debug/Debug.h"
#include "blockforest/Block.h"
#include "blockforest/StructuredBlockForest.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "blockforest/communication/UniformDirectScheme.h"
#include "field/communication/StencilRestrictedMPIDatatypeInfo.h"
#include "field/communication/UniformMPIDatatypeInfo.h"
#include "cuda/communication/GPUPackInfo.h"
#include "cuda/communication/UniformGPUScheme.h"
#include "cuda/communication/MemcpyPackInfo.h"
#include "UniformGridGPU_PackInfo.h"


using namespace walberla;


enum CommunicationSchemeType {
    GPUPackInfo_Baseline = 0,
    GPUPackInfo_Streams = 1,
    UniformGPUScheme_Baseline = 2,
    UniformGPUScheme_Memcpy = 3,
    MPIDatatypes = 4,
    MPIDatatypesFull = 5
};


template<typename StencilType, typename GPUFieldType >
class UniformGridGPU_Communication
{
public:
    explicit UniformGridGPU_Communication(weak_ptr_wrapper<StructuredBlockForest> bf, const BlockDataID & bdId,
                                          CommunicationSchemeType commSchemeType, bool cudaEnabledMPI = false)
        : _commSchemeType(commSchemeType), _cpuCommunicationScheme(nullptr), _gpuPackInfo(nullptr),
          _gpuCommunicationScheme(nullptr), _directScheme(nullptr)
    {
        auto generatedPackInfo = make_shared<pystencils::UniformGridGPU_PackInfo>( bdId );
        auto memcpyPackInfo = make_shared< cuda::communication::MemcpyPackInfo< GPUFieldType > >( bdId );
        auto dataTypeInfo = make_shared< field::communication::StencilRestrictedMPIDatatypeInfo< GPUFieldType, StencilType > >( bdId );
        auto dataTypeInfoFull = make_shared< field::communication::UniformMPIDatatypeInfo<GPUFieldType> >( bdId );

        switch(_commSchemeType)
        {
            case GPUPackInfo_Baseline:
                _gpuPackInfo = make_shared< cuda::communication::GPUPackInfo< GPUFieldType > >( bdId );
                _cpuCommunicationScheme = make_shared< blockforest::communication::UniformBufferedScheme< StencilType > >( bf );
                _cpuCommunicationScheme->addPackInfo( _gpuPackInfo );
                break;
            case GPUPackInfo_Streams:
                _gpuPackInfo = make_shared< cuda::communication::GPUPackInfo< GPUFieldType > >( bdId );
                _cpuCommunicationScheme = make_shared< blockforest::communication::UniformBufferedScheme< StencilType > >( bf );
                _cpuCommunicationScheme->addPackInfo( _gpuPackInfo );
                break;
            case UniformGPUScheme_Baseline:
                _gpuCommunicationScheme = make_shared< cuda::communication::UniformGPUScheme< StencilType > >( bf, cudaEnabledMPI );
                _gpuCommunicationScheme->addPackInfo( generatedPackInfo );
                break;
            case UniformGPUScheme_Memcpy:
                _gpuCommunicationScheme = make_shared< cuda::communication::UniformGPUScheme< StencilType > >( bf, cudaEnabledMPI );
                _gpuCommunicationScheme->addPackInfo( memcpyPackInfo );
                break;
            case MPIDatatypes:
                if( ! cudaEnabledMPI ) {
                    WALBERLA_ABORT("MPI datatype-based communication not possible if no cudaEnabledMPI is available.");
                }
                _directScheme = make_shared< blockforest::communication::UniformDirectScheme< StencilType > >( bf, dataTypeInfo );
                break;
            case MPIDatatypesFull:
                if( ! cudaEnabledMPI ) {
                    WALBERLA_ABORT("MPI datatype-based communication not possible if no cudaEnabledMPI is available.");
                }
                _directScheme = make_shared< blockforest::communication::UniformDirectScheme< StencilType > >( bf, dataTypeInfoFull );
                break;
            default:
                WALBERLA_ABORT("Invalid GPU communication scheme specified!");
        }
    }

    UniformGridGPU_Communication(UniformGridGPU_Communication &) = delete;

    void operator()( cudaStream_t communicationStream = 0 )
    {
        startCommunication( communicationStream );
        wait( communicationStream );
    }

    void startCommunication( cudaStream_t communicationStream )
    {
        switch( _commSchemeType )
        {
            case GPUPackInfo_Streams:
                // Set communication stream to enable asynchronous operations
                // in GPUPackInfo.
                WALBERLA_ASSERT_NOT_NULLPTR( _gpuPackInfo );
                _gpuPackInfo->setCommunicationStream( communicationStream );
                // Start communication using UniformBufferedScheme
                WALBERLA_ASSERT_NOT_NULLPTR( _cpuCommunicationScheme );
                _cpuCommunicationScheme->startCommunication();
                break;
            case GPUPackInfo_Baseline:
                // Start communication using UniformBufferedScheme
                WALBERLA_ASSERT_NOT_NULLPTR( _cpuCommunicationScheme );
                _cpuCommunicationScheme->startCommunication();
                break;
            case UniformGPUScheme_Baseline:
                WALBERLA_ASSERT_NOT_NULLPTR( _gpuCommunicationScheme );
                _gpuCommunicationScheme->startCommunication( communicationStream );
                break;
            case UniformGPUScheme_Memcpy:
                WALBERLA_ASSERT_NOT_NULLPTR( _gpuCommunicationScheme );
                _gpuCommunicationScheme->startCommunication( communicationStream );
                break;
            case MPIDatatypes:
            case MPIDatatypesFull:
                WALBERLA_ASSERT_NOT_NULLPTR( _directScheme );
                _directScheme->startCommunication();
                break;
        }
    }

    void wait( cudaStream_t communicationStream )
    {
        switch( _commSchemeType )
        {
            case GPUPackInfo_Baseline:
                WALBERLA_ASSERT_NOT_NULLPTR( _cpuCommunicationScheme );
                _cpuCommunicationScheme->wait();
                break;
            case GPUPackInfo_Streams:
                WALBERLA_ASSERT_NOT_NULLPTR( _cpuCommunicationScheme );
                _gpuPackInfo->setCommunicationStream( communicationStream );
                _cpuCommunicationScheme->wait();
                break;
            case UniformGPUScheme_Baseline:
                WALBERLA_ASSERT_NOT_NULLPTR( _gpuCommunicationScheme );
                _gpuCommunicationScheme->wait( communicationStream );
                break;
            case UniformGPUScheme_Memcpy:
                WALBERLA_ASSERT_NOT_NULLPTR( _gpuCommunicationScheme );
                _gpuCommunicationScheme->wait( communicationStream );
                break;
            case MPIDatatypes:
            case MPIDatatypesFull:
                WALBERLA_ASSERT_NOT_NULLPTR( _directScheme );
                _directScheme->wait();
                break;
        }
    }

private:
    CommunicationSchemeType _commSchemeType;
    shared_ptr< blockforest::communication::UniformBufferedScheme< StencilType > > _cpuCommunicationScheme;
    shared_ptr< cuda::communication::GPUPackInfo< GPUFieldType > > _gpuPackInfo;
    shared_ptr< cuda::communication::UniformGPUScheme< StencilType > > _gpuCommunicationScheme;
    shared_ptr< blockforest::communication::UniformDirectScheme<StencilType> > _directScheme;
};
