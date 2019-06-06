#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "blockforest/communication/UniformDirectScheme.h"
#include "core/mpi/MPIManager.h"
#include "core/Environment.h"
#include "core/OpenMP.h"
#include "core/mpi/Broadcast.h"
#include "core/math/IntegerFactorization.h"
#include "core/timing/TimingPool.h"
#include "core/waLBerlaBuildInfo.h"
#include "field/communication/StencilRestrictedMPIDatatypeInfo.h"
#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"
#include "field/communication/PackInfo.h"
#include "field/communication/StencilRestrictedPackInfo.h"
#include "field/communication/UniformMPIDatatypeInfo.h"
#include "postprocessing/sqlite/SQLite.h"
#include "python_coupling/CreateConfig.h"
#include "stencil/D3Q7.h"
#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"
#include <functional>

using namespace walberla;
using blockforest::communication::UniformDirectScheme;
using blockforest::communication::UniformBufferedScheme;
using field::communication::UniformMPIDatatypeInfo;
using field::communication::PackInfo;
using field::communication::StencilRestrictedPackInfo;
using field::communication::StencilRestrictedMPIDatatypeInfo;


template<typename Stencil_T>
class SingleMessageBufferedScheme
{
public:
    typedef Stencil_T Stencil;

    SingleMessageBufferedScheme( const weak_ptr_wrapper< StructuredBlockForest > & bf, const int tag = 17953 )
            : blockForest_( bf ), tag_( tag ) {}

    inline void addDataToCommunicate( const shared_ptr< communication::UniformPackInfo > &packInfo )
    {
        tag_ += 1;
        auto newScheme = make_shared< UniformBufferedScheme< Stencil > >( blockForest_, tag_++ );
        newScheme->addDataToCommunicate( packInfo );
        schemes_.push_back( newScheme );
    }

    inline void setLocalMode( const blockforest::LocalCommunicationMode &mode )
    {
        for ( auto &s : schemes_ )
            s->setLocalMOde( mode );
    }

    inline void startCommunication()
    {
        for ( auto &s : schemes_ )
            s->startCommunication();
    }

    inline void wait()
    {
        for ( auto &s : schemes_ )
            s->wait();
    }

private:
    std::vector< shared_ptr< UniformBufferedScheme< Stencil>> > schemes_;
    weak_ptr_wrapper< StructuredBlockForest > blockForest_;
    int tag_;
};


template<typename FieldType, typename Stencil>
void addDataToCommunicate( const shared_ptr< UniformDirectScheme< Stencil > > &scheme, BlockDataID id, uint_t ghostLayers )
{
    scheme->addDataToCommunicate( make_shared< UniformMPIDatatypeInfo< FieldType > >( id, ghostLayers ));
}

template<typename FieldType, typename Scheme_T>
void addDataToCommunicate( const shared_ptr< Scheme_T > &scheme, BlockDataID id, uint_t ghostLayers )
{
    scheme->addDataToCommunicate( make_shared< PackInfo< FieldType > >( id, ghostLayers ));
}


template<typename FieldType, typename Scheme_T>
void addDataToCommunicate( const shared_ptr< Scheme_T > &scheme, BlockDataID id, uint_t ghostLayers, bool )
{
    if ( ghostLayers != 1 )
        scheme->addDataToCommunicate( make_shared< PackInfo< FieldType > >( id, ghostLayers ));
    else
        scheme->addDataToCommunicate( make_shared< StencilRestrictedPackInfo< FieldType, typename Scheme_T::Stencil > >( id ));
}

template<typename FieldType, typename Stencil_T>
void addDataToCommunicate( const shared_ptr< UniformDirectScheme< Stencil_T > > &scheme, BlockDataID id, uint_t ghostLayers, bool )
{
    if ( ghostLayers != 1 )
        scheme->addDataToCommunicate( make_shared< UniformMPIDatatypeInfo< FieldType > >( id, ghostLayers ));
    else
    {
        scheme->addDataToCommunicate( make_shared< StencilRestrictedMPIDatatypeInfo< FieldType, Stencil_T > >( id ));
    }
}

template<typename Scheme1, typename Scheme2>
void addData( const shared_ptr< StructuredBlockForest > &blocks, const config::Config::BlockHandle &configBlock,
              const shared_ptr< Scheme1 > &scheme1, const shared_ptr< Scheme2 > &scheme2,
              uint_t ghostLayers, field::Layout layout )
{
    auto numPdfFields = configBlock.getParameter< uint_t >( "pdf" );
    for ( uint_t i = 0; i < numPdfFields; ++i )
    {
        typedef field::GhostLayerField< real_t, Scheme1::Stencil::Q > Field_T;
        BlockDataID bdId = field::addToStorage< Field_T >( blocks, "pdf", 0.0, layout, ghostLayers );
        addDataToCommunicate< Field_T >( scheme1, bdId, ghostLayers );
        addDataToCommunicate< Field_T >( scheme2, bdId, ghostLayers );
    }


    auto numPdfOptFields = configBlock.getParameter< uint_t >( "pdfOpt" );
    for ( uint_t i = 0; i < numPdfOptFields; ++i )
    {
        typedef field::GhostLayerField< real_t, Scheme1::Stencil::Q > Field_T;
        BlockDataID bdId = field::addToStorage< Field_T >( blocks, "pdfopt", 0.0, layout, ghostLayers );
        addDataToCommunicate< Field_T >( scheme1, bdId, ghostLayers, true );
        addDataToCommunicate< Field_T >( scheme2, bdId, ghostLayers, true );
    }


    auto numVectorFields = configBlock.getParameter< uint_t >( "vector" );
    for ( uint_t i = 0; i < numVectorFields; ++i )
    {
        typedef field::GhostLayerField< real_t, 3 > Field_T;
        BlockDataID bdId = field::addToStorage< Field_T >( blocks, "vector", 0.0, layout, ghostLayers );
        addDataToCommunicate< Field_T >( scheme1, bdId, ghostLayers );
        addDataToCommunicate< Field_T >( scheme2, bdId, ghostLayers );
    }

    auto numScalarFields = configBlock.getParameter< uint_t >( "scalar" );
    for ( uint_t i = 0; i < numScalarFields; ++i )
    {
        typedef field::GhostLayerField< real_t, 1 > Field_T;
        BlockDataID bdId = field::addToStorage< Field_T >( blocks, "scalar", 0.0, layout, ghostLayers );
        addDataToCommunicate< Field_T >( scheme1, bdId, ghostLayers );
        addDataToCommunicate< Field_T >( scheme2, bdId, ghostLayers );
    }
}

template<typename Stencil>
void createCommunication( const shared_ptr< StructuredBlockForest > &blocks,
                          bool buffered, const config::Config::BlockHandle &fieldCfg, uint_t ghostLayers, field::Layout layout,
                          blockforest::LocalCommunicationMode localCommunicationMode, bool singleMessage,
                          std::function< void() > &commStart, std::function< void() > &commWait )
{
    auto directScheme = make_shared< UniformDirectScheme< Stencil > >( blocks, shared_ptr< communication::UniformMPIDatatypeInfo >(), 42 );
    auto bufferedScheme = make_shared< UniformBufferedScheme< Stencil > >( blocks, 4242 );
    auto bufferedSchemeSingle = make_shared< SingleMessageBufferedScheme< Stencil > >( blocks , 24242);

    bufferedScheme->setLocalMode( localCommunicationMode );

    if ( buffered )
    {
        if ( !singleMessage )
        {
            addData( blocks, fieldCfg, directScheme, bufferedScheme, ghostLayers, layout );
            commStart = [=]() { bufferedScheme->startCommunication(); };
            commWait = [=]() { bufferedScheme->wait(); };
        }
        else
        {
            addData( blocks, fieldCfg, directScheme, bufferedSchemeSingle, ghostLayers, layout );
            commStart = [=]() { bufferedSchemeSingle->startCommunication(); };
            commWait = [=]() { bufferedSchemeSingle->wait(); };
        }
    }
    else
    {
        addData( blocks, fieldCfg, directScheme, bufferedScheme, ghostLayers, layout );
        commStart = [=]() { directScheme->startCommunication(); };
        commWait = [=]() { directScheme->wait(); };
    }
}

std::string fromEnv( const char *envVar )
{
    auto env = std::getenv( envVar );
    return env != nullptr ? std::string( env ) : "";
}

int main( int argc, char **argv )
{
    mpi::Environment env( argc, argv );
    int scenarioNr = 0;
    auto mpiManager = mpi::MPIManager::instance();
    for ( auto cfg = python_coupling::configBegin( argc, argv ); cfg != python_coupling::configEnd(); ++cfg )
    {
        if ( mpiManager->isMPIInitialized())
            mpiManager->resetMPI();
        auto config = *cfg;
        auto commCfg = config->getOneBlock( "Communication" );
        auto domainCfg = config->getOneBlock( "Domain" );

        bool cartesianCommunicator = commCfg.getParameter< bool >( "cartesianCommunicator", true );
        if ( !cartesianCommunicator )
            mpiManager->useWorldComm();
        scenarioNr += 1;
        WALBERLA_LOG_INFO_ON_ROOT( "Simulating scenario " << scenarioNr );
        WALBERLA_LOG_INFO_ON_ROOT( *config );

        // ---- Domain Setup ----

        const Vector3< uint_t > cellsPerBlock = domainCfg.getParameter< Vector3< uint_t > >( "cellsPerBlock" );
        const Vector3< real_t > domainWeights = domainCfg.getParameter< Vector3< real_t > >( "domainWeights", Vector3< real_t >( 1.0, 1.0, 1.0 ));
        uint_t blocksPerProcess = domainCfg.getParameter< uint_t >( "blocksPerProcess", 1 );

        auto numProcesses = mpiManager->numProcesses();
        auto processes = math::getFactors3D( uint_c( numProcesses ), domainWeights );
        auto blockDecomposition = math::getFactors3D( uint_c( numProcesses ) * blocksPerProcess, domainWeights );
        auto aabb = AABB( real_t( 0 ), real_t( 0 ), real_t( 0 ),
                          real_c( cellsPerBlock[0] * processes[0] * blocksPerProcess ),
                          real_c( cellsPerBlock[1] * processes[1] * blocksPerProcess ),
                          real_c( cellsPerBlock[2] * processes[2] * blocksPerProcess ));

        auto blocks = blockforest::createUniformBlockGrid( aabb,
                                                           blockDecomposition[0], blockDecomposition[1], blockDecomposition[2],
                                                           cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],
                                                           processes[0], processes[1], processes[2],
                                                           true, true, true, //periodicity
                                                           false // keepGlobalBlockInformation
                                                         );


        // ---- Communication Setup ----
        auto fieldCfg = commCfg.getOneBlock( "Fields" );
        const bool buffered = commCfg.getParameter< bool >( "buffered", true );
        const std::string stencil = commCfg.getParameter< std::string >( "stencil", "D3Q19" );
        const uint_t ghostLayers = commCfg.getParameter< uint_t >( "ghostLayers", 1 );
        const std::string layoutStr = commCfg.getParameter< std::string >( "layout", "fzyx" );
        const std::string localCommModeStr = commCfg.getParameter< std::string >( "localCommunicationMode", "start" );
        const bool singleMessage = commCfg.getParameter< bool >( "singleMessage", false );

        blockforest::LocalCommunicationMode localCommunicationMode;
        if ( localCommModeStr == "start" )
            localCommunicationMode = blockforest::START;
        else if ( localCommModeStr == "wait" )
            localCommunicationMode = blockforest::WAIT;
        else if ( localCommModeStr == "buffer" )
            localCommunicationMode = blockforest::BUFFER;
        else if ( localCommModeStr == "noOptimization" )
            localCommunicationMode = blockforest::NO_OPTIMIZATION;
        else
        {
            WALBERLA_ABORT_NO_DEBUG_INFO( "Unknown localCommunicationMode " << layoutStr << ". Valid values are start, wait, buffer and noOptimization" )
        }


        field::Layout layout;
        if ( layoutStr == "fzyx" )
            layout = field::fzyx;
        else if ( layoutStr == "zyxf" )
            layout = field::zyxf;
        else
        {
            WALBERLA_ABORT_NO_DEBUG_INFO( "Unknown layout string " << layoutStr << ". Valid values are fzyx and zyxf." )
        }

        std::function< void() > commStart;
        std::function< void() > commWait;

        if ( stencil == "D3Q19" )
            createCommunication< stencil::D3Q19 >( blocks, buffered, fieldCfg, ghostLayers, layout, localCommunicationMode, singleMessage, commStart,
                                                   commWait );
        else if ( stencil == "D3Q27" )
            createCommunication< stencil::D3Q27 >( blocks, buffered, fieldCfg, ghostLayers, layout, localCommunicationMode, singleMessage, commStart,
                                                   commWait );
        else if ( stencil == "D3Q7" )
            createCommunication< stencil::D3Q7 >( blocks, buffered, fieldCfg, ghostLayers, layout, localCommunicationMode, singleMessage, commStart, commWait );
        else
        {
            WALBERLA_ABORT_NO_DEBUG_INFO( "Unknown stencil " << stencil << ". Has to be one of D3Q7, D3Q19, D3Q27." )
        }


        // ---- Timing ----
        auto runCfg = config->getOneBlock( "Run" );
        const uint_t warmupIterations = runCfg.getParameter< uint_t >( "warmupIterations", 2 );
              uint_t iterations = runCfg.getParameter< uint_t >( "iterations", 10 );
        const uint_t minIterations = runCfg.getParameter< uint_t >( "minIterations", 2 );
        const uint_t maxIterations = runCfg.getParameter< uint_t >( "maxIterations", 100 );

        const real_t timeForBenchmark = runCfg.getParameter< real_t >( "timeForBenchmark", real_t(-1.0) );
        const uint_t outerIterations = runCfg.getParameter< uint_t >( "outerIterations", 2 );

        const std::string databaseFile = runCfg.getParameter< std::string >( "databaseFile", "FieldCommunication.sqlite" );


        commStart();
        commWait();

        WcTimer warmupTimer;
        warmupTimer.start();
        for ( uint_t warmupCounter = 0; warmupCounter < warmupIterations; ++warmupCounter )
        {
            commStart();
            commWait();
        }
        warmupTimer.end();

        auto estimatedTimePerIteration = warmupTimer.last() / real_c(warmupIterations);
        if( timeForBenchmark > 0 ) {
            iterations = uint_c( timeForBenchmark / estimatedTimePerIteration );
            if( iterations < minIterations )
                iterations = minIterations;
            if( iterations > maxIterations)
                iterations = maxIterations;
        }

        mpi::broadcastObject(iterations);

        WcTimingPool timingPool;
        WALBERLA_MPI_BARRIER();
        WALBERLA_LOG_INFO_ON_ROOT("Running " << outerIterations << " outer iterations of size " << iterations );
        for ( uint_t outerCtr = 0; outerCtr < outerIterations; ++outerCtr )
        {
            timingPool["totalTime"].start();
            for ( uint_t ctr = 0; ctr < iterations; ++ctr )
            {
                timingPool["commStart"].start();
                commStart();
                timingPool["commStart"].end();

                timingPool["commWait"].start();
                commWait();
                timingPool["commWait"].end();
            }
            timingPool["totalTime"].end();
        }

        auto numThreads = omp_get_max_threads();

        auto reducedTimingPool = timingPool.getReduced( timing::REDUCE_TOTAL, 0 );

        WALBERLA_ROOT_SECTION()
        {
            WALBERLA_LOG_RESULT( *reducedTimingPool );

            std::map< std::string, walberla::int64_t > integerProperties;
            std::map< std::string, double > realProperties;
            std::map< std::string, std::string > stringProperties;

            auto databaseBlock = config->getBlock( "Database" );
            if ( databaseBlock )
            {
                for ( auto it = databaseBlock.begin(); it != databaseBlock.end(); ++it )
                    stringProperties[it->first] = it->second;
            }

            realProperties["total_min"] = real_c( timingPool["totalTime"].min()) / real_c( iterations );
            realProperties["total_avg"] = real_c( timingPool["totalTime"].average() / real_c( iterations ));
            realProperties["total_max"] = real_c( timingPool["totalTime"].max() / real_c( iterations ));

            integerProperties["cellsPerBlock0"] = int64_c( cellsPerBlock[0] );
            integerProperties["cellsPerBlock1"] = int64_c( cellsPerBlock[1] );
            integerProperties["cellsPerBlock2"] = int64_c( cellsPerBlock[2] );

            integerProperties["processes0"] = int64_c( processes[0] );
            integerProperties["processes1"] = int64_c( processes[1] );
            integerProperties["processes2"] = int64_c( processes[2] );

            integerProperties["blocks0"] = int64_c( blockDecomposition[0] );
            integerProperties["blocks1"] = int64_c( blockDecomposition[1] );
            integerProperties["blocks2"] = int64_c( blockDecomposition[2] );

            integerProperties["blocksPerProcess"] = int64_c( blocksPerProcess );
            integerProperties["ghostLayers"] = int64_c( ghostLayers );

            integerProperties["fieldsPdf"] = fieldCfg.getParameter< int64_t >( "pdf" );
            integerProperties["fieldsPdfOpt"] = fieldCfg.getParameter< int64_t >( "pdfOpt" );
            integerProperties["fieldsVector"] = fieldCfg.getParameter< int64_t >( "vector" );
            integerProperties["fieldsScalar"] = fieldCfg.getParameter< int64_t >( "scalar" );

            integerProperties["numThreads"] = int64_c( numThreads );
            integerProperties["cartesianCommunicator"] = mpiManager->hasCartesianSetup();

            integerProperties["warmupIterations"] = int64_c( warmupIterations );
            integerProperties["iterations"] = int64_c( iterations );
            integerProperties["outerIterations"] = int64_c( outerIterations );
            integerProperties["buffered"] = int64_c( buffered );
            integerProperties["singleMessage"] = int64_c( singleMessage );

            stringProperties["stencil"] = stencil;
            stringProperties["layout"] = layoutStr;
            stringProperties["localCommunicationMode"] = localCommModeStr;

            stringProperties["SLURM_CLUSTER_NAME"] = fromEnv( "SLURM_CLUSTER_NAME" );
            stringProperties["SLURM_CPUS_ON_NODE"] = fromEnv( "SLURM_CPUS_ON_NODE" );
            stringProperties["SLURM_CPUS_PER_TASK"] = fromEnv( "SLURM_CPUS_PER_TASK" );
            stringProperties["SLURM_JOB_ACCOUNT"] = fromEnv( "SLURM_JOB_ACCOUNT" );
            stringProperties["SLURM_JOB_ID"] = fromEnv( "SLURM_JOB_ID" );
            stringProperties["SLURM_JOB_CPUS_PER_NODE"] = fromEnv( "SLURM_JOB_CPUS_PER_NODE" );
            stringProperties["SLURM_JOB_NAME"] = fromEnv( "SLURM_JOB_NAME" );
            stringProperties["SLURM_JOB_NUM_NODES"] = fromEnv( "SLURM_JOB_NUM_NODES" );
            stringProperties["SLURM_NTASKS"] = fromEnv( "SLURM_NTASKS" );
            stringProperties["SLURM_NTASKS_PER_CORE"] = fromEnv( "SLURM_NTASKS_PER_CORE" );
            stringProperties["SLURM_NTASKS_PER_NODE"] = fromEnv( "SLURM_NTASKS_PER_NODE" );
            stringProperties["SLURM_NTASKS_PER_SOCKET"] = fromEnv( "SLURM_NTASKS_PER_SOCKET" );
            stringProperties["SLURM_TASKS_PER_NODE"] = fromEnv( "SLURM_TASKS_PER_NODE" );

            stringProperties["buildMachine"] = std::string( WALBERLA_BUILD_MACHINE );
            stringProperties["gitVersion"] = std::string( WALBERLA_GIT_SHA1 );
            stringProperties["buildType"] = std::string( WALBERLA_BUILD_TYPE );
            stringProperties["compilerFlags"] = std::string( WALBERLA_COMPILER_FLAGS );

            auto runId = postprocessing::storeRunInSqliteDB( databaseFile, integerProperties, stringProperties, realProperties );
            postprocessing::storeTimingPoolInSqliteDB( databaseFile, runId, timingPool, "TimingRoot" );
            postprocessing::storeTimingPoolInSqliteDB( databaseFile, runId, *reducedTimingPool, "TimingReduced" );
        }

    }

    return 0;
}