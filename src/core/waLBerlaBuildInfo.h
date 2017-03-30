//======================================================================================================================
/*!
 *  \file   waLBerlaBuildInfo.h
 *  \brief  Functions that provide information about the current build
 */
//======================================================================================================================


#pragma once 

namespace walberla {
namespace core {
namespace buildinfo {
   
const char * gitSHA1();
const char * buildType();
const char * compilerFlags();
const char * buildMachine();
const char * sourceDir();
const char * binaryDir();

} // namespace buildinfo
} // namespace core
} // namespace walberla 

#define WALBERLA_GIT_SHA1       walberla::core::buildinfo::gitSHA1()
#define WALBERLA_BUILD_TYPE     walberla::core::buildinfo::buildType()
#define WALBERLA_COMPILER_FLAGS walberla::core::buildinfo::compilerFlags()
#define WALBERLA_BUILD_MACHINE  walberla::core::buildinfo::buildMachine()
#define WALBERLA_SOURCE_DIR     walberla::core::buildinfo::sourceDir()
#define WALBERLA_BUILD_DIR      walberla::core::buildinfo::binaryDir()

