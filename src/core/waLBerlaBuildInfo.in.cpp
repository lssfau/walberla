//======================================================================================================================
/*!
 *  \file   waLBerlaBuildInfo.in.cpp
 *  \brief  Configured by CMake, contains information about the current build
 */
//======================================================================================================================


namespace walberla {
namespace core {
namespace buildinfo {
   
const char * gitSHA1()       { return "@WALBERLA_GIT_SHA1@";       }
const char * buildType()     { return "@WALBERLA_BUILD_TYPE@";     }
const char * compilerFlags() { return "@WALBERLA_COMPILER_FLAGS@"; }
const char * buildMachine()  { return "@WALBERLA_BUILD_MACHINE@";  }
const char * sourceDir()     { return "@walberla_SOURCE_DIR@";     }
const char * binaryDir()     { return "@walberla_BINARY_DIR@";     }

} // namespace buildinfo
} // namespace core
} // namespace walberla 

