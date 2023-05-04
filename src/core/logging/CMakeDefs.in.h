//======================================================================================================================
/*!
 *  \file   CMakeDefs.in.h
 *  \brief  Definitions for core/logging configured by cmake
 */
//======================================================================================================================

#pragma once

#cmakedefine WALBERLA_LOGLEVEL_WARNING
#cmakedefine WALBERLA_LOGLEVEL_INFO
#cmakedefine WALBERLA_LOGLEVEL_PROGRESS
#cmakedefine WALBERLA_LOGLEVEL_DETAIL
#cmakedefine WALBERLA_LOGLEVEL_TRACING

#define WALBERLA_LOGLEVEL ${WALBERLA_LOGLEVEL}
#define WALBERLA_LOGLEVEL_STRING "${WALBERLA_LOGLEVEL}"