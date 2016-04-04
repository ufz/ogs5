/**
 * \file BuildInfo.h
 * 24/08/2010 LB Initial implementation
 * #defines which gets set through CMake
 */

#ifndef BUILDINFO_H
#define BUILDINFO_H

#include <string>

namespace BuildInfo
{
extern const std::string GIT_COMMIT_INFO;
extern const std::string GIT_BRANCH_INFO;
extern const std::string BUILD_TIMESTAMP;
extern const std::string CMAKE_SYSTEM;
extern const std::string CMAKE_SYSTEM_PROCESSOR;
extern const std::string CMAKE_CXX_COMPILER;
extern const std::string GCC_VERSION;
extern const std::string CMAKE_GENERATOR;
extern const std::string CMAKE_BUILD_TYPE;
extern const std::string CMAKE_CXX_FLAGS;
extern const std::string CMAKE_CXX_FLAGS_RELEASE;
extern const std::string CMAKE_CXX_FLAGS_DEBUG;
extern const std::string SOLVER_PACKAGE_NAME;

extern const std::string SOURCEPATH;
extern const std::string BUILDPATH;
extern const std::string TESTDATAPATH;
extern const std::string OGS_VERSION;
extern const std::string OGS_DATE;
extern const std::string OGS_EXECUTABLE;
extern const std::string PUT_TMP_DIR_IN;
}

#endif // BUILDINFO_H
