/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_CL_TESTS_H
#define OGS_CL_TESTS_H

#include "BuildInfo.h"

#include "gtest.h"

#include <cstdlib>
#include <sstream>
#include <iostream>

#include <unistd.h>

namespace {

std::string tmpStr = ( BuildInfo::OGS_EXECUTABLE ); // passed by CMakeLists.txt


// The fixture for testing class Ogs.
class OgsTest : public ::testing::Test {
 protected:
  OgsTest() {
    // Set-up work.
    // see if system() works
    if ( ! system(NULL) )
      {
	std::cout << "System call failed." << std::endl;
	exit(1);
      }
  }
  void callOgs( std::string, char*, int sze );//char);
};

  void OgsTest::callOgs( std::string callStr, char* outChar, int sze)
  {
    int outfd[2];
    int infd[2];

    pipe(outfd); // ogs writes to
    pipe(infd); //  gtest reads from

    pid_t pid = fork();

    if( ! pid )    // the child process
      {

	close(STDOUT_FILENO);     // close the file descriptors inherited
	close(STDIN_FILENO);      //   from the parent

	dup2(outfd[0], STDIN_FILENO);  // reassign the new fd
	dup2(infd[1], STDOUT_FILENO);

	close(outfd[0]);          // close all uneeded fd
	close(outfd[1]);
	close(infd[0]);
	close(infd[1]);

	tmpStr.append( callStr );
	system( tmpStr.c_str() );  // call ogs here
      }
    else  // in the parent process
      {

	close(outfd[0]); // close ends of the pipes used by the child
	close(infd[1]);
	close(outfd[1]); // not sending to child's STDIN so close it too

	outChar[ read(infd[0], outChar, sze ) ] = 0; // Read from child’s stdout

	//printf("%s",outChar);   //   for debugging
	close(infd[0]);
      }
  }


  TEST_F(OgsTest, ShowsDisplay)
  {
    char input[187];
    int sze = sizeof(input);
    std::string strCall = " a";
    OgsTest::callOgs( strCall, input, sze );
    EXPECT_EQ("\n\
          ###################################################\n\
          ##                                               ##\n\
          ##               OpenGeoSys-Project              ##\n",
	      std::string( input ) );
  }

  TEST_F(OgsTest, ShowsBuildInfo)
  {
    char input[15];
    int sze = sizeof(input);
    std::string strCall = " -b";
    OgsTest::callOgs( strCall, input, sze );
    EXPECT_EQ("ogs version: 5.",
	      std::string( input ) );
  }

  TEST_F(OgsTest, ShowsHelp)
  {
    char input[101];
    int sze = sizeof(input);
    std::string strCall = " --help";
    OgsTest::callOgs( strCall, input, sze );
    std::stringstream ss (std::stringstream::in | std::stringstream::out);

    ss << "Usage: ogs [MODEL_ROOT] [OPTIONS]\n"
       << "Where OPTIONS are:\n"
       << "  -h [--help]       print this message and exit" << std::endl;

    EXPECT_EQ( ss.str(),
	       std::string( input ) );
  }

  TEST_F(OgsTest, ShowsVersion)
  {
    char input[2];
    int sze = sizeof(input);
    std::string strCall = " --version";
    OgsTest::callOgs( strCall, input, sze );
    EXPECT_EQ("5.",
	      std::string( input ) );
  }

}  // namespace
/*
int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
*/
#endif //OGS_CL_TESTS_H
