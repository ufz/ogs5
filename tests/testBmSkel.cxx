/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_MINI_BM_TESTS_H
#define OGS_MINI_BM_TESTS_H

using namespace std;

/* Adjust verbosity here */
//#define NDEBUG
//#define NINFO
#define NWARNING
#define NERROR
#include "BuildInfo.h"
#include "logging.h"

#include "gtest.h"

#include <cstdlib>
#include <fstream>   // for ifstream
#include <sstream>   // are these needed?
#include <iostream>  // are these needed?
#include <stdio.h>   // for remove()

#include <unistd.h> // also used for rmdir

#include <sys/types.h> // for opendir()
#include <dirent.h>    // for opendir()
#include <errno.h>

namespace {

std::string tmpStr = ( BuildInfo::OGS_EXECUTABLE );

class MinBMTest : public ::testing::Test {
  /**
     Test fixture for minimal models.

     Copies a model into a temporary directory, runs the model, tests various
     aspects of the output, and deletes the temporary directory.

     The rationale is to exercise the code on a smaller scale than the
     'benchmark tests'.
  */
  std::string TmpDirectory; // path to a temporary directory
public:
  char *tmpDirectory;

protected:

  MinBMTest()
  {
    /** Set-up work. */
    // see if system() works
    if ( ! system(NULL) )
      {
	std::cout << "System call failed." << std::endl;
	exit(1);
      }

    TmpDirectory.append( BuildInfo::PUT_TMP_DIR_IN ); // passed by CMakeLists.txt
    PRINT_DEBUG( BuildInfo::PUT_TMP_DIR_IN )
    // TODO: put this somewhere configurable
    TmpDirectory.append( "data/tmpTestOutput_XXXXXX" );
    // convert to char * for mkdtemp()
    tmpDirectory  = ( char * )( TmpDirectory.c_str() );

    // make a temporary directory to put working files in
    char * newDir;
    newDir = mkdtemp( tmpDirectory );

    if( newDir == NULL) {
       std::string errStr = ( "Error creating directory : " +
			      std::string( strerror(errno) ) );
       PRINT_ERROR( errStr << std::endl );
      exit(1); }

    PRINT_DEBUG( tmpDirectory );
  }

  void copyModelToTmpDir( std::string modelDir )
  {
    // open source directory
    DIR * bdir = opendir( modelDir.c_str() );
    std::string fullPath;

    if( ! bdir )
      {
	PRINT_WARNING( "Directory " << modelDir << " not opened for copying\n"
		       << std::endl; );
	  // << std::endl;
      }
    else
      {
	dirent * aDirEnt;
	aDirEnt = readdir( bdir );

	while( aDirEnt )
	  {
	    std::string fname = aDirEnt->d_name;
	    if( fname != "." && fname != ".." )
	    {
	      std::string toFpath = tmpDirectory;
	      toFpath +=  "/" + fname;
	      std::string fromFpath = modelDir + "/" + fname;
	      PRINT_INFO( "Copying file " << fromFpath  << " to " << toFpath
			  << std::endl );
	      // copy each file to the temporary directory
	      std::ifstream f1( fromFpath.c_str(), std::fstream::binary );
	      std::ofstream f2( toFpath.c_str(),
				std::fstream::trunc | std::fstream::binary );
	      f2 << f1.rdbuf();
	    }
	    // get the next file (or NULL)
	    aDirEnt = readdir( bdir );
	  }
	closedir( bdir );
      }
  }

  void callOgs( std::string callStr )
  {
    /** Call ogs with callStr as the argument and wait until it has finished.
     */
    /* This block should work but for these models the bug arising from the
       absolute path to the model in:

       `/path/to/ogs /path/to/model > /dev/null`

       causes a segfault later on. As a workaround, the call is changed to:

       `cd /path/to; /path/to/ogs model > /dev/null`

       but this is a bug that should be fixed.  Thought to work already for
       MULTIPHASE_FLOW models.

    string baseStr = tmpStr;
    baseStr.append( callStr );
    PRINT_DEBUG( baseStr );
    system( baseStr.c_str() ); // call ogs here
    */
    callStr = ""; // this is to avoid compiler warning 'unused parameter..'
    std::string TmpDirectory = tmpDirectory;
    string baseStr = "cd " + TmpDirectory + "; " + tmpStr + " a > /dev/null";
    PRINT_DEBUG( baseStr );
    system( baseStr.c_str() );  // call ogs here
  }

  void callOgs( std::string callStr, char* outChar, int sze)
  {
    /** Call ogs with callStr as the argument on a parent process and read sze
	characters of the echo output into outChar via a child process.

	(Sometimes, if the echo output ends in endl rather then \n, not all
	output is read.)
    */

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
	PRINT_DEBUG( tmpStr.c_str() );
	system( tmpStr.c_str() );  // call ogs here
      }
    else  // in the parent process
      {

	close(outfd[0]); // close ends of the pipes used by the child
	close(infd[1]);
	close(outfd[1]); // not sending to child's STDIN so close it too

	outChar[ read(infd[0], outChar, sze ) ] = 0; // Read child’s stdout

	//printf("%s",outChar);   //   for debugging
	close(infd[0]);
      }
  }

  virtual void SetUp() {
    // Code here will be called immediately after the constructor (right
    // before each test).
    // copy the model files to the temporaray directory

  }

  virtual void TearDown() {
    // Code here will be called immediately after each test (right
    // before the destructor).

    // Remove the temporary directory and contents

    DIR * adir = opendir( tmpDirectory );
    std::string fullPath;

    if( ! adir )
      {
	PRINT_WARNING( "Directory " << tmpDirectory
		       << " not opened for deleting"
		       << std::endl );
      }
    else
      {
	dirent * aDirEnt;
	aDirEnt = readdir( adir );
	PRINT_INFO( "Reading directory " << tmpDirectory << std::endl );
	while( aDirEnt )
	  {
	    std::string fname = aDirEnt->d_name;
	    if( fname != "." && fname != ".." )
	      {
		std::string toFpath = tmpDirectory;
		toFpath += "/" + fname;
		PRINT_INFO( "Deleting file from " << toFpath << std::endl );
		if( remove( toFpath.c_str() ) != 0 ) {
		  // untested
		  std::string errStr = ( "Error deleting file : " +
					 std::string( strerror(errno) ) );
		  PRINT_ERROR( errStr << std::endl; );
		}
		else
		  { PRINT_INFO( "File successfully deleted" << std::endl; );}
	      }
		// get the next file (or NULL)
		aDirEnt = readdir( adir );
	  }
	closedir( adir );
	remove( tmpDirectory );
      }
  }
};

  TEST_F(MinBMTest, MostBasicModelEverMade)
  {
    // setup the model for running
	char result[256];
	strcpy(result,(BuildInfo::SOURCEPATH).c_str());
	strcat(result,"/tests/data/bmskel");
    copyModelToTmpDir( result );

    // run ogs, ignoring echo output
    // need to make "ogs /path/to/tmpdir/a"
    std::string TmpDirectory = tmpDirectory;
    // ( the dev/null is just to turn off the echo) > /dev/null
    std::string callStr = " " + TmpDirectory + "/a > /dev/null";
    callOgs( callStr );

    // get the results
    std::string toFpath = tmpDirectory;
    toFpath += "/a_time_POINT0_HEAT_TRANSPORT.tec";
    PRINT_DEBUG( std::cout << "toFpath ==== " << toFpath );
    std::ifstream ifs( toFpath.c_str() );

    // read the whole file into a string
    std::string gfs( ( std::istreambuf_iterator< char > ( ifs ) ),
		     std::istreambuf_iterator< char > () );

    // this is what should have resulted
    std::string ans = " TITLE = \"Time curves in points\"\n";
    ans += " VARIABLES = \"TIME \"  \"TEMPERATURE1\" \n";
    ans += " ZONE T=\"POINT=POINT0\"\n";
    ans += "0.000000000000e+00 2.730000000000e+02 \n";
    ans += "1.000000000000e+00 2.740000000000e+02 \n";

    // test
    EXPECT_EQ( ans, gfs);
  }

}  // namespace
/*
int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
*/
#endif //OGS_MINI_BM_TESTS_H
