#!/usr/bin/env bash

SOURCE_LOCATION=`pwd`
SOURCE_LOCATION="$SOURCE_LOCATION/../.."

# Parse options
while getopts "a:d:m" opt; do
	case $opt in
		a)
			source $SOURCE_LOCATION/scripts/base/architecture_option_win.sh
			;;
		d)
			BUILD_LOCATION="$SOURCE_LOCATION/$OPTARG"
			;;
		m)
			mpi=true
			;;
		\?)
			echo "Invalid option: -$OPTARG"
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument."
			exit 1
			;;
	esac
done

source $SOURCE_LOCATION/scripts/base/check_architecture_option_win.sh
source $SOURCE_LOCATION/scripts/base/check_and_cleanup_build_directory.sh
source $SOURCE_LOCATION/scripts/base/configure_compiler.sh

# Executables will be copied to Release/
mkdir -p Release

BUILD_ARGS=""
if [ "$OSTYPE" != 'msys' ]; then
	BUILD_ARGS="-- -j $NUM_PROCESSORS"
fi

# Return code. Will be set to 1 (error) when a build failed
returncode=0

build_dir=build_fileconverter

# Cleanup
rm -rf $build_dir

# Create build directory
mkdir -p $build_dir && cd $build_dir

# Run CMake
cmake -D$config_cmake=ON -DOGS_BUILD_UTILITIES=ON -DCMAKE_BUILD_TYPE=Release $cmake_args -G "$CMAKE_GENERATOR" $SOURCE_LOCATION
cmake -G "$CMAKE_GENERATOR" $SOURCE_LOCATION

# Build
cmake --build . --config Release --target OGSFileConverter $BUILD_ARGS

# Remember exit code
if [ "${?}" -ne "0" ] ; then
	returncode=1
fi

cd ..

echo "exit code is ${returncode}"
exit ${returncode}
