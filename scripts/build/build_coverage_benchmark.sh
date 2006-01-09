#!/usr/bin/env bash

SOURCE_LOCATION=`pwd`
SOURCE_LOCATION="$SOURCE_LOCATION/../.."

# Parse options
while getopts "a:d:" opt; do
	case $opt in
		a)
			source $SOURCE_LOCATION/scripts/base/architecture_option_win.sh
			;;
		d)
			BUILD_LOCATION="$SOURCE_LOCATION/$OPTARG"
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

CMAKE_ARGS=""
BUILD_ARGS=""

if [ "$OSTYPE" != 'msys' ]; then
	BUILD_ARGS="-- -j $NUM_PROCESSORS"
fi

cmake $CMAKE_ARGS -DOGS_DONT_USE_QT=ON -DOGS_BUILD_TESTS=ON -DOGS_COVERAGE=ON \
	"$CMAKE_GENERATOR" $SOURCE_LOCATION
cmake $SOURCE_LOCATION

cmake --build . --config Debug $BUILD_ARGS
cmake --build . --config Debug --target benchmark_coverage $BUILD_ARGS

cd "$SOURCE_LOCATION/scripts/build"

if [ "$OSTYPE" == 'msys' ]; then
	exit 0
fi