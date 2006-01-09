#!/usr/bin/env bash

CMAKE_ARGS=""
BUILD_ARGS=""

SOURCE_LOCATION=`pwd`
SOURCE_LOCATION="$SOURCE_LOCATION/../.."

# Parse options
while getopts "a:d:o:" opt; do
	case $opt in
		a)
			source $SOURCE_LOCATION/scripts/base/architecture_option_win.sh
			;;
		d)
			BUILD_LOCATION="$SOURCE_LOCATION/$OPTARG"
			;;
		o)
			CMAKE_ARGS="$OPTARG"
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

if [ "$OSTYPE" == 'msys' ]; then
	CMAKE_ARGS="$CMAKE_ARGS -DOGS_PACKAGING_ZIP=ON"
fi
if [ "$OSTYPE" != 'msys' ]; then
	BUILD_ARGS="-- -j $NUM_PROCESSORS"
fi

cmake $CMAKE_ARGS -DOGS_PACKAGING=ON -DOGS_NO_EXTERNAL_LIBS=ON \
  -DCMAKE_BUILD_TYPE=Release -G "$CMAKE_GENERATOR" $SOURCE_LOCATION
cmake $SOURCE_LOCATION

cmake --build . --config Release --target package $BUILD_ARGS

cd "$SOURCE_LOCATION/scripts/build"

if [ "$OSTYPE" == 'msys' ]; then
	exit 0
fi
