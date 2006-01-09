#!/usr/bin/env bash

SOURCE_LOCATION=`pwd`
SOURCE_LOCATION="$SOURCE_LOCATION/../.."

# Parse options
while getopts "a:d:m" opt; do
	case $opt in
		a)
			source "$SOURCE_LOCATION/scripts/base/architecture_option_win.sh"
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

source "$SOURCE_LOCATION/scripts/base/check_architecture_option_win.sh"
source "$SOURCE_LOCATION/scripts/base/check_and_cleanup_build_directory.sh"
source "$SOURCE_LOCATION/scripts/base/configure_compiler.sh"

# Executables will be copied to Release/
mkdir -p Release

BUILD_ARGS=""
if [ "$OSTYPE" != 'msys' ]; then
	BUILD_ARGS="-- -j $NUM_PROCESSORS"
fi

# Return code. Will be set to 1 (error) when a build failed
returncode=0

exe_dir="bin"
# Set configurations
if [ "$OSTYPE" == 'msys' ]; then
	configs=("FEM" "SP" "PQC")                                 # Windows
	exe_dir="bin/Release"
elif [[ "$OSTYPE" == darwin* ]]; then
	configs=("FEM" "SP")                                       # Mac
elif [[ "$mpi" == true ]]; then
	configs=("MPI" "PETSC" "PETSC_GEMS")                       # MPI Linux
	export CC=mpicc
	export CXX=mpic++
else
	configs=("FEM" "SP" "GEMS" "PQC" "BRNS" "MKL" "LIS")       # Linux
fi

# Iterate over configurations
for config in ${configs[@]}
do
	cmake_args=""
	if [ "$config" = "FEM" ] ; then
		cmake_args="-DOGS_BUILD_TESTS=ON"
		config_cmake="OGS_FEM"
		exe_name="ogs"
		build_dir="build_fem"
	else
		config_low=$( echo "$config"|tr -s '[:upper:]' '[:lower:]' )
		config_cmake="OGS_FEM_$config"
		exe_name="ogs_$config_low"
		build_dir="build_$config_low"
	fi

	if [ "$HOSTNAME" = "dev2.intern.ufz.de" ] ; then
		cmake_args="$cmake_args -DBoost_USE_STATIC_LIBS=OFF"
	fi

	# Cleanup
	rm -rf $build_dir

	# Create build directory
	mkdir -p $build_dir && cd $build_dir

	# Run CMake
	cmake -D$config_cmake=ON -DOGS_DONT_USE_QT=ON -DCMAKE_BUILD_TYPE=Release $cmake_args -G "$CMAKE_GENERATOR" "$SOURCE_LOCATION"
	cmake -G "$CMAKE_GENERATOR" "$SOURCE_LOCATION"

	# Build
	cmake --build . --config Release $BUILD_ARGS

	# Remember exit code
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi

	# Copy executable
	cp $exe_dir/ogs ../Release/$exe_name

	cd ..

done

echo "exit code is ${returncode}"
exit ${returncode}
