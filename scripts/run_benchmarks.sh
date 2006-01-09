#!/usr/bin/env bash

SOURCE_LOCATION=`pwd`
SOURCE_LOCATION="$SOURCE_LOCATION/.."
source "$SOURCE_LOCATION/scripts/base/configure_compiler.sh"

# Return code. Will be set to 1 (error) when a build failed
returncode=0

# Parse options
while getopts "a:ed:" opt; do
	case $opt in
		e)
			RUN_EXCEEDING=true
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

# Clean results directory
rm -rf "$SOURCE_LOCATION/../benchmarks/results/*.html"

# Goto sources directory
cd "$SOURCE_LOCATION" >/dev/null

if [ -d ".svn" ] || [ -d "../.svn" ]; then
	# Get svn information
	svn info > "$BUILD_LOCATION/svnInfo.txt"
elif [ -d ".git" ]; then
	# Get git information
	git log HEAD~1..HEAD > "$BUILD_LOCATION/svnInfo.txt"
else
	echo "Aborting: Version information not found."
	exit 1
fi

# Run FEM benchmarks
cd "$BUILD_LOCATION/build_fem"
if  [ $RUN_EXCEEDING ]; then
	ctest -R 'EXCEED' -E 'JOD|Tests|FILE' --output-on-failure -j $NUM_PROCESSORS
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi
	ctest -R 'EXCEEDING_FILECOMPARE' -E 'JOD'
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi

	cd "$BUILD_LOCATION/build_brns"
	ctest -R 'EXCEED' -E 'Tests|FILE' --output-on-failure -j $NUM_PROCESSORS
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi
	ctest -R 'EXCEEDING_FILECOMPARE'
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi

	cd "$BUILD_LOCATION/build_pqc"
	ctest -R 'EXCEED' -E 'Tests|FILE' --output-on-failure -j $NUM_PROCESSORS
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi
	ctest -R 'EXCEEDING_FILECOMPARE'
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi

	cd "$BUILD_LOCATION/build_ipqc"
	ctest -R 'EXCEED' -E 'Tests|FILE' --output-on-failure -j $NUM_PROCESSORS
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi
	ctest -R 'EXCEEDING_FILECOMPARE'
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi

	cd "$BUILD_LOCATION/build_gems"
	ctest -R 'EXCEED' -E 'Tests|FILE' --output-on-failure -j $NUM_PROCESSORS
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi
	ctest -R 'EXCEEDING_FILECOMPARE'
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi

	cd "$BUILD_LOCATION/build_petsc"
	ctest -R 'EXCEED' -E 'Tests|FILE' --output-on-failure -j $NUM_PROCESSORS
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi
	ctest -R 'EXCEEDING_FILECOMPARE'
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi

	cd "$BUILD_LOCATION/build_petsc_gems"
	ctest -R 'EXCEED' -E 'Tests|FILE' --output-on-failure -j $NUM_PROCESSORS
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi
	ctest -R 'EXCEEDING_FILECOMPARE'
	if [ "${?}" -ne "0" ] ; then
		returncode=1
	fi

else
	# Don´t abort on errors
	set +e >/dev/null

	ctest -E 'Tests|FILE|EXCEED' --output-on-failure -j $NUM_PROCESSORS > ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt

	cd "$BUILD_LOCATION/build_brns"
	ctest -E 'FILE|EXCEED|Tests' --output-on-failure -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt

	cd "$BUILD_LOCATION/build_pqc"
	ctest -E 'FILE|EXCEED|Tests' --output-on-failure -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt

	cd "$BUILD_LOCATION/build_ipqc"
	ctest -E 'FILE|EXCEED|Tests' --output-on-failure -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt

	cd "$BUILD_LOCATION/build_gems"
	ctest -E 'FILE|EXCEED|Tests' --output-on-failure -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt

	cd "$BUILD_LOCATION/build_petsc"
	ctest -E 'FILE|EXCEED|Tests' --output-on-failure -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt

	cd "$BUILD_LOCATION/build_petsc_gems"
	ctest -E 'FILE|EXCEED|Tests' --output-on-failure -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt

	cd "$BUILD_LOCATION/build_mpi"
	ctest -E 'FILE|EXCEED|Tests' --output-on-failure -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt

	cd "$BUILD_LOCATION/build_lis"
	ctest -E 'FILE|EXCEED|Tests' --output-on-failure -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt

	cd "$BUILD_LOCATION/build_mkl"
	ctest -E 'FILE|EXCEED|Tests' --output-on-failure -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt

	# Print results
	cat ../benchOut.txt

	cd "$SOURCE_LOCATION/scripts"
	# Send emails on errors
	FILESIZE=$(stat -c %s "$BUILD_LOCATION/svnInfo.txt")
	if [ "$FILESIZE" > "0" ] ; then
	  echo "Running process_benchmark_job.rb"
	  cd process_benchmark_job
	  ruby process_benchmark_job.rb "$BUILD_LOCATION/svnInfo.txt" "$BUILD_LOCATION/benchOut.txt" $HUDSON_EMAIL $1
	  cd ..
	fi

	set -e >/dev/null
fi
cd "$BUILD_LOCATION" >/dev/null

echo "exit code is ${returncode}"
exit ${returncode}
