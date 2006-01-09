cmake_minimum_required(VERSION 2.8.3)

project(IPhreeqc)

set(SOURCES
	src/IPhreeqc.cpp
	src/IPhreeqcLib.cpp
	src/CSelectedOutput.cpp
	src/fwrap.cpp
	src/Var.c

	src/phreeqcpp/ChartHandler.cpp
	src/phreeqcpp/ChartObject.cpp
	src/phreeqcpp/CurveObject.cpp
	src/phreeqcpp/ExchComp.cxx
	src/phreeqcpp/Exchange.cxx
	src/phreeqcpp/GasComp.cxx
	src/phreeqcpp/GasPhase.cxx
	src/phreeqcpp/ISolution.cxx
	src/phreeqcpp/ISolutionComp.cxx
	src/phreeqcpp/Keywords.cpp
	src/phreeqcpp/KineticsComp.cxx
	src/phreeqcpp/NameDouble.cxx
	src/phreeqcpp/NumKeyword.cxx
	src/phreeqcpp/PBasic.cpp
	src/phreeqcpp/PHRQ_base.cxx
	src/phreeqcpp/PHRQ_io.cpp
	src/phreeqcpp/PHRQ_io_output.cpp
	src/phreeqcpp/PPassemblage.cxx
	src/phreeqcpp/PPassemblageComp.cxx
	src/phreeqcpp/Parser.cxx
	src/phreeqcpp/Phreeqc.cpp
	src/phreeqcpp/Pressure.cxx
	src/phreeqcpp/Reaction.cxx
	src/phreeqcpp/ReadClass.cxx
	src/phreeqcpp/SS.cxx
	src/phreeqcpp/SSassemblage.cxx
	src/phreeqcpp/SScomp.cxx
	src/phreeqcpp/Solution.cxx
	src/phreeqcpp/SolutionIsotope.cxx
	src/phreeqcpp/StorageBin.cxx
	src/phreeqcpp/StorageBinList.cpp
	src/phreeqcpp/Surface.cxx
	src/phreeqcpp/SurfaceCharge.cxx
	src/phreeqcpp/SurfaceComp.cxx
	src/phreeqcpp/System.cxx
	src/phreeqcpp/Temperature.cxx
	src/phreeqcpp/Use.cpp
	src/phreeqcpp/Utils.cxx
	src/phreeqcpp/advection.cpp
	src/phreeqcpp/basicsubs.cpp
	src/phreeqcpp/cl1.cpp
	src/phreeqcpp/cvdense.cpp
	src/phreeqcpp/cvode.cpp
	src/phreeqcpp/cxxKinetics.cxx
	src/phreeqcpp/cxxMix.cxx
	src/phreeqcpp/dense.cpp
	src/phreeqcpp/dumper.cpp
	src/phreeqcpp/dw.cpp
	src/phreeqcpp/gases.cpp
	src/phreeqcpp/input.cpp
	src/phreeqcpp/integrate.cpp
	src/phreeqcpp/inverse.cpp
	src/phreeqcpp/isotopes.cpp
	src/phreeqcpp/kinetics.cpp
	src/phreeqcpp/mainsubs.cpp
	src/phreeqcpp/model.cpp
	src/phreeqcpp/nvector.cpp
	src/phreeqcpp/nvector_serial.cpp
	src/phreeqcpp/parse.cpp
	src/phreeqcpp/phqalloc.cpp
	src/phreeqcpp/pitzer.cpp
	src/phreeqcpp/pitzer_structures.cpp
	src/phreeqcpp/prep.cpp
	src/phreeqcpp/print.cpp
	src/phreeqcpp/read.cpp
	src/phreeqcpp/readtr.cpp
	src/phreeqcpp/runner.cpp
	src/phreeqcpp/sit.cpp
	src/phreeqcpp/smalldense.cpp
	src/phreeqcpp/spread.cpp
	src/phreeqcpp/step.cpp
	src/phreeqcpp/structures.cpp
	src/phreeqcpp/sundialsmath.cpp
	src/phreeqcpp/tally.cpp
	src/phreeqcpp/tidy.cpp
	src/phreeqcpp/transport.cpp
	src/phreeqcpp/utilities.cpp
	src/phreeqcpp/SelectedOutput.cpp
	src/phreeqcpp/UserPunch.cpp
)

if(WIN32)
	set(SOURCES ${SOURCES}
		src/fwrap2.cpp
		src/fwrap3.cpp
		src/fwrap4.cpp
		src/fwrap5.cpp
		src/fwrap6.cpp
		src/fwrap7.cpp
		src/fwrap8.cpp
	)

	add_definitions(
		-DUSE_PHRQ_ALLOC
		-D_CRT_SECURE_NO_DEPRECATE
		-D_SCL_SECURE_NO_WARNINGS
		-D_MBCS
	)
endif()

add_definitions(
	-DSWIG_SHARED_OBJ
	-DPHREEQC_CPP
	-DPHREEQC_CLASS
	-DPIC
)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${CMAKE_CURRENT_SOURCE_DIR}/src/phreeqcpp
)

add_library( iphreeqc STATIC ${SOURCES} )

include (CheckIncludeFiles)
CHECK_INCLUDE_FILES(sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES(sys/stat.h HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES(stdlib.h HAVE_STDLIB_H)
CHECK_INCLUDE_FILES(string.h HAVE_STRING_H)
CHECK_INCLUDE_FILES(memory.h HAVE_MEMORY_H)
CHECK_INCLUDE_FILES(strings.h HAVE_STRINGS_H)
CHECK_INCLUDE_FILES(inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES(stdint.h HAVE_STDINT_H)
CHECK_INCLUDE_FILES(unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILES(dlfcn.h HAVE_DLFNC_H)
CHECK_INCLUDE_FILES(float.h HAVE_FLOAT_H)
CHECK_INCLUDE_FILES(limits.h HAVE_LIMITS_H)
CHECK_INCLUDE_FILES(stddef.h HAVE_STDDEF_H)

include (CheckFunctionExists)
CHECK_FUNCTION_EXISTS(memmove HAVE_MEMMOVE)
CHECK_FUNCTION_EXISTS(memset HAVE_MEMSET)
CHECK_FUNCTION_EXISTS(strchr HAVE_STRCHR)
CHECK_FUNCTION_EXISTS(strcspn HAVE_STRCSPN)
CHECK_FUNCTION_EXISTS(strtol HAVE_STRTOL)
CHECK_FUNCTION_EXISTS(floor HAVE_FLOOR)
CHECK_FUNCTION_EXISTS(pow HAVE_POW)

add_definitions(
	-DHAVE_SYS_TYPES_H=${HAVE_SYS_TYPES_H}
	-DHAVE_SYS_STAT_H=${HAVE_SYS_STAT_H}
	-DHAVE_STDLIB_H=${HAVE_STDLIB_H}
	-DHAVE_STRING_H=${HAVE_STRING_H}
	-DHAVE_MEMORY_H=${HAVE_MEMORY_H}
	-DHAVE_STRINGS_H=${HAVE_STRINGS_H}
	-DHAVE_INTTYPES_H=${HAVE_INTTYPES_H}
	-DHAVE_STDINT_H=${HAVE_STDINT_H}
	-DHAVE_UNISTD_H=${HAVE_UNISTD_H}
	-DHAVE_DLFNC_H=${HAVE_DLFNC_H}
	-DHAVE_FLOAT_H=${HAVE_FLOAT_H}
	-DHAVE_LIMITS_H=${HAVE_LIMITS_H}
	-DHAVE_STDDEF_H=${HAVE_STDDEF_H}
	-DHAVE_MEMMOVE=${HAVE_MEMMOVE}
	-DHAVE_MEMSET=${HAVE_MEMSET}
	-DHAVE_STRCHR=${HAVE_STRCHR}
	-DHAVE_STRCSPN=${HAVE_STRCSPN}
	-DHAVE_STRTOL=${HAVE_STRTOL}
	-DHAVE_FLOOR=${HAVE_FLOOR}
	-DHAVE_POW=${HAVE_POW}
)