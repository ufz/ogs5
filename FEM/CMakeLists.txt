set( HEADERS
	BoundaryCondition.h
	burgers.h
	PhysicalConstant.h
	Constrained.h
	conversion_rate.h
	DistributionInfo.h
	DUMUX.h
	Eclipse.h
	eos.h
	fem_ele.h
	fem_ele_std.h
	fem_ele_vec.h
	FEMCondition.h
	FEMEnums.h
	femlib.h
	files0.h
	GeoInfo.h
	InitialCondition.h
	invariants.h
	LinearFunctionData.h
	mathlib.h
	matrix_class.h
	minkley.h
	Output.h
	pcs_dm.h
	problem.h
	ProcessInfo.h
	prototyp.h
	rf_bc_new.h
	rf_fct.h
	rf_fluid_momentum.h
	rf_ic_new.h
	rf_kinreact.h
	rf_mfp_new.h
	rf_mmp_new.h
	rf_msp_new.h
	rf_node.h
	rf_num_new.h
	rf_out_new.h
	rf_pcs.h
	rf_random_walk.h
	rf_react.h
	rf_react_int.h
	rf_st_new.h
	rf_tim_new.h
	rfmat_cp.h
	solver.h
	SourceTerm.h
	SparseMatrixDOK.h
	Stiff_Bulirsch-Stoer.h
	tools.h
	vtk.h
	OutputTools.h
	ShapeFunctionPool.h
)

set( SOURCES
	BoundaryCondition.cpp
	burgers.cpp
	CAP_IO.cpp
	conversion_rate.cpp
	DistributionInfo.cpp
	DUMUX.cpp
	Eclipse.cpp
	eos.cpp
	fem_ele.cpp
	fem_ele_std.cpp
	fem_ele_std1.cpp
	fem_ele_std_tes.cpp
	fem_ele_std_tneq.cpp
	fem_ele_vec.cpp
	FEMCondition.cpp
	FEMEnums.cpp
	femlib.cpp
	files0.cpp
	GeoInfo.cpp
	InitialCondition.cpp
	invariants.cpp
	LinearFunctionData.cpp
	mathlib.cpp
	matrix_class.cpp
	minkley.cpp
	Output.cpp
	pcs_dm.cpp
	problem.cpp
	ProcessInfo.cpp
	rf_bc_new.cpp
	rf_fct.cpp
	rf_fluid_momentum.cpp
	rf_ic_new.cpp
	rf_kinreact.cpp
	rf_mfp_new.cpp
	rf_mmp_new.cpp
	rf_msp_new.cpp
	rf_node.cpp
	rf_num_new.cpp
	rf_out_new.cpp
	rf_pcs.cpp
	rf_random_walk.cpp
	rf_react.cpp
	rf_react_cap.cpp
	rf_react_int.cpp
	rf_st_new.cpp
	rf_tim_new.cpp
	rfmat_cp.cpp
	SourceTerm.cpp
	SparseMatrixDOK.cpp
	Stiff_Bulirsch-Stoer.cpp
	tools.cpp
	vtk.cpp
	ShapeFunctionPool.cpp
)

if(NOT OGS_LSOLVER STREQUAL PETSC)
	set( SOURCES ${SOURCES} par_ddc.h par_ddc.cpp)
endif()



if(OGS_LSOLVER STREQUAL RF)
	set( HEADERS ${HEADERS} solver.h matrix_routines.h)
	set( SOURCES ${SOURCES} solver.cpp matrix_routines.cpp)
elseif(OGS_LSOLVER STREQUAL PETSC)
	set( SOURCES ${SOURCES} rf_pcs1.cpp fct_mpi.h fct_mpi.cpp)
elseif(OGS_LSOLVER STREQUAL SP)
	set( SOURCES ${SOURCES} equation_class.h equation_class.cpp )
	if (PARALLEL_USE_MPI)
		set(HEADERS ${HEADERS} SplitMPI_Communicator.h )
		set(SOURCES ${SOURCES} SplitMPI_Communicator.cpp )
	endif ()
endif()

if(OGS_CONFIG STREQUAL JFNK)
	set( SOURCES ${SOURCES} rf_pcs1.cpp )
endif()

if(OGS_CHEMSOLVER STREQUAL BRNS)
	set( SOURCES ${SOURCES} rf_REACT_BRNS.h rf_REACT_BRNS.cpp )
endif()

if(OGS_CHEMSOLVER STREQUAL GEMS)
	set( SOURCES ${SOURCES} rf_REACT_GEM.h rf_REACT_GEM.cpp )
endif()

include_directories(
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/FEMIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO
	${CMAKE_SOURCE_DIR}/GCC
	${CMAKE_SOURCE_DIR}/GEM
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MSH
)

# Create the library
add_library( FEM STATIC ${HEADERS} ${SOURCES} )

target_link_libraries(FEM
	Base
	FileIO
	GCC
	GEO
	MSH
)

if(LIS)
	include_directories(${LIS_INCLUDE_DIR})
	target_link_libraries( FEM ${LIS_LIBRARIES} )
	if (WIN32)
		set_target_properties( FEM PROPERTIES LINK_FLAGS "/NODEFAULTLIB:libcmt" )
	endif()
endif()

if(MKL)
	include_directories(${MKL_INCLUDES})
	target_link_libraries( FEM ${MKL_LIBRARIES} )
endif()

if(OGS_CHEMSOLVER STREQUAL PQC)
	include_directories(${CMAKE_SOURCE_DIR}/ThirdParty/PQC)
	target_link_libraries(FEM PQC)
endif()

if(OGS_CHEMSOLVER STREQUAL IPQC)
	include_directories(${CMAKE_SOURCE_DIR}/ThirdParty/IPQC/src)
endif ()

if(PARALLEL_USE_MPI)
	target_link_libraries(FEM ${MPI_CXX_LIBRARIES})
endif()

if (WIN32)
	if (OGS_CHEMSOLVER STREQUAL GEMS)
		target_link_libraries( FEM GEM )
	endif ()
	if (OGS_CHEMSOLVER STREQUAL BRNS)
		target_link_libraries( FEM BrnsDll )
	endif ()
endif ()

if (UNIX)
	if (OGS_CHEMSOLVER STREQUAL BRNS)
		target_link_libraries( FEM ${CMAKE_DL_LIBS} ${LAPACK_LIBRARIES} )
	endif ()
endif ()

if (OGS_CHEMSOLVER STREQUAL CAP)
	include_directories (  ${CMAKE_SOURCE_DIR}/CMP )
	target_link_libraries (
		FEM
		GEO
		MSH
		CMP
	)
	if (UNIX)
		target_link_libraries (
			FEM
			ChemAppCS_64
			EChemAppS_64
		)
	endif ()
endif ()

if(OGS_USE_CVODE)
	target_link_libraries(FEM
		sundials_cvode_static
		sundials_nvecserial_static
	)
endif()
