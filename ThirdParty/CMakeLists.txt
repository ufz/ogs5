add_subdirectory(Brent)
if(OGS_USE_CVODE)
	set(EXAMPLES_ENABLE OFF)
	add_subdirectory(cvode)
	set_target_properties(sundials_cvode_shared sundials_nvecserial_shared PROPERTIES MACOSX_RPATH ON)
endif()

if(OGS_CHEMSOLVER STREQUAL IPQC)
	add_subdirectory(IPQC)
endif()

if(OGS_CHEMSOLVER STREQUAL GEMS)
	add_subdirectory(GEM)
endif()

if(OGS_CHEMSOLVER STREQUAL PQC)
	add_subdirectory ( PQC )
endif()
