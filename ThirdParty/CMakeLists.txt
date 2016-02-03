add_subdirectory(Brent)
if(OGS_USE_CVODE)
	add_subdirectory(cvode)
	set_target_properties(sundials_cvode_shared sundials_nvecserial_shared PROPERTIES MACOSX_RPATH ON)
endif()

if(OGS_CONFIG STREQUAL IPQC)
	add_subdirectory(IPQC)
endif()

if(OGS_CONFIG STREQUAL GEMS OR OGS_CONFIG STREQUAL PETSC_GEMS)
	add_subdirectory(GEM)
endif()

if(OGS_CONFIG STREQUAL PQC)
	add_subdirectory ( PQC )
endif()
