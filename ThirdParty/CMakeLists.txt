add_subdirectory(Brent)
if(OGS_USE_CVODE)
	add_subdirectory(cvode)
	set_target_properties(sundials_cvode_shared sundials_nvecserial_shared PROPERTIES MACOSX_RPATH ON)
endif()

if(OGS_FEM_IPQC)
	add_subdirectory(IPQC)
endif()

if(OGS_FEM_GEMS OR OGS_FEM_PETSC_GEMS)
	add_subdirectory(GEM)
endif()

if(OGS_FEM_PQC)
	add_subdirectory ( PQC )
endif()
