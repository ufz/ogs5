include(test/Benchmark)

Benchmark(AUTHOR PV_WW
	PATH PETSc/hm1_1Dbeam/hm1_1Dbeam
	CONFIG PETSC
	NUM_PROCESSORS 2
	OUTPUT_FILES
		hm1_1Dbeam_domain_DEFORMATION_hex_0.tec
		hm1_1Dbeam_domain_DEFORMATION_hex_1.tec
		hm1_1Dbeam_domain_LIQUID_FLOW_hex_0.tec
		hm1_1Dbeam_domain_LIQUID_FLOW_hex_1.tec
)

Benchmark(AUTHOR PV_WW
	PATH PETSc/tm1_3Dorigin/tm1_3Dorigin
	CONFIG PETSC
	NUM_PROCESSORS 2
	OUTPUT_FILES
		tm1_3Dorigin_domain_hex_0.tec
		tm1_3Dorigin_domain_hex_1.tec
)

Benchmark(AUTHOR PV_WW
	PATH PETSc/m1_3Dload/m1_3Dload
	CONFIG PETSC
	NUM_PROCESSORS 4
	OUTPUT_FILES
		m1_3Dload_domain_DEFORMATION_hex_0.tec
		m1_3Dload_domain_DEFORMATION_hex_1.tec
		m1_3Dload_domain_DEFORMATION_hex_2.tec
		m1_3Dload_domain_DEFORMATION_hex_3.tec
)
