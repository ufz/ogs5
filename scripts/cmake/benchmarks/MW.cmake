include(test/Benchmark)

Benchmark(AUTHOR MW
	PATH H_us/Darcy/unconf_WO_rch/unconf
	CONFIG FEM
	RUNTIME 16
	OUTPUT_FILES
		unconf_2.vtu
		unconf_ply_bottom_t1.tec
)

Benchmark(AUTHOR MW
	PATH H_us/Darcy/unconf_W_rch/unconf
	CONFIG FEM
	RUNTIME 11
	OUTPUT_FILES
		unconf_2.vtu
		unconf_ply_bottom_t1.tec
)

Benchmark(AUTHOR MW
	PATH DENSITY-DEPENDENT_FLOW/goswami-clement/constrBC_PressAsHead_tri/HM
	CONFIG FEM
	RUNTIME 361
	OUTPUT_FILES
		HM0047.vtk
		HM0122.vtk
		HM0161.vtk
)

Benchmark(AUTHOR MW
	PATH DENSITY-DEPENDENT_FLOW/goswami-clement/wholeBC_PressAsHead_tri/HM
	CONFIG FEM
	RUNTIME 373
	OUTPUT_FILES
		HM0042.vtk
		HM0119.vtk
		HM0159.vtk
)
