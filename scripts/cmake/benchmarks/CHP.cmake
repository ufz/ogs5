include(test/Benchmark)

Benchmark(AUTHOR CHP
	PATH FLUID_MOMENTUM/1d_line/1d_line
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES 1d_line_1.vtu
)

Benchmark(AUTHOR CHP
	PATH FLUID_MOMENTUM/1d_quad/1d_quad
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES 1d_quad_1.vtu
)

Benchmark(AUTHOR CHP
	PATH FLUID_MOMENTUM/1d_tri/1d_tri
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES 1d_tri_1.vtu
)

Benchmark(AUTHOR CHP
	PATH FLUID_MOMENTUM/1d_hex/1d_hex
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES 1d_hex_1.vtu
)

Benchmark(AUTHOR CHP
	PATH FLUID_MOMENTUM/1d_pri/1d_pri
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES 1d_pri_1.vtu
)

Benchmark(AUTHOR CHP
	PATH FLUID_MOMENTUM/1d_tet/1d_tet
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES 1d_tet_1.vtu
)

Benchmark(AUTHOR CHP
	PATH FLUID_MOMENTUM/1d_pyra/1d_pyra
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES 1d_pyra_1.vtu
)

Benchmark(AUTHOR CHP
	PATH Anisotropy/permeability/soil_layer
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES soil_layer_30.vtu
)

Benchmark(AUTHOR CHP
	PATH Anisotropy/moleculardiffusion/soil_layer
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES soil_layer_30.vtu
)

Benchmark(AUTHOR CHP
	PATH MULTIPHASE/KueperProblem-PcPnw/kueper
	CONFIG FEM
	RUNTIME 51
	OUTPUT_FILES kueper_5.vtu
)

Benchmark(AUTHOR CHP
	PATH DENSITY-DEPENDENT_FLOW/Elder/elder
	CONFIG FEM
	RUNTIME 252
	OUTPUT_FILES elder_12.vtu
)

Benchmark(AUTHOR CHP
	PATH MULTIPHASE/KueperProblem-PS/kueper
	CONFIG LIS
	RUNTIME 9999
	OUTPUT_FILES kueper_5.vtu
)

Benchmark(AUTHOR CHP
	PATH MULTIPHASE/BuckleyLeverett/h2_line
	CONFIG LIS
	RUNTIME 9999
	OUTPUT_FILES h2_line0080.vtk
)

Benchmark(AUTHOR CHP
	PATH MULTIPHASE/McWhorterProblem/h2_line
	CONFIG MKL
	RUNTIME 1540
	OUTPUT_FILES h2_line_165.vtu
)
