include(test/Benchmark)

Benchmark(AUTHOR CL
	PATH TES/heat_of_reaction/reaction_heat
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 4
	OUTPUT_FILES
		reaction_heat_250.vtu
		reaction_heat_750.vtu
)

Benchmark(AUTHOR CL
	PATH TES/friction/Friction_test
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 8
	OUTPUT_FILES Friction_test_140.vtu
)

Benchmark(AUTHOR CL
	PATH TES/Convection_2D/conv_2D_cart
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 132
	OUTPUT_FILES conv_2D_cart_11830.vtu
)

Benchmark(AUTHOR CL
	PATH TES/Convection_2D_axi/conv_2D_axi
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 122
	OUTPUT_FILES conv_2D_axi_11830.vtu
)
