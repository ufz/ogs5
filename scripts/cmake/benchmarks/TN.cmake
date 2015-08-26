include(test/Benchmark)

Benchmark(AUTHOR TN
	PATH T2HC/heat_transfer/heat_transfer
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 1
	OUTPUT_FILES
		heat_transfer_0.vtu
		heat_transfer_50.vtu
		heat_transfer_100.vtu
		heat_transfer_500.vtu
)

Benchmark(AUTHOR TN
	PATH T2HC/heat_of_reaction/reaction_heat
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 2
	OUTPUT_FILES
		reaction_heat_250.vtu
		reaction_heat_750.vtu
)

Benchmark(AUTHOR TN
	PATH T2HC/friction/Friction_test
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 14
	OUTPUT_FILES Friction_test_140.vtu
)

Benchmark(AUTHOR TN
	PATH T2HC/Convection_2D/conv_2D_cart
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 232
	OUTPUT_FILES conv_2D_cart_11830.vtu
)

Benchmark(AUTHOR TN
	PATH T2HC/Convection_2D_axi/conv_2D_axi
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 183
	OUTPUT_FILES conv_2D_axi_11830.vtu
)
