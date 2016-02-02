include(test/Benchmark)

Benchmark(AUTHOR FC
	PATH C/monod/rt1
	CONFIG BRNS
	RUNTIME 175
	OUTPUT_FILES rt1_ply_PLY_X1_t0.tec
)

Benchmark(AUTHOR FC
	PATH C/1d_degradation_network/ce
	CONFIG BRNS
	RUNTIME 9999
	OUTPUT_FILES ce_domain_line.tec
)
