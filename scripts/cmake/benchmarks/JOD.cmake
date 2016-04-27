include(test/Benchmark)

Benchmark(AUTHOR JOD
	PATH OVERLAND_FLOW/gian_quad
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES
		gian_quad_domain_OVERLAND_FLOW_quad.tec
		gian_quad_time_POINT4_OVERLAND_FLOW.tec
		gian_quad_time_POINT5_OVERLAND_FLOW.tec
)

Benchmark(AUTHOR JOD
	PATH OVERLAND_FLOW/gian_tri
	CONFIG FEM
	RUNTIME 2
	OUTPUT_FILES
			gian_tri_domain_OVERLAND_FLOW_tri.tec
			gian_tri_time_POINT4_OVERLAND_FLOW.tec
			gian_tri_time_POINT5_OVERLAND_FLOW.tec
)

Benchmark(AUTHOR JOD
	PATH OVERLAND_FLOW/govin_line
	CONFIG FEM
	RUNTIME 52
	OUTPUT_FILES
		govin_line_domain_OVERLAND_FLOW_line.tec
		govin_line_time_POINT0_OVERLAND_FLOW.tec
)

Benchmark(AUTHOR JOD
	PATH OVERLAND_FLOW/govin_quad
	CONFIG FEM
	RUNTIME 171
	OUTPUT_FILES
		govin_quad_domain_OVERLAND_FLOW_quad.tec
		govin_quad_time_POINT0_OVERLAND_FLOW.tec
		govin_quad_time_POINT1_OVERLAND_FLOW.tec
)

Benchmark(AUTHOR JOD
	PATH OVERLAND_FLOW/Wool_quad
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES Wool_quad_time_POINT1_OVERLAND_FLOW.tec
)
