include(test/Benchmark)

Benchmark(AUTHOR TK
	PATH GROUNDWATER_FLOW/uc_quad
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES uc_quad_domain_quad.tec
)

Benchmark(AUTHOR TK
	PATH GROUNDWATER_FLOW/uc_pris
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES uc_pris_domain_pris.tec
)

Benchmark(AUTHOR TK
	PATH GROUNDWATER_FLOW/uc_tri
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES uc_tri_domain_GROUNDWATER_FLOW_tri.tec
)

Benchmark(AUTHOR TK
	PATH GROUNDWATER_FLOW/q_quad
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES q_quad_time_POINT4.tec
)

Benchmark(AUTHOR TK
	PATH GROUNDWATER_FLOW/q_hex
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES q_hex_time_POINT4.tec
)

Benchmark(AUTHOR TK
	PATH H/sat_1D/H_sat_flow_1d
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES H_sat_flow_1d_domain_LIQUID_FLOW_line.tec
)

Benchmark(AUTHOR TK
	PATH H/sat_2D/H_sat_flow_K_ortho
	CONFIG FEM
	RUNTIME 26
	OUTPUT_FILES H_sat_flow_K_ortho_domain_tri.tec
)

Benchmark(AUTHOR TK
	PATH H/HetGWFlow/2D/2D1P-GWFlow
	CONFIG FEM
	RUNTIME 62
	OUTPUT_FILES 2D1P-GWFlow_domain_quad.tec
)

Benchmark(AUTHOR TK
	PATH H/Theis_1D/h_quad_axisym
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES h_quad_axisym_time_POINT5_LIQUID_FLOW.tec
)

Benchmark(AUTHOR TK
	PATH H/Theis_2D/Thies_quad_2d
	CONFIG FEM
	RUNTIME 14
	OUTPUT_FILES Thies_quad_2d_time_POINT5_GROUNDWATER_FLOW.tec
)

Benchmark(AUTHOR TK
	PATH H/HetGWFlow/3D/3D_HGW
	CONFIG FEM
	RUNTIME 203
	OUTPUT_FILES 3D_HGW_domain_hex.tec
)

Benchmark(AUTHOR TK
	PATH H/Theis/GWF_Theis_1-5D/GWF_Theis
	CONFIG FEM
	RUNTIME 3
	OUTPUT_FILES
		GWF_Theis_domain_GROUNDWATER_FLOW_line.tec
		GWF_Theis_time_OBS_GROUNDWATER_FLOW.tec
)

Benchmark(AUTHOR TK
	PATH H/Theis/GWF_Theis_2-5D/GWF_Theis
	CONFIG FEM
	RUNTIME 6
	OUTPUT_FILES
		GWF_Theis_domain_GROUNDWATER_FLOW_quad.tec
		GWF_Theis_time_OBS_GROUNDWATER_FLOW.tec
)

Benchmark(AUTHOR TK
	PATH H/Theis/GWF_Theis_2D/GWF_Theis_2d
	CONFIG FEM
	RUNTIME 15
	OUTPUT_FILES
		GWF_Theis_2d_domain_GROUNDWATER_FLOW_quad.tec
		GWF_Theis_2d_domain_GROUNDWATER_FLOW_tri.tec
		GWF_Theis_2d_time_OBS_GROUNDWATER_FLOW.tec
)

Benchmark(AUTHOR TK
	PATH H/Theis/GWF_Theis_3D/GWF_Theis_3D
	CONFIG FEM
	RUNTIME 33
	OUTPUT_FILES
		GWF_Theis_3D_domain_GROUNDWATER_FLOW_hex.tec
		GWF_Theis_3D_domain_GROUNDWATER_FLOW_pris.tec
		GWF_Theis_3D_time_OBS_GROUNDWATER_FLOW.tec
)

Benchmark(AUTHOR TK
	PATH H/Theis/LF_Theis_1-5D/LF_Theis
	CONFIG FEM
	RUNTIME 3
	OUTPUT_FILES
		LF_Theis_domain_LIQUID_FLOW_line.tec
		LF_Theis_time_OBS_LIQUID_FLOW.tec
)

Benchmark(AUTHOR TK
	PATH H/Theis/LF_Theis_2-5D/LF_Theis
	CONFIG FEM
	RUNTIME 7
	OUTPUT_FILES
		LF_Theis_domain_LIQUID_FLOW_quad.tec
		LF_Theis_time_OBS_LIQUID_FLOW.tec
)

Benchmark(AUTHOR TK
	PATH H/Theis/LF_Theis_2D/LF_Theis_2D
	CONFIG FEM
	RUNTIME 16
	OUTPUT_FILES
		LF_Theis_2D_domain_LIQUID_FLOW_quad.tec
		LF_Theis_2D_domain_LIQUID_FLOW_tri.tec
		LF_Theis_2D_time_OBS_LIQUID_FLOW.tec
)

Benchmark(AUTHOR TK
	PATH H/Theis/LF_Theis_3D/LF_Theis_3D
	CONFIG FEM
	RUNTIME 30
	OUTPUT_FILES
		LF_Theis_3D_domain_LIQUID_FLOW_hex.tec
		LF_Theis_3D_domain_LIQUID_FLOW_pris.tec
		LF_Theis_3D_time_OBS1_LIQUID_FLOW.tec
		LF_Theis_3D_time_OBS2_LIQUID_FLOW.tec
)
