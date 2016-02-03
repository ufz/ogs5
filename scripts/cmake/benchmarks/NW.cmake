include(test/Benchmark)

Benchmark(AUTHOR NW
	PATH T/t3d/t3d
	CONFIG FEM
	OUTPUT_FILES t3d_time_POINT12.tec
	RUNTIME 2
)

Benchmark(AUTHOR NW
	PATH T/1d_thermal_expansion/exp1
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES
		exp1_time_POINT_LEFT.tec
		exp1_time_POINT_RIGHT.tec
)

Benchmark(AUTHOR NW
	PATH TM/tm_01_3Du
	CONFIG FEM
	RUNTIME 4
	OUTPUT_FILES tm_01_3Du_domain_hex.tec
)

Benchmark(AUTHOR NW
	PATH TM/tm2d
	CONFIG FEM
	RUNTIME 2
	OUTPUT_FILES
		tm2d_domain_quad.tec
		tm2d_time_POINT2.tec
		tm2d_time_POINT3.tec
		tm2d_time_POINT6.tec
		tm2d_time_POINT7.tec
		tm2d_time_POINT8.tec
)

Benchmark(AUTHOR NW
	PATH TM/TM_axi
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES
		TM_axi_domain_tri.tec
		TM_axi_ply_PLY_B_t1.tec
)

Benchmark(AUTHOR NW
	PATH NUMERICS/DISCRETE_FEATURES/InclinedFeature/H_incline_45r_line
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES
		H_incline_45r_line_domain_ele.tec
		H_incline_45r_line_domain_line.tec
		H_incline_45r_line_ply_PLY_0_t0.tec
)

Benchmark(AUTHOR NW
	PATH NUMERICS/DISCRETE_FEATURES/InclinedFeature/H_incline_45r_quad
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES
		H_incline_45r_quad_domain_ele.tec
		H_incline_45r_quad_domain_quad.tec
		H_incline_45r_quad_ply_PLY_1_t0.tec
		H_incline_45r_quad.pvd
		H_incline_45r_quad_0.vtu
		H_incline_45r_quad_1.vtu
	)

Benchmark(AUTHOR NW
	PATH NUMERICS/DISCRETE_FEATURES/Lauwerier/Lauwerier
	CONFIG FEM
	RUNTIME 6
	OUTPUT_FILES
		Lauwerier_ply_V2_t0.tec
		Lauwerier_ply_FRACTURE_t1.tec
)

Benchmark(AUTHOR NW
	PATH NUMERICS/SUPG/T_adv_diff_steady_SUPG_line
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES T_adv_diff_steady_SUPG_line_ply_PLY_0_t0.tec
)

Benchmark(AUTHOR NW
	PATH NUMERICS/FEM_FCT/mass_adv_line
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES mass_adv_line_ply_PLY_0_t0.tec
)

Benchmark(AUTHOR NW
	PATH TM/tm_02_3Du
	CONFIG FEM
	RUNTIME 10
	OUTPUT_FILES tm_02_3Du_domain_hex.tec
)

Benchmark(AUTHOR NW
	PATH TM/tm3d
	CONFIG FEM
	RUNTIME 18
	OUTPUT_FILES
		tm3d_domain_tet.tec
		tm3d_time_POINT12.tec
)

Benchmark(AUTHOR NW
	PATH THM/thm_decov
	CONFIG FEM
	RUNTIME 442
	OUTPUT_FILES
		thm_decov_domain_tri.tec
		thm_decov_ply_H_PROFILE_t1.tec
		thm_decov_ply_V_PROFILE_t2.tec
		thm_decov_time_POINT1.tec
		thm_decov_time_POINT10.tec
		thm_decov_time_POINT11.tec
		thm_decov_time_POINT12.tec
		thm_decov_time_POINT13.tec
		thm_decov_time_POINT15.tec
		thm_decov_time_POINT18.tec
		thm_decov_time_POINT2.tec
		thm_decov_time_POINT6.tec
		thm_decov_time_POINT7.tec
		thm_decov_time_POINT8.tec
		thm_decov_time_POINT9.tec
)

Benchmark(AUTHOR NW
	PATH NUMERICS/SUPG/T_adv_diff_transient_SUPG_line
	CONFIG FEM
	RUNTIME 48
	OUTPUT_FILES
		T_adv_diff_transient_SUPG_line_ply_PLY_0_t1.tec
		T_adv_diff_transient_SUPG_line_time_POINT1.tec
)
