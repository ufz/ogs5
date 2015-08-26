include(test/Benchmark)

Benchmark(AUTHOR UJG
	PATH HM/hm_tri
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	 2
	OUTPUT_FILES
		hm_tri_domain_DEFORMATION_FLOW_tri.tec
		hm_tri_ply_RIGHT_t1_DEFORMATION_FLOW.tec
		hm_tri_time_POINT0_DEFORMATION_FLOW.tec
)

Benchmark(AUTHOR UJG
	PATH HM/hm_foot_tri
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 3
	OUTPUT_FILES
		hm_foot_tri_domain_DEFORMATION_FLOW_tri.tec
		hm_foot_tri_ply_ex2_Out_axis_t1_DEFORMATION_FLOW.tec
		hm_foot_tri_ply_ex2_out_Bedge_t2_DEFORMATION_FLOW.tec
)

Benchmark(AUTHOR UJG
	PATH HM/hm_cc_tri_s
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 2
	OUTPUT_FILES
		HM/hm_cc_tri_s_domain_DEFORMATION_FLOW_tri.tec
		HM/hm_cc_tri_s_time_POINT2_DEFORMATION_FLOW.tec
)

Benchmark(AUTHOR UJG
	PATH HM/hm_dyn_tri
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 4
	OUTPUT_FILES
		hm_dyn_tri_domain_DEFORMATION_FLOW_tri.tec
		hm_dyn_tri_ply_ex2_Out_axis_t1_DEFORMATION_FLOW.tec
		hm_dyn_tri_ply_ex2_out_Bedge_t2_DEFORMATION_FLOW.tec
)

Benchmark(AUTHOR UJG
	PATH HM/hm_unsat
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 4
	OUTPUT_FILES
		hm_unsat_domain_quad.tec
		hm_unsat_ply_RIGHT_t1.tec
		hm_unsat_time_POINT0.tec
)

Benchmark(AUTHOR UJG
	PATH M/m_cc_tri_s
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 127
	OUTPUT_FILES
		m_cc_tri_s_domain_tri.tec
		m_cc_tri_s_time_POINT2.tec
)


Benchmark(AUTHOR UJG
	PATH M/m_cc_quad_s
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 88
	OUTPUT_FILES
		m_cc_quad_s_domain_quad.tec
		m_cc_quad_s_time_POINT2.tec
)

Benchmark(AUTHOR UJG
	PATH M/m_ssy_quad
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 2
	OUTPUT_FILES
		m_ssy_quad_domain_quad.tec
		m_ssy_quad_ply_WMCORE_center_t1.tec
		m_ssy_quad_time_POINT3.tec
)

Benchmark(AUTHOR UJG
	PATH M/m_brick_l
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 3
	OUTPUT_FILES
		m_brick_l_domain_DEFORMATION_tet.tec
		m_brick_l_time_POINT9_DEFORMATION.tec
)

Benchmark(AUTHOR UJG
	PATH M/m_sdc
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 29
	OUTPUT_FILES
		m_sdc_domain_tri.tec
		m_sdc_time_POINT3.tec
		m_sdc_time_POINT6.tec
)

Benchmark(AUTHOR UJG
	PATH M/m_excav
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 5
	OUTPUT_FILES
		m_excav_domain_tri.tec
		m_excav_ply_H_PROFILE_t1.tec
		m_excav_ply_LEFT_t2.tec
		m_excav_time_POINT1.tec
		m_excav_time_POINT2.tec
		m_excav_time_POINT7.tec
		m_excav_time_POINT8.tec
		m_excav_time_POINT9.tec
		m_excav_time_POINT10.tec
		m_excav_time_POINT11.tec
		m_excav_time_POINT12.tec
		m_excav_time_POINT13.tec
		m_excav_time_POINT15.tec
)

Benchmark(AUTHOR UJG
	PATH M/m_drift_init
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 1
	OUTPUT_FILES
		m_drift_init_domain_tri.tec
		m_drift_init_ply_PLY_11_t1.tec
		m_drift_init_ply_PLY_13_t2.tec
		m_drift_init_ply_PLY_14_t3.tec
		m_drift_init_ply_PLY_15_t4.tec
)

Benchmark(AUTHOR UJG
	PATH M/m_mises
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 6
	OUTPUT_FILES
		m_mises_domain_tri.tec
		m_mises_time_POINT0.tec
		m_mises_time_POINT5.tec
)

Benchmark(AUTHOR UJG
	PATH M/m_e_transiso_2D
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME1
	OUTPUT_FILES
		m_e_transiso_2D_domain_quad.tec
		m_e_transiso_2D_time_POINT0.tec
		m_e_transiso_2D_time_POINT1.tec
		m_e_transiso_2D_time_POINT2.tec
		m_e_transiso_2D_time_POINT3.tec
)

Benchmark(AUTHOR UJG
	PATH M/creep/m_crp_tri
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 3
	OUTPUT_FILES
		creep/m_crp_tri_domain_tri.tec
		creep/m_crp_tri_time_POINT2.tec
)

Benchmark(AUTHOR UJG
	PATH M/elasticity/m_e_transiso_3D
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 2
	OUTPUT_FILES
		m_e_transiso_3D_domain_hex.tec
		m_e_transiso_3D_time_POINT1.tec
		m_e_transiso_3D_time_POINT2.tec
		m_e_transiso_3D_time_POINT3.tec
		m_e_transiso_3D_time_POINT5.tec
		m_e_transiso_3D_time_POINT6.tec
		m_e_transiso_3D_time_POINT7.tec
)

Benchmark(AUTHOR UJG
	PATH HM/hm_foot_tet
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 63
	OUTPUT_FILES 32
		hm_foot_tet_domain_DEFORMATION_FLOW_tet.tec
		hm_foot_tet_ply_OUT_AXIS_t1_DEFORMATION_FLOW.tec
)

Benchmark(AUTHOR UJG
	PATH M/m_drift
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 72
	OUTPUT_FILES
		m_drift_domain_quad.tec
		m_drift_ply_PLY_11_t1_DEFORMATION.tec
		m_drift_ply_PLY_13_t2_DEFORMATION.tec
		m_drift_ply_PLY_14_t3_DEFORMATION.tec
		m_drift_ply_PLY_15_t4_DEFORMATION.tec
)

Benchmark(AUTHOR UJG
	PATH M/m_triax_lubby1
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 72
	OUTPUT_FILES
		m_triax_lubby1_domain_DEFORMATION_tri.tec
		m_triax_lubby1_time_POINT0_DEFORMATION.tec
		m_triax_lubby1_time_POINT1_DEFORMATION.tec
		m_triax_lubby1_time_POINT2_DEFORMATION.tec
		m_triax_lubby1_time_POINT3_DEFORMATION.tec
)

Benchmark(AUTHOR UJG
	PATH M/creep/m_crp_bgr
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 800
	OUTPUT_FILES
		m_crp_bgr_domain_tri.tec
		m_crp_bgr_time_POINT1.tec
)

Benchmark(AUTHOR UJG
	PATH M/creep/uc_creep01
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 305
	OUTPUT_FILES
		uc_creep01_domain_tri.tec
		uc_creep01_time_POINT1.tec
)

Benchmark(AUTHOR UJG
	PATH M/elasticity/M_e_displacement_3Du
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 55
	OUTPUT_FILES M_e_displacement_3Du_domain_hex.tec)

Benchmark(AUTHOR UJG
	PATH M/elasticity/M_e_stress_3Du
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 55
	OUTPUT_FILES M_e_stress_3Du_domain_hex.tec)

Benchmark(AUTHOR UJG
	PATH M/creep/m_triax_lubby2
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 174
	OUTPUT_FILES
		m_triax_lubby2_domain_DEFORMATION_tri.tec
		m_triax_lubby2_time_POINT0_DEFORMATION.tec
		m_triax_lubby2_time_POINT1_DEFORMATION.tec
		m_triax_lubby2_time_POINT2_DEFORMATION.tec
		m_triax_lubby2_time_POINT3_DEFORMATION.tec
)
