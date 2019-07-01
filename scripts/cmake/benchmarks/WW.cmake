include(test/Benchmark)

Benchmark(AUTHOR WW
	PATH H2/BuckleyLeverett/buck
	CONFIG FEM
	NUM_PROCESSORS 1
	RUNTIME 19
	OUTPUT_FILES buck_ply_OUT_t0.tec
)

Benchmark(AUTHOR WW
	PATH H2/Liakopoulos/Line/h2_Liako
	CONFIG FEM
	NUM_PROCESSORS 1
	RUNTIME 1
	OUTPUT_FILES
		h2_Liako_time_POINT0.tec
		h2_Liako_ply_left_t0.tec
)

Benchmark(AUTHOR WW
	PATH H2/Liakopoulos/Quad/h2_Liako
	CONFIG FEM
	NUM_PROCESSORS 1
	RUNTIME 2
	OUTPUT_FILES
	h2_Liako_time_POINT0.tec
	h2_Liako_domain_quad.tec
	h2_Liako_ply_left_t1.tec
)

Benchmark(AUTHOR WW
	PATH H2/Liakopoulos/Tri/h2_Liako
	CONFIG FEM
	NUM_PROCESSORS 1
	RUNTIME 2
	OUTPUT_FILES
	h2_Liako_time_POINT0.tec
	h2_Liako_domain_tri.tec
	h2_Liako_ply_left_t1.tec
)

Benchmark(AUTHOR WW
	PATH H2/McWhorter/mcwt
	CONFIG FEM
	NUM_PROCESSORS 1
	RUNTIME 36
	OUTPUT_FILES
		mcwt_domain_quad.tec
		mcwt_time_POINT4.tec
		mcwt_ply_TOP_t0.tec
)

Benchmark(AUTHOR WW
	PATH H2/McWhorter/1D/mcwt
	CONFIG FEM
	NUM_PROCESSORS 1
	RUNTIME 59
	OUTPUT_FILES
		mcwt_ply_Profile_t0.tec
		mcwt_domain_line.tec
		mcwt_time_POINT0.tec
)

Benchmark(AUTHOR WW
	PATH H2/Liakopoulos/Tet/h2_Liako
	CONFIG FEM
	NUM_PROCESSORS 1
	RUNTIME 19
	OUTPUT_FILES
		h2_Liako_time_POINT0.tec
		h2_Liako_ply_left_t1.tec
		h2_Liako_domain_tet.tec
)

Benchmark(AUTHOR WW
	PATH H2/Liakopoulos/Hex/h2_Liako
	CONFIG FEM
	NUM_PROCESSORS 1
	RUNTIME 19
	OUTPUT_FILES
		h2_Liako_time_POINT0.tec
		h2_Liako_ply_left_t1.tec
		h2_Liako_domain_hex.tec
)

Benchmark(AUTHOR WW
	PATH HM/excavation/BoreholeExcavation/borehole_excav
	CONFIG FEM
	NUM_PROCESSORS 1
	RUNTIME 40
	OUTPUT_FILES
		borehole_excav_ply_horizon_t2.tec
		borehole_excav_ply_vertikal_t3.tec
)


Benchmark(AUTHOR WW
	PATH TH2M/H2M_TEP/w_exp
	CONFIG FEM
	NUM_PROCESSORS 1
	RUNTIME 8
	OUTPUT_FILES
		w_exp_ply_LEFT_t9_DEFORMATION.tec
		w_exp_time_POINT6_DEFORMATION.tec
		w_exp_domain_MULTI_PHASE_FLOW_tri.tec
		w_exp_time_POINT5_MULTI_PHASE_FLOW.tec
		w_exp_time_POINT4_MULTI_PHASE_FLOW.tec
		w_exp_time_POINT5_DEFORMATION.tec
		w_exp_time_POINT4_DEFORMATION.tec
		w_exp_domain_DEFORMATION_tri.tec
		w_exp_time_POINT6_MULTI_PHASE_FLOW.tec
		w_exp_ply_LEFT_t8_MULTI_PHASE_FLOW.tec
)

Benchmark(AUTHOR WW
	PATH TH2M/Newton/th2m_quad
	CONFIG FEM
	NUM_PROCESSORS 1
	RUNTIME 19
	OUTPUT_FILES
		th2m_quad_domain_DEFORMATION_H2_quad.tec
		th2m_quad_time_POINT6_DEFORMATION_H2.tec
		th2m_quad_time_POINT7_DEFORMATION_H2.tec
)

Benchmark(AUTHOR WW
	PATH TH2M/th2m_quad
	CONFIG FEM
	NUM_PROCESSORS 1
	RUNTIME 140
	OUTPUT_FILES
		th2m_quad_domain_DEFORMATION_quad.tec
		th2m_quad_time_POINT6_HEAT_TRANSPORT.tec
		th2m_quad_time_POINT6_DEFORMATION.tec
		th2m_quad_domain_HEAT_TRANSPORT_quad.tec
		th2m_quad_time_POINT6_MULTI_PHASE_FLOW.tec
		th2m_quad_domain_MULTI_PHASE_FLOW_quad.tec
)

Benchmark(AUTHOR WW
	PATH M/creep/ImplicitBGRa/ImplicitBGRa
	CONFIG FEM
	NUM_PROCESSORS 1
	RUNTIME 140
	OUTPUT_FILES
		ImplicitBGRa_time_POINT1.tec
)

Benchmark(AUTHOR WW
	PATH M/creep/CreepAfterExcavationImplicitBGRa/creep_after_excavation
	CONFIG FEM
	NUM_PROCESSORS 1
	RUNTIME 320
	OUTPUT_FILES
		creep_after_excavation_time_point_top.tec
)

Benchmark(AUTHOR WW
	PATH PETSc/h_tri/h_tri
	CONFIG PETSC
	NUM_PROCESSORS 3
	RUNTIME 200
	OUTPUT_FILES
		h_tri_domain_tri_0.tec
		h_tri_domain_tri_1.tec
		h_tri_domain_tri_2.tec
)

Benchmark(AUTHOR WW
	PATH PETSc/McWhorter/mcwt
	CONFIG PETSC
	NUM_PROCESSORS 4
	RUNTIME 178
	OUTPUT_FILES
		mcwt_domain_MULTI_PHASE_FLOW_quad_0.tec
		mcwt_domain_MULTI_PHASE_FLOW_quad_1.tec
		mcwt_domain_MULTI_PHASE_FLOW_quad_2.tec
		mcwt_domain_MULTI_PHASE_FLOW_quad_3.tec
)

Benchmark(AUTHOR WW
	PATH PETSc/Richards/h_us_quad
	CONFIG PETSC
	NUM_PROCESSORS 4
	RUNTIME 217
	OUTPUT_FILES
		h_us_quad_domain_RICHARDS_FLOW_quad_0.tec
		h_us_quad_domain_RICHARDS_FLOW_quad_1.tec
		h_us_quad_domain_RICHARDS_FLOW_quad_2.tec
		h_us_quad_domain_RICHARDS_FLOW_quad_3.tec
)

Benchmark(AUTHOR WW
	PATH PETSc/hm_tri/hm_tri
	CONFIG PETSC
	NUM_PROCESSORS 4
	RUNTIME 200
	OUTPUT_FILES
		hm_tri_domain_DEFORMATION_tri_0.tec
		hm_tri_domain_DEFORMATION_tri_1.tec
		hm_tri_domain_DEFORMATION_tri_2.tec
		hm_tri_domain_DEFORMATION_tri_3.tec
		hm_tri_domain_LIQUID_FLOW_tri_0.tec
		hm_tri_domain_LIQUID_FLOW_tri_1.tec
		hm_tri_domain_LIQUID_FLOW_tri_2.tec
		hm_tri_domain_LIQUID_FLOW_tri_3.tec
)

Benchmark(AUTHOR WW
	PATH PETSc/KueperProblem-PS/kueper
	CONFIG PETSC
	NUM_PROCESSORS 3
	RUNTIME 183
	OUTPUT_FILES kueper_time_POINT4_PS_GLOBAL.tec
)

Benchmark(AUTHOR WW
	PATH PETSc/m_tri/m_tri
	CONFIG PETSC
	NUM_PROCESSORS 2
	RUNTIME 216
	OUTPUT_FILES
		m_tri_domain_DEFORMATION_quad_0.tec
		m_tri_domain_DEFORMATION_quad_1.tec
)

Benchmark(AUTHOR WW
	PATH PETSc/T_tri/t_tri
	CONFIG PETSC
	NUM_PROCESSORS 4
	RUNTIME 181
	OUTPUT_FILES t_tri_time_POINT4.tec
)

Benchmark(AUTHOR WW
	PATH MPI/h_tri/h_tri
	CONFIG MPI
	NUM_PROCESSORS 5
	RUNTIME 210
	OUTPUT_FILES h_tri_domain_tri_0.tec
)

Benchmark(AUTHOR WW
	PATH MPI/HM/hm
	CONFIG MPI
	NUM_PROCESSORS 2
	RUNTIME 190
	OUTPUT_FILES hm_domain_quad_0.tec
)

Benchmark(AUTHOR WW
	PATH MPI/m_tri/m_tri
	CONFIG MPI
	NUM_PROCESSORS 2
	RUNTIME 210
	OUTPUT_FILES m_tri_domain_tri_0.tec
)

Benchmark(AUTHOR WW
	PATH MPI/thm_quad/thm_quad
	CONFIG MPI
	NUM_PROCESSORS 4
	RUNTIME 204
	OUTPUT_FILES thm_quad_domain_quad_0.tec
)

Benchmark(AUTHOR WW
	PATH MPI/McWhorter/mcwt
	CONFIG MPI
	NUM_PROCESSORS 4
	RUNTIME 205
	OUTPUT_FILES
		mcwt_ply_TOP_t0.tec
		mcwt_time_POINT4.tec
)

Benchmark(AUTHOR WW
	PATH MPI/urach/urach
	CONFIG MPI
	NUM_PROCESSORS 4
	RUNTIME 209
	OUTPUT_FILES
		urach_domain_tet_0.tec
		urach_ply_OUT_t1_HEAT_TRANSPORT.tec
)

Benchmark(AUTHOR WW
	PATH MPI/THM/thm_decov
	CONFIG MPI
	NUM_PROCESSORS 3
	RUNTIME 209
	OUTPUT_FILES
		thm_decov_time_POINT1.tec
		thm_decov_time_POINT2.tec
		thm_decov_time_POINT9.tec
		thm_decov_time_POINT18.tec
)
