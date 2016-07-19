include(test/Benchmark)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/1D_isofrac/1d_isofrac_AS
#	CONFIG FEM
#	RUNTIME 4
#	OUTPUT_FILES 1d_isofrac_AS_ply_OUT_LINE_t0.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/1D_isofrac/1d_isofrac
#	CONFIG FEM
#	RUNTIME 6
#	OUTPUT_FILES 1d_isofrac_ply_OUT_LINE_t0.tec
#)

Benchmark(AUTHOR CB
	PATH GROUNDWATER_FLOW/Transient_flow/trans_bd_homo
	CONFIG FEM
	RUNTIME 4
	OUTPUT_FILES trans_bd_homo_domain_GROUNDWATER_FLOW_quad.tec
)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/FG_3ports/rt1
#	CONFIG FEM
#	RUNTIME 34
#	OUTPUT_FILES rt1_domain_quad.tec
#)

Benchmark(AUTHOR CB
	PATH C/hetk+n+restart/2D1P_transport
	CONFIG FEM
	RUNTIME 47
	OUTPUT_FILES
		2D1P_transport_domain_quad.tec
		2D1P_transport_time_POINT10.tec
		2D1P_transport_time_POINT11.tec
)

Benchmark(AUTHOR CB
	PATH C/ChemAppCoupling/cation_exchange/cmp9
	CONFIG FEM
	RUNTIME 91
	OUTPUT_FILES cmp9_time_POINT2.tec
)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/ChemAppCoupling/slow_kin_cap_act_precalc/cmp8
#	CONFIG FEM
#	RUNTIME 134
#	OUTPUT_FILES
#		cmp8_Activities.dump
#		cmp8_EquilibriumConstants.dump
#		cmp8_ply_OUT_LINE_t1.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/ChemAppCoupling/Wagrien_BatchEqui_PVLE_CAP/wagrien_1D
#	CONFIG FEM
#	RUNTIME 169
#	OUTPUT_FILES
#		wagrien_1D_domain_ele.tec
#		wagrien_1D_ply_OUT_LINE_t0.tec
#)

# require phreeqc executable
# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/Engesgaard/2Kin/slow_kin_pqc/pds
#	CONFIG FEM
#	RUNTIME 2
#	OUTPUT_FILES pds_ply_OUT_LINE_t1.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/Engesgaard/2Kin/slow_kin_pqc_krc/pds
#	CONFIG FEM
#	RUNTIME 2
#	OUTPUT_FILES
#		pds_ply_OUT_LINE_t1.tec
#		pds_Activities.dump
#		pds_EquilibriumConstants.dump
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/Engesgaard/Additional_pqc_output/pds
#	CONFIG FEM
#	RUNTIME 2
#	OUTPUT_FILES pds_ply_OUT_LINE_t1.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/Engesgaard/equi/calcite_pqc_volume/pds
#	CONFIG FEM
#	RUNTIME 2
#	OUTPUT_FILES pds_ply_OUT_LINE_t1.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/Engesgaard/equi/calcite_pqc_weight/pds
#	CONFIG FEM
#	RUNTIME 2
#	OUTPUT_FILES pds_ply_OUT_LINE_t1.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/Engesgaard/Kin/fast_kin_pqc/pds
#	CONFIG FEM
#	RUNTIME 2
#	OUTPUT_FILES pds_ply_OUT_LINE_t1.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/Engesgaard/Kin/fast_kin_pqc_krc/pds
#	CONFIG FEM
#	RUNTIME 453
#	OUTPUT_FILES
#		pds_ply_OUT_LINE_t1.tec
#		pds_Activities.dump
#		pds_EquilibriumConstants.dump
#)

Benchmark(AUTHOR CB
	PATH C/Engesgaard/Kin/slow_kin_pqc/pds
	CONFIG FEM
	RUNTIME 2
	OUTPUT_FILES pds_ply_OUT_LINE_t1.tec
)

Benchmark(AUTHOR CB
	PATH C/Engesgaard/Kin/slow_kin_pqc_krc/pds
	CONFIG FEM
	RUNTIME 430
	OUTPUT_FILES
		pds_ply_OUT_LINE_t1.tec
		pds_Activities.dump
		pds_EquilibriumConstants.dump
)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/NAPL-dissolution/1D_NAPL-diss_dissolution/1D_TPF_resS_trans
#	CONFIG FEM
#	RUNTIME 94
#	OUTPUT_FILES 1D_TPF_resS_trans_time_POINT17.tec
#)

Benchmark(AUTHOR CB
	PATH C/NAPL-dissolution/1D_NAPL-diss_flow/hex/1D_TPF_resS_flow
	CONFIG FEM
	RUNTIME 2
	OUTPUT_FILES 1D_TPF_resS_flow_ply_OUT_LINE_t0.tec
)

Benchmark(AUTHOR CB
	PATH C/NAPL-dissolution/1D_NAPL-diss_flow/lin/1D_TPF_resS_flow
	CONFIG FEM
	RUNTIME 2
	OUTPUT_FILES 1D_TPF_resS_flow_ply_OUT_LINE_t0.tec
)

Benchmark(AUTHOR CB
	PATH C/NAPL-dissolution/1D_NAPL-diss_flow/pri/1D_TPF_resS_flow
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES 1D_TPF_resS_flow_ply_OUT_LINE_t0.tec
)

Benchmark(AUTHOR CB
	PATH C/NAPL-dissolution/1D_NAPL-diss_flow/quad/1D_TPF_resS_flow
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES 1D_TPF_resS_flow_ply_OUT_LINE_t0.tec
)

Benchmark(AUTHOR CB
	PATH C/NAPL-dissolution/1D_NAPL-diss_flow/tet/1D_TPF_resS_flow
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES 1D_TPF_resS_flow_ply_OUT_LINE_t0.tec
)

Benchmark(AUTHOR CB
	PATH C/NAPL-dissolution/1D_NAPL-diss_flow/tri/1D_TPF_resS_flow
	CONFIG FEM
	RUNTIME 1
	OUTPUT_FILES 1D_TPF_resS_flow_ply_OUT_LINE_t0.tec
)

Benchmark(AUTHOR CB
	PATH C/NAPL-dissolution/1D_NAPL-diss_trans/hex/1D_TPF_resS_trans
	CONFIG FEM
	RUNTIME 58
	OUTPUT_FILES
		1D_TPF_resS_trans_ply_OUT_LINE_t1.tec
		1D_TPF_resS_trans_time_POINT8.tec
)

Benchmark(AUTHOR CB
	PATH C/NAPL-dissolution/1D_NAPL-diss_trans/lin/1D_TPF_resS_trans
	CONFIG FEM
	RUNTIME 4
	OUTPUT_FILES
		1D_TPF_resS_trans_ply_OUT_LINE_t1.tec
		1D_TPF_resS_trans_time_POINT8.tec
)

Benchmark(AUTHOR CB
	PATH C/NAPL-dissolution/1D_NAPL-diss_trans/pri/1D_TPF_resS_trans
	CONFIG FEM
	RUNTIME 24
	OUTPUT_FILES
		1D_TPF_resS_trans_ply_OUT_LINE_t1.tec
		1D_TPF_resS_trans_time_POINT8.tec
)

Benchmark(AUTHOR CB
	PATH C/NAPL-dissolution/1D_NAPL-diss_trans/quad/1D_TPF_resS_trans
	CONFIG FEM
	RUNTIME 8
	OUTPUT_FILES
		1D_TPF_resS_trans_ply_OUT_LINE_t1.tec
		1D_TPF_resS_trans_time_POINT8.tec
)

Benchmark(AUTHOR CB
	PATH C/NAPL-dissolution/1D_NAPL-diss_trans/tet/1D_TPF_resS_trans
	CONFIG FEM
	RUNTIME 31
	OUTPUT_FILES
		1D_TPF_resS_trans_ply_OUT_LINE_t1.tec
		1D_TPF_resS_trans_time_POINT8.tec
)

Benchmark(AUTHOR CB
	PATH C/NAPL-dissolution/1D_NAPL-diss_trans/tri/1D_TPF_resS_trans
	CONFIG FEM
	RUNTIME 6
	OUTPUT_FILES
		1D_TPF_resS_trans_ply_OUT_LINE_t1.tec
		1D_TPF_resS_trans_time_POINT8.tec
)

# require phreeqc executable
# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/Poro-Perm_Update/Lagneau_Batch/pds
#	CONFIG FEM
#	RUNTIME 13
#	OUTPUT_FILES
#		pds_node_porosities.dump
#		pds_ply_OUT_LINE_t0.tec
#		pds_time_POINT2.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/ReactDeact/1D_isofrac_deac1/1d_isofrac
#	CONFIG FEM
#	RUNTIME 7
#	OUTPUT_FILES
#		1d_isofrac_ply_OUT_LINE_t0.tec
#		1d_isofrac_Deactivated_nodes.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/ReactDeact/1D_isofrac_deac1/1d_isofrac_AS
#	CONFIG FEM
#	RUNTIME 4
#	OUTPUT_FILES
#		1d_isofrac_AS_ply_OUT_LINE_t0.tec
#		1d_isofrac_AS_Deactivated_nodes.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/ReactDeact/1D_isofrac_deac2/1d_isofrac
#	CONFIG FEM
#	RUNTIME 7
#	OUTPUT_FILES
#		1d_isofrac_ply_OUT_LINE_t0.tec
#		1d_isofrac_Deactivated_nodes.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/ReactDeact/1D_isofrac_deac2/1d_isofrac_AS
#	CONFIG FEM
#	RUNTIME 4
#	OUTPUT_FILES
#		1d_isofrac_AS_ply_OUT_LINE_t0.tec
#		1d_isofrac_AS_Deactivated_nodes.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/ReactDeact/1D_isofrac_deac3/1d_isofrac
#	CONFIG FEM
#	RUNTIME 26
#	OUTPUT_FILES
#		1d_isofrac_ply_OUT_LINE_t0.tec
#		1d_isofrac_Deactivated_nodes.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/ReactDeact/1D_isofrac_deac3/1d_isofrac_AS
#	CONFIG FEM
#	RUNTIME 5
#	OUTPUT_FILES
#		1d_isofrac_AS_ply_OUT_LINE_t0.tec
#		1d_isofrac_AS_Deactivated_nodes.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH C/ReactDeact/FG_3ports_deac_3/rt1
#	CONFIG FEM
#	RUNTIME 41
#	OUTPUT_FILES
#		rt1_domain_quad.tec
#		rt1_Deactivated_nodes.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH T/CopyValueBHE/CopyValue
#	CONFIG FEM
#	RUNTIME 21
#	OUTPUT_FILES
#		CopyValue_domain_quad.tec
#		CopyValue_time_P12.tec
#		CopyValue_time_P13.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH TCR/Temperature_BacGrowth/bact_growth_new
#	CONFIG FEM
#	RUNTIME 44
#	OUTPUT_FILES bact_growth_new_domain_line.tec
#)

Benchmark(AUTHOR CB
	PATH TCR/Temperature_Diff/HBr_10C_Diff_new
	CONFIG FEM
	RUNTIME 36
	OUTPUT_FILES HBr_10C_Diff_new_ply_PLY_Model_Area_t0.tec
)

Benchmark(AUTHOR CB
	PATH TCR/Temperature_Disp/HBr_Disp
	CONFIG FEM
	RUNTIME 15
	OUTPUT_FILES HBr_Disp_ply_PLY_BC_RIGHT_t1.tec
)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH TCR/Temperature_NAPLdiss_Csat/TCE_10C_new
#	CONFIG FEM
#	RUNTIME 26
#	OUTPUT_FILES
#		TCE_10C_new_ply_Modellgebiet_t0.tec
#		TCE_10C_new_time_POINT4.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR CB
#	PATH TCR/Temperature_NAPLdiss_PhaseDiff/TCE_10C_new
#	CONFIG FEM
#	RUNTIME 23
#	OUTPUT_FILES
#		TCE_10C_new_ply_Modellgebiet_t0.tec
#		TCE_10C_new_time_POINT6.tec
#)
