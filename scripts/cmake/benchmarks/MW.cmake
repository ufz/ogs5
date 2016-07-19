include(test/Benchmark)

Benchmark(AUTHOR MW
	PATH H_us/Darcy/unconf_WO_rch/unconf
	CONFIG FEM
	RUNTIME 16
	OUTPUT_FILES
		unconf_2.vtu
		unconf_ply_bottom_t1.tec
)

Benchmark(AUTHOR MW
	PATH H_us/Darcy/unconf_W_rch/unconf
	CONFIG FEM
	RUNTIME 11
	OUTPUT_FILES
		unconf_2.vtu
		unconf_ply_bottom_t1.tec
)

# Deactivated because results are not confirmed
#Benchmark(AUTHOR MW
#	PATH DENSITY-DEPENDENT_FLOW/goswami-clement/constrBC_PressAsHead_tri/HM
#	CONFIG FEM
#	RUNTIME 361
#	OUTPUT_FILES
#		HM_47.vtu
#		HM_122.vtu
#		HM_161.vtu
#)

# Deactivated because results are not confirmed
#Benchmark(AUTHOR MW
#	PATH DENSITY-DEPENDENT_FLOW/goswami-clement/wholeBC_PressAsHead_tri/HM
#	CONFIG FEM
#	RUNTIME 373
#	OUTPUT_FILES
#		HM_42.vtu
#		HM_119.vtu
#		HM_159.vtu
#)
