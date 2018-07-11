Benchmark(AUTHOR YS
	PATH RWPT/Harter/colloid_t
	CONFIG FEM
	RUNTIME 2
	OUTPUT_FILES
		colloid_t_LIQUID_FLOW60.vtu
		colloid_t_RWPT_60.particles.vtk
		colloid_t_time_BC_DOWN_RANDOM_WALK.tec
)

Benchmark(AUTHOR YS
	PATH RWPT/2DGrains/2d_grains
	CONFIG FEM
	RUNTIME 39
	OUTPUT_FILES 
		2d_grains_GROUNDWATER_FLOW60.vtu
		2d_grains_RWPT_particles60.vtp
)

Benchmark(AUTHOR YS
	PATH RWPT/3DGrain/3d_grain
	CONFIG FEM
	RUNTIME 43
	OUTPUT_FILES
		3d_grain_GROUNDWATER_FLOW20.vtu
		3d_grain_RWPT_particles20.vtp
)

Benchmark(AUTHOR YS
	PATH RWPT/HomoCube/3DRWPTCubTet
	CONFIG FEM
	RUNTIME 27
	OUTPUT_FILES
	3DRWPTCubTet_GROUNDWATER_FLOW120.vtu
	3DRWPTCubTet_RWPT_particles120.vtp
)

Benchmark(AUTHOR YS
	PATH RWPT/Forchheimer/forchheimer_rwpt
	CONFIG FEM
	RUNTIME 4
	OUTPUT_FILES
		forchheimer_rwpt_GROUNDWATER_FLOW60.vtu
		forchheimer_rwpt_RWPT_particles60.vtp
		forchheimer_rwpt_domain_ele_GROUNDWATER_FLOW.tec
		forchheimer_rwpt_time_BC_DOWN_RANDOM_WALK.tec
)
