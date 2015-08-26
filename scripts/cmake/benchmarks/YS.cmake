Benchmark(AUTHOR YS
	PATH RWPT/Harter/colloid_t
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 2
	OUTPUT_FILES
		colloid_t_LIQUID_FLOW60.vtu
		colloid_t_RWPT_60.particles.vtk
		colloid_t_time_BC_DOWN_RANDOM_WALK.tec
)

Benchmark(AUTHOR YS
	PATH RWPT/2DGrains/2d_grains
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 39
	OUTPUT_FILES 2d_grains_GROUNDWATER_FLOW60.vtu
)

Benchmark(AUTHOR YS
	PATH RWPT/3DGrain/3d_grain
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 43
	OUTPUT_FILES
		3d_grain_GROUNDWATER_FLOW20.vtu
		3d_grain_RWPT_20.particles.vtk
)

Benchmark(AUTHOR YS
	PATH RWPT/HomoCube/3DRWPTCubTet
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 27
	OUTPUT_FILES
	3DRWPTCubTet_GROUNDWATER_FLOW120.vtu
	3DRWPTCubTet_RWPT_120.particles.vtk
)

Benchmark(AUTHOR YS
	PATH RWPT/Forchheimer/forchheimer_rwpt
	REQUIRED_CMAKE_OPTIONS OGS_FEM
	RUNTIME 4
	OUTPUT_FILES
		forchheimer_rwpt_GROUNDWATER_FLOW60.vtu
		forchheimer_rwpt_RWPT_60.particles.vtk
		forchheimer_rwpt_domain_ele_GROUNDWATER_FLOW.tec
		forchheimer_rwpt_time_BC_DOWN_RANDOM_WALK.tec
)
