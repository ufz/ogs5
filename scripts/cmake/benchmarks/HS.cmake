include(test/Benchmark)

Benchmark(AUTHOR HS
	PATH C/1d_analyt/1d_1
	CONFIG FEM
	RUNTIME 21
	OUTPUT_FILES 1d_1_domain_tet.tec
)

Benchmark(AUTHOR HS
	PATH C/2d_analyt/2d_1
	CONFIG FEM
	RUNTIME 22
	OUTPUT_FILES 2d_1_domain_tri.tec
)

# TODO ODE solver
#Benchmark(AUTHOR HS
#	PATH C/1d_xylene_degradation/h2_line
#	CONFIG FEM
#	RUNTIME 9
#	OUTPUT_FILES h2_line_domain_NO_PCS_line.tec
#)

# TODO ODE solver
#Benchmark(AUTHOR HS
#	PATH C/1d_TCEaufEisen/1d_TCE_Ion
#	CONFIG FEM
#	RUNTIME 104
#	OUTPUT_FILES 1d_TCE_Ion_domain_line.tec
#)

Benchmark(AUTHOR HS
	PATH C/decay/HC_decay_1Du
	CONFIG FEM
	RUNTIME 2
	OUTPUT_FILES HC_decay_1Du_domain_line.tec
)

Benchmark(AUTHOR HS
	PATH C/diffusion/Diff_HTO_test
	CONFIG FEM
	RUNTIME 68
	OUTPUT_FILES Diff_HTO_test_domain_line.tec
)

Benchmark(AUTHOR HS
	PATH C/diffusion/diff_aniso
	CONFIG FEM
	RUNTIME 5
	OUTPUT_FILES diff_aniso_domain_tri.tec
)

Benchmark(AUTHOR HS
	PATH C/sorption_decay/HC_decay_sorp_henry_1Du
	CONFIG FEM
	RUNTIME 3
	OUTPUT_FILES HC_decay_sorp_henry_1Du_domain_line.tec
)

Benchmark(AUTHOR HS
	PATH C/sorption/Freundlich/HC_sorp_Freundl_1D
	CONFIG FEM
	RUNTIME 4
	OUTPUT_FILES HC_sorp_Freundl_1D_domain_line.tec
)

Benchmark(AUTHOR HS
	PATH C/sorption/Henry/HC_sorp_henry_1D
	CONFIG FEM
	RUNTIME 3
	OUTPUT_FILES HC_sorp_henry_1D_domain_line.tec
)

Benchmark(AUTHOR HS
	PATH C/sorption/Langmuir/HC_sorp_langmuir_1D
	CONFIG FEM
	RUNTIME 3
	OUTPUT_FILES HC_sorp_langmuir_1D_domain_line.tec
)

Benchmark(AUTHOR HS
	PATH C/calcite_pqc/pds
	CONFIG PQC
	RUNTIME 226
	OUTPUT_FILES pds_domain_line.tec
)

Benchmark(AUTHOR HS
	PATH C/calcite_gems/calcite
	CONFIG GEMS
	RUNTIME 232
	OUTPUT_FILES calcite0210.vtk
)

Benchmark(AUTHOR HS
	PATH C/comedy2d/cement2d
	CONFIG GEMS
	RUNTIME 327
	OUTPUT_FILES cement2d0030.vtk
)

Benchmark(AUTHOR HS
	PATH C/HAYEKIT/ab1d
	CONFIG GEMS
	RUNTIME 230
	OUTPUT_FILES ab1d0100.vtk
)

Benchmark(AUTHOR HS
	PATH PETSc/TransLay2d/lag2d
	CONFIG PETSC_GEMS
	RUNTIME 210
	NUM_PROCESSORS 4
	OUTPUT_FILES lag2d0001.vtk
)

Benchmark(AUTHOR HS
	PATH PETSc/ConcreteCrack/decal
	CONFIG PETSC_GEMS
	RUNTIME 114
	NUM_PROCESSORS 4
	OUTPUT_FILES decal0008.vtk
)
