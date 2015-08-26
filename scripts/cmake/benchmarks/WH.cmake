include(test/Benchmark)

Benchmark(AUTHOR WH
	PATH C/IPhreeqcCoupling/isotope_fractionation/1d_isofrac
	REQUIRED_CMAKE_OPTIONS OGS_FEM_IPQC
	RUNTIME 458
	OUTPUT_FILES 1d_isofrac_domain_line.tec
)
