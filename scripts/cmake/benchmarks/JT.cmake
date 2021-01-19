include(test/Benchmark)

Benchmark(AUTHOR JT
	PATH H/HetGWFlow/aniso/het_test
	CONFIG FEM
	RUNTIME 5
	OUTPUT_FILES
            het_test_domain_ele_LIQUID_FLOW_tri.tec
            het_test_domain_LIQUID_FLOW_tri.tec
)