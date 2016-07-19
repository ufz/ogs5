include(test/Benchmark)

Benchmark(AUTHOR BG
	PATH ECLIPSE_DUMUX/1phase_radialflow_1phase_transport/1pf_1pt
	CONFIG FEM
	OUTPUT_FILES 1pf_1pt_domain_hex.tec
	RUNTIME 16
)

Benchmark(AUTHOR BG
	PATH ECLIPSE_DUMUX/2phase_flow_2phase_tracertransport/2pf_2pt
	CONFIG FEM
	OUTPUT_FILES
		2pf_2pt_domain_hex.tec
		2pf_2pt_time_POINT16.tec
		2pf_2pt_time_POINT17.tec
	RUNTIME 28
)

Benchmark(AUTHOR BG
	PATH ECLIPSE_DUMUX/2phase_flow_radialmodel/2pf_radialmodel
	CONFIG FEM
	OUTPUT_FILES 2pf_radialmodel_domain_hex.tec
	RUNTIME 6
)

Benchmark(AUTHOR BG
	PATH ECLIPSE_DUMUX/kinetic_CO2phase_generation_E100/CO2phase_gen_E100
	CONFIG FEM
	OUTPUT_FILES
		CO2phase_gen_E100_domain_ele.tec
		CO2phase_gen_E100_domain_hex.tec
	RUNTIME 10
)

Benchmark(AUTHOR BG
	PATH ECLIPSE_DUMUX/kinetic_CO2phase_generation_E300/CO2phase_gen_E300
	CONFIG FEM
	OUTPUT_FILES
		CO2phase_gen_E300_domain_ele.tec
		CO2phase_gen_E300_domain_hex.tec
	RUNTIME 10
)

Benchmark(AUTHOR BG
	PATH C/2d_Cl_transport_Clay/Nuklidtransport
	CONFIG FEM
	OUTPUT_FILES
		Nuklidtransport_domain_tri.tec
		Nuklidtransport_POLYLINE_BC_up_TIM_MASS_TRANSPORT.tec
		Nuklidtransport_POLYLINE_PLY_Source_TIM_MASS_TRANSPORT.tec
	RUNTIME 11
)
