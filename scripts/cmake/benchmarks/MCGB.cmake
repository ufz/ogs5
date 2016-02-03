include(test/Benchmark)

Benchmark(AUTHOR MCGB
	PATH T/2units2faults/2units2faults
	CONFIG FEM
	RUNTIME 324
	OUTPUT_FILES
		2units2faults_time_POINT15.tec
		2units2faults_time_POINT16.tec
		2units2faults_time_POINT17.tec
)

Benchmark(AUTHOR MCGB
	PATH T/2units2faults/FCT/2units2faults_FCT
	CONFIG FEM
	RUNTIME 452
	OUTPUT_FILES
		2units2faults_FCT_time_POINT15.tec
		2units2faults_FCT_time_POINT16.tec
		2units2faults_FCT_time_POINT17.tec
)
