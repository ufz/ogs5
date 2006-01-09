# Source files

set( SOURCES
	basic.c
	basicsubs.c
	cl1.c
	cvdense.c
	cvode.c
	dense.c
     	dw.c
	input.c
	integrate.c
	isotopes.c
	kinetics.c
	mainsubs.c
	model.c
	nvector.c
	nvector_serial.c
	output.c
	p2clib.c
	parse.c
	phqalloc.c
	phreeqc_files.c
	pitzer.c
	pitzer_structures.c
	pqcint.c
	prep.c
	print.c
	read.c
	smalldense.c
	spread.c
	step.c
	structures.c
	sundialsmath.c
	tally.c
	tidy.c
	utilities.c
	advection.c
	transport.c
	inverse.c
	readtr.c
)

# Header files
set( HEADERS
     cvdense.h
     cvode.h
     dense.h
     global.h
     input.h
     kinetics.h	
     nvector.h     
     nvector_serial.h
     output.h
     p2c.h
     phqalloc.h
     phreeqc.h
     phrqproto.h
     phrqtype.h
     pitzer.h
     pqcint.h
     smalldense.h
     step.h	
     sundialsmath.h
     sundialstypes.h
)

include_directories (
	.
)


# Create the library
add_library( PQC STATIC
	      ${SOURCES}
	      ${HEADERS}
 )

