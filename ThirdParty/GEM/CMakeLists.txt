# Source files

set( HEADERS
	databr.h
	datach.h
	gdatastream.h
	io_arrays.h
	jama_cholesky.h
	jama_lu.h
	m_const.h
	m_param.h
	ms_multi.h
	node.h
	num_methods.h
        s_kinmet.h
	s_solmod.h
        s_sorpmod.h
	tnt.h
	tnt_array1d.h
	tnt_array2d.h
	tnt_i_refvec.h
	v_user.h
	verror.h
)

set ( SOURCES
	gdatastream.cpp
	io_arrays.cpp
	ipm_chemical2.cpp
	ipm_chemical3.cpp
	ipm_chemical.cpp
	ipm_main.cpp
	ipm_simplex.cpp
	ms_multi_file.cpp
	ms_multi_format.cpp
	ms_param.cpp
	node.cpp
	node_format.cpp
	num_methods.cpp
	s_kinmet.cpp
	s_solmod2.cpp
	s_solmod3.cpp
	s_solmod4.cpp
	s_solmod5.cpp
	s_solmod6.cpp
	s_solmod.cpp
	s_sorpmod.cpp
)

# Create the library
add_library( GEM STATIC 
	     ${SOURCES}
             ${HEADERS} 
)
include_directories(
        ${CMAKE_SOURCE_DIR}/GEM
)


##set_target_properties(GEM PROPERTIES LINKER_LANGUAGE CXX)

if (WIN32)
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPMGEMPLUGIN")
endif ()

if (UNIX)
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPMGEMPLUGIN -D__unix")
endif ()
