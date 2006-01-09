# This is Googletest 1.6 compressed to 2 source files

if (HAVE_PTHREADS)
  add_definitions (-DGTEST_HAS_PTHREAD=1)
else ()
  add_definitions (-DGTEST_HAS_PTHREAD=0)
endif ()

add_library (GTest gtest.h gtest-all.cc)

set_target_properties(GTest PROPERTIES FOLDER Testing)