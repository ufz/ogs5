file(GLOB Eigen_Geometry_SRCS "*.h")

install(FILES
  ${Eigen_Geometry_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/Geometry COMPONENT Devel
  )

add_subdirectory(arch)
