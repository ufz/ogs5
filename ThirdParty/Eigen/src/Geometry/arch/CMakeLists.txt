file(GLOB Eigen_Geometry_arch_SRCS "*.h")

install(FILES
  ${Eigen_Geometry_arch_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/Geometry/arch COMPONENT Devel
  )
