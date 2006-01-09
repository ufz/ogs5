file(GLOB Eigen_Eigen2Support_SRCS "*.h")

install(FILES
  ${Eigen_Eigen2Support_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/Eigen2Support COMPONENT Devel
  )

add_subdirectory(Geometry)