file(GLOB Eigen_LU_SRCS "*.h")

install(FILES 
  ${Eigen_LU_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/LU COMPONENT Devel
  )

add_subdirectory(arch)
