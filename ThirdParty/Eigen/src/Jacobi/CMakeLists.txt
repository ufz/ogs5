file(GLOB Eigen_Jacobi_SRCS "*.h")

install(FILES
  ${Eigen_Jacobi_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/Jacobi COMPONENT Devel
  )
