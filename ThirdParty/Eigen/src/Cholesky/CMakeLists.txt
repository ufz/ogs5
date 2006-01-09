file(GLOB Eigen_Cholesky_SRCS "*.h")

install(FILES
  ${Eigen_Cholesky_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/Cholesky COMPONENT Devel
  )
