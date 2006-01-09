file(GLOB Eigen_SVD_SRCS "*.h")

install(FILES
  ${Eigen_SVD_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/SVD COMPONENT Devel
  )
