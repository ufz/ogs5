file(GLOB Eigen_Sparse_SRCS "*.h")

install(FILES 
  ${Eigen_Sparse_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/Sparse COMPONENT Devel
  )
