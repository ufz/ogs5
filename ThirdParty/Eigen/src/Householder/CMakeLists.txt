file(GLOB Eigen_Householder_SRCS "*.h")

install(FILES
  ${Eigen_Householder_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/Householder COMPONENT Devel
  )
