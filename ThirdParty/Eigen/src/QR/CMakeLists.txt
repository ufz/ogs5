file(GLOB Eigen_QR_SRCS "*.h")

install(FILES
  ${Eigen_QR_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/QR COMPONENT Devel
  )
