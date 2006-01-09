file(GLOB Eigen_misc_SRCS "*.h")

install(FILES
  ${Eigen_misc_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/misc COMPONENT Devel
  )
