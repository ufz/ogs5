file(GLOB Eigen_plugins_SRCS "*.h")

install(FILES
  ${Eigen_plugins_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/plugins COMPONENT Devel
  )
