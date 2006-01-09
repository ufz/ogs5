file(GLOB Eigen_StlSupport_SRCS "*.h")

install(FILES
  ${Eigen_StlSupport_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/StlSupport COMPONENT Devel
  )
