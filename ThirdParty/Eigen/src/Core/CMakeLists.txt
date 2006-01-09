file(GLOB Eigen_Core_SRCS "*.h")

install(FILES
  ${Eigen_Core_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/Core COMPONENT Devel
  )

add_subdirectory(products)
add_subdirectory(util)
add_subdirectory(arch)
