file(GLOB Eigen_Core_arch_SSE_SRCS "*.h")

install(FILES
  ${Eigen_Core_arch_SSE_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/Core/arch/SSE COMPONENT Devel
)
