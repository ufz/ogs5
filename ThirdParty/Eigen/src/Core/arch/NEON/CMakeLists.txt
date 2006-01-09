file(GLOB Eigen_Core_arch_NEON_SRCS "*.h")

install(FILES
  ${Eigen_Core_arch_NEON_SRCS}
  DESTINATION ${INCLUDE_INSTALL_DIR}/Eigen/src/Core/arch/NEON COMPONENT Devel
)
