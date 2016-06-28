include(InstallRequiredSystemLibraries)
set(CPACK_GENERATOR ZIP)
if(NOT CMAKE_CROSSCOMPILING)
	set(CPACK_GENERATOR NSIS ZIP)
endif()

file(TO_NATIVE_PATH "${CMAKE_SOURCE_DIR}/scripts/cmake/OGS_Logo_Installer.bmp" BACKGROUND_IMAGE)
set(CPACK_PACKAGE_ICON ${BACKGROUND_IMAGE})
set(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_DESCRIPTION_SUMMARY}")
set(CPACK_NSIS_CONTACT "info@opengeosys.org")
set(CPACK_NSIS_MODIFY_PATH OFF)
set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
set(CPACK_NSIS_MENU_LINKS
	"bin" "Executables folder"
	"http://www.opengeosys.org" "Website"
	"https://github.com/ufz/ogs5" "Source code on GitHub"
)
