include(InstallRequiredSystemLibraries)
set(CPACK_GENERATOR ZIP)

return() # NSIS disabled for now
set(CPACK_GENERATOR NSIS ZIP)
# There is a bug in NSI that does not handle full unix paths properly. Make
# sure there is at least one set of four (4) backlasshes.
set(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/scripts/cmake\\\\OGS_Logo_Installer.bmp")
#set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\MyExecutable.exe")
set(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_DESCRIPTION_SUMMARY}")
set(CPACK_NSIS_HELP_LINK "https:\\\\\\\\geosys.ufz.de")
set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\geosys.ufz.de")
set(CPACK_NSIS_CONTACT "lars.bilke@ufz.de")
set(CPACK_NSIS_MODIFY_PATH ON)
#set(CPACK_NSIS_MENU_LINKS "https://geosys.ufz.de/trac" "OGS Wiki")
message (STATUS "Packaging set to NSIS")
