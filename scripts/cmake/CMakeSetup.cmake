
# Set additional CMake modules path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/scripts/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/scripts/cmake/cmake")

if(CMAKE_VERSION VERSION_GREATER 3.1)
	cmake_policy(SET CMP0054 NEW)
endif()

# Suppress warning on setting policies
cmake_policy(SET CMP0011 OLD)

# Load addional modules
include(UseBackportedModules)
include(OptionRequires)
include(CppcheckTargets)

# Adds useful macros and variables
include( scripts/cmake/Macros.cmake )

# Provide a way for Visual Studio Express users to turn OFF the new FOLDER
# organization feature. Default to ON for non-Express users. Express users must
# explicitly turn off this option to build CMake in the Express IDE...
option(CMAKE_USE_FOLDERS "Enable folder grouping of projects in IDEs." ON)
mark_as_advanced(CMAKE_USE_FOLDERS)
set_property(GLOBAL PROPERTY USE_FOLDERS ${CMAKE_USE_FOLDERS})