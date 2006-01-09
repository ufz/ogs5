# Build scripts #

These scripts are intended to be used from a shell (Git Bash on Windows).

## build_gui.sh ##

Builds the GUI configuration.

Parameters:

- `-d PATH` - the path of the build relative to source directory
- `-a ARCH` - the compiler architecture: `x32` or `x64` (Windows only)

A typical call would look like this:

		cd /scripts/build
		./build_gui.sh -d ../build_gui -a x64


## build_configs.sh ##

Builds all FEM configurations.

Parameters:

- `-d PATH` - the path of the build relative to the source directory
- `-a ARCH` - the compiler architecture: `x32` or `x64` (Windows only)

This script will create subdirectories for each configuration. All executables are copied over to a `Release` directory.

When called like this:

		cd /scripts/build
		./build_configs.sh -d ../configs -a x64

You get the following directory structure:

		/sources
		/configs
		    /build_fem
		    /build_sp
		    /build_...
		    /Release
		      ogs.exe
		      ogs_sp.exe
		      ogs_...