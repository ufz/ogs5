README_lib.txt

Potsdam, Germany,  April 2009

This package ("libphreeqc") is based on PHREEQC version 2.15.0-2697,
modified in order to be compiled as static or dynamically linked
library and offering an interface for coupling with other programs.

Please read PHREEQC's User Rights Notice (NOTICE.TXT), which is part
of the original PHREEQC source code and also attached to this
distribution. Note that the software included in the present
distribution is Free Software and the redistribution/modification is
subjected to the very same conditions stated in the PHREEQC's User
Rights Notice. 

1. MODIFICATIONS FROM ORIGINAL PHREEQC

The major modification is the possibility of suppressing file-based
input/output, which is then replaced by a string buffer (the input)
and a normal array of double for the output (currently only the PUNCH
directive is implemented). 

I tried to implement the interface with as few interventions on the
original code as possible. A summary of differences with respect to
the original distribution:
1) new files 
   - phreeqc.h: header to be included from the extern calling function.
   - pqcint.c : contains the source code of the interface function.
   - pqcint.h : header needed for a couple of new functions who need
                particular visibility.
2) files with important modifications
  - Makefile
  - input.c
  - mainsubs.c
  - phreeqc_files.c
3) files with minor modifications
  - globals.h
  - basic.c
  - print.c
4) deleted files
  - main.c
5) coding conventions
  - Every difference from the original source code of PHREEQC has been
    pointed out directly in the source code with comments starting
    with the '// MDL:" tag.
  - The names of all new functions and any new defined global variable
    start with 'P' or 'P_'.

Please note that, this version of libphreeqc EXPECTS NO COMMENT ('#') NOR
EMPTY LINES IN THE INPUT.

2. COMPILING

The code is ansi-c, so virtually every compiler should be able to
compile it. The provided Makefile was successfully tested under both
linux and windows systems - for the latter is needed at least a subset
of the GNU tools to be installed on the system. Compiling under Visual
Studio C++ is still ongoing and a vproj file will be provided in
further releases.

The compilation using "gcc" and "make" just needs the command "make"
to be issued. For the DOS version it will be "make dos". The compiled
library will then be "libphreeqc.a" for linux and "libphreeqc.lib" for
windows. 

I added in the source code the preprocessor directive MDL_DEBUG in
order to make the library more verbose about what it's being done. The
target "deb" defined in the Makefile adds debugging flags to the
compiled code and also activates such additional output.

2.1 COMPILING UNDER WINDOWS WITH Rtools

A sufficient environment (and to my knowledge, the smallest) for the
compilation of libphreeqc with GNU tools under windows is provided by
the rtools package:

http://www.murdoch-sutherland.com/Rtools/Rtools29.exe

Once this package installed (supposed under C:\Rtools), just open a
DOS window, then set the PATH environment variable:

C:\PATH_TO_PQC>cd PATH\TO\PQC
PATH=C:\Rtools\bin;C:\Rtools\perl\bin;C:\Rtools\MinGW;%PATH% 
C:\PATH\TO\PQC>make clean
C:\PATH\TO\PQC>make dos

The last command will compile the static library "libphreeqc.lib"
specifying the "dos" target. Once this done, one can proceed to
compiling GeoSys the usual way.

NOTE: the pre-compiled library is provided as binary in this release,
in order to facilitate the compilation of the new interface in GeoSys
(rf4_pqc.vcproj)



Marco De Lucia 
delucia@gfz-potsdam.de



Last modified: $Id: README_lib.txt 15 2009-04-23 12:01:11Z delucia $
