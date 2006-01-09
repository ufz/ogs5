:: Setup Visual Studio environment
if %1!==! goto msvc2005

if "%1" == "msvc2005" (
  :msvc2005
  call "%VS80COMNTOOLS%\..\..\VC\bin\vcvars32.bat"
  set generator="Visual Studio 8 2005"
)

if "%1" == "msvc2008" (
  call "%VS90COMNTOOLS%\..\..\VC\bin\vcvars32.bat"
  set generator="Visual Studio 9 2008"
)

echo Building with %generator%