:: Goto sources directory
cd ..

:: Cleanup
rd /S /Q cpp
del /S /Q cpd.xml

:: Documentation
doxygen.exe

:: Code duplication
java -cp c:\pmd-4.2.5\lib\pmd-4.2.5.jar net.sourceforge.pmd.cpd.CPD --minimum-tokens 100 --language cpp --format net.sourceforge.pmd.cpd.XMLRenderer --files Base --files GEO --files FEM --files FileIO --files GEO --files MathLib --files MSH --files OGS --files OGSProject --files PQC --files Qt\Base --files Qt\DataView --files Qt\Gui --files Qt\StationView --files Qt\VtkVis >> cpd.xml

:: CPP Check
mkdir cpp
cd cpp
cppcheck --xml --enable=all ../Base 2> Base.xml
::cppcheck --xml --enable=all ../FEM 2> FEM.xml
cppcheck --xml --enable=all ../FileIO 2> FileIO.xml
cppcheck --xml --enable=all ../GEO 2> GEO.xml
cppcheck --xml --enable=all ../MathLib 2> MathLib.xml
cppcheck --xml --enable=all ../MSH 2> MSH.xml
cppcheck --xml --enable=all ../OGS 2> OGS.xml
cppcheck --xml --enable=all ../OGSProject 2> OGSProject.xml
cppcheck --xml --enable=all ../Qt/Base 2> QtBase.xml
cppcheck --xml --enable=all ../Qt/DataView 2> QtDataView.xml
cppcheck --xml --enable=all ../Qt/Gui 2> QtGui.xml
cppcheck --xml --enable=all ../Qt/StationView 2> QtStationView.xml
cppcheck --xml --enable=all ../Qt/VtkVis 2> QtVtkVis.xml
cppcheck --xml --enable=all ../Qt/VtkAct 2> QtVtkAct.xml