cd styleguide
for dir in "Base" "FEM" "FileIO" "GEO" "MathLib" "MSH" "MSHGEOTOOLS" "OGS" "OGSProject" "Qt" "UTL"
do
	./uncrustify.sh ../../"$dir" h
	./uncrustify.sh ../../"$dir" cpp
done
cd ..