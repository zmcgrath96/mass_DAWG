# make this script executable with
# $> chmod u+x build.sh
# then run with ./build.sh
echo "Compiling the c++ source code..."
cd src
make clean 
make
cd ../tests
make 
cd ..
echo "Done"