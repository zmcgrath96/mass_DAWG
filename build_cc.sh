# make this script executable with
# $> chmod u+x build.sh
# then run with ./build.sh
cd src
make clean 
make
cd ../tests
make 
cd ..