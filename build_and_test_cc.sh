# make this script executable with
# $> chmod u+x build_and_test_cc.sh
# then run with ./build_and_test_cc.sh
chmod u+x build_cc.sh 
./build_cc.sh 
echo "Running c++ tests..."
./tests/testmain
