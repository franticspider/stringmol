#!/bin/bash


#g++ main.cpp -o travis_gcc_cpp98
#./travis_gcc_cpp98
#cppcheck --quiet --error-exitcode=1 main.cpp

echo "Should compile now, but.. it's a TODO!"




echo "============================"
echo "Checking with cppcheck"
cd src
# the --check-config argument seems to find standard library headers
#cppcheck --error-exitcode=1 --force --check-config --enable=all *.cpp
echo "----------------------------"
echo "Checking config with cppcheck"
cppcheck --error-exitcode=1 --force --check-config --suppress=missingIncludeSystem .
echo "----------------------------"
echo "Checking code with cppcheck"
cppcheck --error-exitcode=1 --force --enable=all --inline-suppr .
cd ../
echo ""

# NB: use "// cppcheck-suppress unusedFunction" before functions to suppress warnings
#           //cppcheck-suppress invalidscanf_libc"

echo "============================"
echo "Running Tests.  Please Wait."
#g++ -Wall Shapes-Catch-Testing-Example/Source/Shapes-Catch-Testing-Example.cpp Shapes-Catch-Testing-Example/Source/Implementation/*.cpp -o main
#g++ -Wall Shapes-Catch-Testing-Example/Test/*.cpp Shapes-Catch-Testing-Example/Source/Implementation/*.cpp -o test
cd tests
echo "  compiling..."
g++ -std=gnu++11 -Wall -o test  *.cpp ../release/mt19937-2.o ../release/randutil.o 
cd ..
echo ""
echo "  now testing.."
./tests/test


echo "  cleaning up.."
rm -f rng.txt
